"""Each packaged target's side of its consumer's contract, checked without the consumer.

Whether the consumers really read these files is checked in the consumer image; here the
written files must carry what each consumer's reader looks for.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import polars as pl
import pytest
from apb2.api import (
    AnnotationCompiler,
    ParsedLevels,
    ParseRuleCompiler,
    read_parsed_levels,
    write_parsed_levels,
)

from apb_export.api import Exporter
from conftest import (
    EXPERIMENTS,
    GROUP_OF,
    GROUPS,
    PEPTIDES,
    PRECURSORS,
    RUNS,
    DiannInput,
    MaxQuantInput,
    apb_part,
    assert_rewrite_identical,
    export_sources,
    level_part,
    observed,
    write_design,
    write_diann,
    write_maxquant,
)

# validate_prolfquapp_anndata: the keys an "lfqdata" artifact must carry.
PROLFQUAPP_KEYS = {
    "artifact_type",
    "schema_version",
    "source_software",
    "analysis_configuration",
    "protein_annotation",
}


def _parsed(diann: DiannInput) -> ParsedLevels:
    return (
        ParseRuleCompiler.detect(
            diann.report, diann.log, requested_levels=("peptidoform_ion", "protein")
        )
        .compile()
        .parse()
    )


def _annotated(parsed: ParsedLevels, column: str = "condition") -> ParsedLevels:
    """As ``apb2 annotate`` leaves a result: one more column beside every level's run key."""
    levels = {}
    for name, level in parsed.levels.items():
        frame = level.obs.frame.with_columns(pl.Series(column, ["A", "A", "B", "B"]))
        levels[name] = replace(level, obs=replace(level.obs, frame=frame))
    return replace(parsed, levels=levels)


def _written(
    target: str, parsed: ParsedLevels, folder: Path, abundance: str | None = None
) -> ad.AnnData | md.MuData:
    """Export, write the file the target reads, and read it back as the tool does."""
    export = Exporter(target, abundance)
    path = folder / f"{target}{export.extension}"
    export.write(export.export(parsed), path)
    return md.read_h5mu(path) if path.suffix == ".h5mu" else ad.read_h5ad(path)


def _anndata(value: object) -> ad.AnnData:
    assert isinstance(value, ad.AnnData)
    return value


def _frame(value: object) -> pd.DataFrame:
    assert isinstance(value, pd.DataFrame)
    return value


def test_proteopy_gets_proteins_by_sample(diann: DiannInput, tmp_path: Path) -> None:
    adata = _anndata(_written("proteopy", _annotated(_parsed(diann)), tmp_path))
    obs, var = _frame(adata.obs), _frame(adata.var)

    assert list(obs["sample_id"]) == list(adata.obs_names) == list(RUNS)
    proteins = list(var["protein_id"])
    assert proteins == list(adata.var_names)
    assert len(proteins) == len(set(proteins)), "protein ids must be unique"
    assert all(isinstance(protein, str) and protein for protein in proteins), "and present"
    assert "peptide_id" not in var.columns, "a peptide_id column would make it peptide data"
    assert list(obs.columns) == ["sample_id", "condition"], "apb2 annotate's columns pass"
    assert not np.isinf(np.asarray(adata.X)).any()


def test_proteopy_is_written_by_apb2_and_reads_back(diann: DiannInput, tmp_path: Path) -> None:
    exported = Exporter("proteopy").export(_annotated(_parsed(diann)))
    assert isinstance(exported, ParsedLevels), "ProteoPy's export is a result APB2 writes"
    first, again = tmp_path / "first.h5ad", tmp_path / "again.h5ad"

    write_parsed_levels(exported, first)
    back = read_parsed_levels(first)
    write_parsed_levels(back, again)

    a, b = ad.read_h5ad(first), ad.read_h5ad(again)
    assert _frame(a.obs).equals(_frame(b.obs)) and _frame(a.var).equals(_frame(b.var))
    np.testing.assert_array_equal(np.asarray(a.X), np.asarray(b.X))
    assert (_frame(a.obs).index.name, _frame(a.var).index.name) == ("sample_id", "protein_id")
    record = back.levels["protein"].metadata["export"]
    assert isinstance(record, dict) and isinstance(record["provenance"], dict)
    assert record["provenance"]["source_level"] == "protein", "the export's record survives"
    assert "export" in back.metadata


def test_prolfquapp_finds_its_configuration_and_columns(diann: DiannInput, tmp_path: Path) -> None:
    adata = _anndata(_written("prolfqua", _annotated(_parsed(diann)), tmp_path))
    obs, var = _frame(adata.obs), _frame(adata.var)
    meta = adata.uns["prolfquapp"]
    config = meta["analysis_configuration"]

    assert set(meta) >= PROLFQUAPP_KEYS
    assert meta["artifact_type"] == "lfqdata"
    assert meta["source_software"] == "DIA-NN"
    assert list(config["hierarchy_keys"]) == ["protein_Id", "peptide_Id", "precursor_Id"]
    assert set(config["hierarchy"]) == set(config["hierarchy_keys"])
    assert set(config["hierarchy"]) <= set(var.columns)
    assert {config["sample_name"], config["file_name"]} <= set(obs.columns)
    assert list(obs[config["sample_name"]]) == list(adata.obs_names)
    assert dict(config["factors"]) == {"condition": "condition"}
    assert list(config["work_intensity"]) == [adata.uns["X_layer_name"]]
    assert list(meta["layer_names"]) == [
        "intensity",
        "qValue",
        "pg_qValue",
        "pg_qValue_experiment",
        "pep",
    ]
    assert "qValue" in adata.layers
    first = var.iloc[0]
    keys = ("protein_Id", "peptide_Id", "precursor_Id")
    assert adata.var_names[0] == config["sep"].join(str(first[key]) for key in keys)


def test_prolfquapp_factors_come_from_each_result_alone(diann: DiannInput) -> None:
    export = Exporter("prolfqua")
    parsed = _parsed(diann)
    export.export(_annotated(parsed))

    path = diann.report.with_name("again.h5ad")
    export.write(export.export(parsed), path)
    again = ad.read_h5ad(path)

    assert again.uns["prolfquapp"]["analysis_configuration"]["factors"] == {}


def _per_cell(adata: ad.AnnData, value: Callable[[int, int], float]) -> np.ndarray:
    """``value(run, precursor)`` at each cell of a prolfqua export; NaN where none was reported."""
    groups = [group for _, group in PRECURSORS]
    precursors = [groups.index(group) for group in _frame(adata.var)["protein_Id"]]
    runs = [RUNS.index(run) for run in adata.obs_names]
    return np.array([[value(r, p) if observed(r, p) else np.nan for p in precursors] for r in runs])


def _every_cell(adata: ad.AnnData, value: Callable[[int], float]) -> np.ndarray:
    """``value(precursor)`` in every run, reported or not: one experiment-wide value per group."""
    groups = [group for _, group in PRECURSORS]
    precursors = [groups.index(group) for group in _frame(adata.var)["protein_Id"]]
    return np.array([[value(p) for p in precursors] for _run in adata.obs_names])


def test_prolfquapp_reads_protein_group_q_values_through_each_precursors_group(
    diann: DiannInput,
) -> None:
    adata = _anndata(_written("prolfqua", _parsed(diann), diann.report.parent))

    run_wise = _per_cell(adata, lambda run, precursor: 0.0001 * (precursor + 1) * (run + 1))
    # DIA-NN's global protein-group q-value is one value per group, kept in var.
    experiment = _every_cell(adata, lambda precursor: 0.002 * (precursor + 1) / 2)
    np.testing.assert_allclose(np.asarray(adata.layers["pg_qValue"]), run_wise)
    np.testing.assert_allclose(
        np.asarray(adata.layers["pg_qValue_experiment"]),
        experiment,
        err_msg="without match-between-runs the global q-value stands in for the library's",
    )


def test_prolfquapp_prefers_the_library_q_value_with_match_between_runs(tmp_path: Path) -> None:
    adata = _anndata(
        _written("prolfqua", _parsed(write_diann(tmp_path, library_q=0.003)), tmp_path)
    )

    library = _every_cell(adata, lambda precursor: 0.003 * (precursor + 1) / 2)
    np.testing.assert_allclose(np.asarray(adata.layers["pg_qValue_experiment"]), library)


def test_prolfquapp_without_a_protein_level_has_no_protein_group_q_values(
    diann: DiannInput,
) -> None:
    ion_only = ParseRuleCompiler.detect(
        diann.report, diann.log, requested_levels=("peptidoform_ion",)
    )
    adata = _anndata(_written("prolfqua", ion_only.compile().parse(), diann.report.parent))

    assert {"pg_qValue", "pg_qValue_experiment"}.isdisjoint(adata.layers)
    assert "qValue" in adata.layers


def test_prolfquapp_reads_each_precursors_pep(diann: DiannInput) -> None:
    adata = _anndata(_written("prolfqua", _parsed(diann), diann.report.parent))

    pep = _per_cell(adata, lambda _run, precursor: 0.01 * (precursor + 1))
    np.testing.assert_allclose(np.asarray(adata.layers["pep"]), pep)


def test_prolfquapp_reads_the_intensity_without_normalization(diann: DiannInput) -> None:
    adata = _anndata(_written("prolfqua", _parsed(diann), diann.report.parent))

    assert export_sources(adata)["layers.intensity"] == {
        "location": "layers",
        "name": "Precursor_Quantity",
        "catalog": "abundance",
        "kind": "intensity",
        "normalization": "none",
    }, "prolfquapp normalizes intensities itself"


def test_prolfquapp_logs_no_warning_with_an_un_normalized_layer(
    diann: DiannInput, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING, logger="apb_catalog"):
        _written("prolfqua", _parsed(diann), diann.report.parent)

    assert "SHOULD field" not in caplog.text


def test_prolfquapp_warns_and_takes_the_primary_layer_without_an_un_normalized_one(
    diann: DiannInput, caplog: pytest.LogCaptureFixture
) -> None:
    """The abundance set's SHOULD: every result should carry an un-normalized intensity."""
    parsed = _parsed(diann)
    ion = parsed.levels["peptidoform_ion"]
    layers = {name: layer for name, layer in ion.layers.items() if name != "Precursor_Quantity"}
    without = replace(parsed, levels={"peptidoform_ion": replace(ion, layers=layers)})

    with caplog.at_level(logging.WARNING, logger="apb_catalog"):
        adata = _anndata(_written("prolfqua", without, diann.report.parent))

    assert "SHOULD field kind=intensity, normalization=none on peptidoform_ion" in caplog.text
    assert export_sources(adata)["layers.intensity"] == {
        "location": "layers",
        "name": ion.primary_layer_name,
    }


def test_a_named_abundance_layer_replaces_the_rules_intensity(diann: DiannInput) -> None:
    adata = _anndata(
        _written("prolfqua", _parsed(diann), diann.report.parent, "Precursor_Normalised")
    )

    assert export_sources(adata)["layers.intensity"] == {
        "location": "layers",
        "name": "Precursor_Normalised",
    }


def _maxquant(maxquant: MaxQuantInput) -> ad.AnnData:
    """The prolfqua export of a MaxQuant folder, parsed into the levels prolfqua reads."""
    levels = Exporter("prolfqua").levels
    compiler = ParseRuleCompiler.detect(maxquant.folder, maxquant.params, requested_levels=levels)
    return _anndata(_written("prolfqua", compiler.compile().parse(), maxquant.folder.parent))


def test_prolfquapp_reads_maxquant_from_its_peptide_level(tmp_path: Path) -> None:
    adata = _maxquant(write_maxquant(tmp_path))
    var = _frame(adata.var)
    config = adata.uns["prolfquapp"]["analysis_configuration"]

    assert level_part(adata)["export"]["provenance"]["source_level"] == "peptide", (
        "peptides.txt, as prolfquapp reads it"
    )
    assert list(var.columns) == [
        "feature_id",
        "protein_Id",
        "IDcolumn",
        "fasta.id",
        "description",
        "peptide_Id",
        "nr_children",
    ]
    assert list(config["hierarchy_keys"]) == ["protein_Id", "peptide_Id"]
    assert set(config["hierarchy"]) == set(config["hierarchy_keys"])
    assert dict(zip(var["peptide_Id"], var["protein_Id"], strict=True)) == {
        sequence: GROUPS[GROUP_OF[row]][0] for row, (sequence, _, _) in enumerate(PEPTIDES)
    }, "each peptide names the group its razor protein leads"
    assert list(adata.obs_names) == list(EXPERIMENTS)
    assert {"qValue", "pg_qValue"}.isdisjoint(adata.layers), "MaxQuant reports neither per run"


def test_prolfquapp_links_maxquant_peptides_to_the_group_their_razor_protein_leads(
    tmp_path: Path,
) -> None:
    adata = _maxquant(write_maxquant(tmp_path))
    peptides = [sequence for sequence, _, _ in PEPTIDES]
    rows = [peptides.index(sequence) for sequence in _frame(adata.var)["peptide_Id"]]
    group_q = {
        peptide: 0.001 * (group + 1)
        for group, (_, members) in enumerate(GROUPS)
        for peptide in members
    }

    experiment_wide = np.array([[group_q[row] for row in rows]] * len(EXPERIMENTS))
    np.testing.assert_allclose(
        np.asarray(adata.layers["pg_qValue_experiment"]),
        experiment_wide,
        err_msg="each peptide takes the q-value of the group its razor protein leads",
    )
    pep = np.array([[0.01 * (row + 1) for row in rows]] * len(EXPERIMENTS))
    np.testing.assert_allclose(np.asarray(adata.layers["pep"]), pep)


def test_prolfquapp_reads_maxquant_evidence_alone_from_its_ion_level(tmp_path: Path) -> None:
    adata = _maxquant(write_maxquant(tmp_path, peptides=False))

    assert level_part(adata)["export"]["provenance"]["source_level"] == "peptidoform_ion"
    assert "precursor_Id" in _frame(adata.var).columns
    assert "pg_qValue_experiment" not in adata.layers


@pytest.mark.parametrize(("target", "protein"), [("prolfqua", "protein_Id")])
def test_precursors_without_a_protein_are_left_out(
    tmp_path: Path, target: str, protein: str
) -> None:
    """prolfquapp needs a protein per feature; PEAKS and DIA-NN leave some without."""
    maxquant = write_maxquant(tmp_path, peptides=False)
    levels = Exporter(target).levels
    compiler = ParseRuleCompiler.detect(maxquant.folder, maxquant.params, requested_levels=levels)
    parsed = compiler.compile().parse()
    ion = parsed.levels["peptidoform_ion"]
    column = ion.var.roles["protein_assignment"]
    unassigned = pl.when(pl.int_range(pl.len()) == 0).then(None).otherwise(pl.col(column))
    var = replace(ion.var, frame=ion.var.frame.with_columns(unassigned.alias(column)))
    result = replace(parsed, levels={**parsed.levels, "peptidoform_ion": replace(ion, var=var)})

    adata = _anndata(_written(target, result, tmp_path))

    assert adata.n_vars == ion.var.frame.height - 1
    assert adata.X is not None and adata.X.shape == (adata.n_obs, adata.n_vars)
    sources = level_part(adata)["export"]["provenance"]["sources"]
    assert sources[f"var.{protein}"]["dropped"] == "1"


def test_alphapepttools_levels_link_by_their_ids(diann: DiannInput, tmp_path: Path) -> None:
    mdata = _written("alphapepttools", _parsed(diann), tmp_path)
    assert isinstance(mdata, md.MuData)
    precursors, proteins = _anndata(mdata["precursors"]), _anndata(mdata["proteins"])
    precursor_var = _frame(precursors.var)

    assert set(mdata.mod) == {"precursors", "proteins"}, "DIA-NN gives APB2 no peptide level"
    assert precursor_var.index.name == "precursor_id"
    assert _frame(proteins.var).index.name == "proteins"
    assert {"sequence", "proteins"} <= set(precursor_var.columns)
    assert set(precursor_var["proteins"]) <= set(proteins.var_names), "every precursor links"


def _maxquant_levels(maxquant: MaxQuantInput, target: str) -> ParsedLevels:
    levels = Exporter(target).levels
    compiler = ParseRuleCompiler.detect(maxquant.folder, maxquant.params, requested_levels=levels)
    return compiler.compile().parse()


def test_alphapepttools_aligns_maxquant_levels_on_the_raw_files(tmp_path: Path) -> None:
    """evidence.txt keys runs by raw file, peptides.txt by experiment; one maps to the other."""
    parsed = _maxquant_levels(write_maxquant(tmp_path), "alphapepttools")
    mdata = _written("alphapepttools", parsed, tmp_path)

    assert isinstance(mdata, md.MuData)
    assert set(mdata.mod) == {"precursors", "peptides", "proteins"}
    raw_files = [f"raw_{name}" for name in EXPERIMENTS]
    assert list(mdata.obs_names) == raw_files, "one axis of four runs, not eight"
    for name in mdata.mod:
        assert list(_frame(mdata[name].obs)["raw_name"]) == raw_files
    assert_rewrite_identical(tmp_path / "alphapepttools.h5mu", tmp_path / "again.h5mu")


def test_alphapepttools_leaves_out_levels_it_cannot_align(tmp_path: Path) -> None:
    """Two raw files per experiment, as with fractions: only the evidence's precursors remain."""
    parsed = _maxquant_levels(write_maxquant(tmp_path), "alphapepttools")
    ion = parsed.levels["peptidoform_ion"]
    fractions = ion.obs.frame.with_columns(pl.Series("Experiment", ["S1", "S1", "S2", "S2"]))
    ion = replace(ion, obs=replace(ion.obs, frame=fractions))
    mdata = _written(
        "alphapepttools",
        replace(parsed, levels={**parsed.levels, "peptidoform_ion": ion}),
        tmp_path,
    )

    assert isinstance(mdata, md.MuData)
    assert set(mdata.mod) == {"precursors"}
    record = apb_part(mdata.uns["apb"])["export"]["provenance"]
    assert record["skipped_levels"] == ["peptide", "protein"]


def test_qfeatures_assays_link_on_the_coarser_levels_key(tmp_path: Path) -> None:
    """addAssayLink joins each assay's rowData column to the next coarser assay's key."""
    parsed = _maxquant_levels(write_maxquant(tmp_path), "qfeatures")
    mdata = _written("qfeatures", parsed, tmp_path)
    assert isinstance(mdata, md.MuData)
    precursors, peptides, proteins = (
        _frame(_anndata(mdata[name]).var) for name in ("precursors", "peptides", "proteins")
    )

    assert set(precursors["sequence"]) <= set(peptides.index), "every precursor links"
    assert set(peptides["proteins"].dropna()) <= set(proteins.index), "every peptide links"


def test_prolfquapp_file_is_written_again_unchanged(diann: DiannInput, tmp_path: Path) -> None:
    _written("prolfqua", _annotated(_parsed(diann)), tmp_path)

    assert_rewrite_identical(tmp_path / "prolfqua.h5ad", tmp_path / "again.h5ad")


def test_an_annotation_named_like_a_declared_column_is_refused(diann: DiannInput) -> None:
    clash = _annotated(_parsed(diann), column="sample_id")

    with pytest.raises(ValueError, match=r"writes columns \['sample_id'\] twice"):
        Exporter("proteopy").export(clash)


MSSTATS_COLUMNS = [
    "Run",
    "Condition",
    "BioReplicate",
    "PrecursorId",
    "ProteinName",
    "PeptideSequence",
    "PrecursorCharge",
    "FragmentIon",
    "ProductCharge",
    "IsotopeLabelType",
    "Intensity",
    "DetectionQValue",
    "PEP",
]


def _msstats(diann: DiannInput, folder: Path) -> pl.DataFrame:
    """The MSstats table of the synthetic report annotated with ``write_design``."""
    parsed = AnnotationCompiler().compile(write_design(folder)).parse(_parsed(diann)).annotate()
    export = Exporter("msstats")
    path = folder / "msstats.csv"
    export.write(export.export(parsed.parsed), path)
    return pl.read_csv(path, null_values="NA")


def test_msstats_gets_its_columns_one_row_per_precursor_and_run(
    diann: DiannInput, tmp_path: Path
) -> None:
    table = _msstats(diann, tmp_path)
    design = {run: (run[4], run[4:]) for run in RUNS}

    assert table.columns == MSSTATS_COLUMNS
    assert table.height == diann.cells
    for row in table.select("Run", "Condition", "BioReplicate").unique().iter_rows():
        assert design[row[0]] == (row[1], row[2])
    assert set(table.get_column("PrecursorCharge")) == {2}
    assert table.get_column("FragmentIon").null_count() == table.height
    assert table.get_column("ProductCharge").null_count() == table.height
    assert set(table.get_column("IsotopeLabelType")) == {"L"}
    assert "AAC[UNIMOD:4]LLK" in set(table.get_column("PeptideSequence")), "ProForma"


def _sdrf(folder: Path) -> Path:
    """An SDRF with one characteristic and two factor values per run."""
    path = folder / "samples.sdrf.tsv"
    header = (
        "source name\tcharacteristics[organism]\tcomment[data file]"
        "\tfactor value[Condition]\tfactor value[spiked compound]"
    )
    rows = [
        f"s{index}\tHomo sapiens\t{run}.raw\t{run[4]}\tyeast{1 - index % 2}"
        for index, run in enumerate(RUNS)
    ]
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return path


def _explorede_files(parsed: ParsedLevels, folder: Path) -> dict[str, ad.AnnData]:
    """Export for exploreDE and read each level's file, by its name in the directory."""
    export = Exporter("explorede")
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "explorede.h5ads"
    export.write(export.export(parsed), path)
    return {file.stem: ad.read_h5ad(file) for file in sorted(path.glob("*.h5ad"))}


def _explorede(
    parsed: ParsedLevels, annotation: Path, folder: Path, level: str = "protein"
) -> ad.AnnData:
    annotated = AnnotationCompiler().compile(annotation).parse(parsed).annotate().parsed
    return _explorede_files(annotated, folder)[level]


def test_explorede_writes_exploreDE_1_0_without_differential_expression(
    diann: DiannInput, tmp_path: Path
) -> None:
    adata = _explorede(_parsed(diann), _sdrf(tmp_path), tmp_path)
    obs, var, uns = _frame(adata.obs), _frame(adata.var), adata.uns

    assert (uns["schema_version"], uns["omics_type"], uns["feature_id_kind"]) == (
        "exploreDE/1.0",
        "proteomics",
        "protein",
    )
    assert (list(uns["hierarchy"]["levels"]), uns["hierarchy"]["unit"]) == (["protein"], "protein")
    assert (uns["de"], uns["enrichment"], uns["enrichment_support"]) == ({}, {}, {})
    assert uns["provenance"]["source_format"] == "apb-export"
    assert "prolfquapp" not in uns, "it would send the file to exploreDE's prolfquapp reader"
    assert not adata.varm, "no contrasts"
    assert list(obs.index) == list(obs["sample_label"]) == list(RUNS)
    factors = [name for name in obs.columns if name.startswith("factor_")]
    assert factors == ["factor_condition", "factor_spiked_compound"], "SDRF factor values only"
    assert isinstance(obs["factor_condition"].dtype, pd.CategoricalDtype)
    assert list(obs["factor_condition"]) == ["A", "A", "B", "B"]
    assert list(obs["factor_spiked_compound"].cat.categories) == ["yeast0", "yeast1"]
    assert "characteristics_organism" in obs.columns, "other columns pass unprefixed"
    assert list(var.index) == list(var["label_protein_id"])
    assert var.index.is_unique
    labels = dict(zip(var.index, var["label_gene_symbol"], strict=True))
    assert labels["sp|P2|P2_HUMAN"] == "P2_HUMAN", "a UniProt member by its entry name"
    assert labels["P5;CON__P6"] == "P5", "a group by its first member"


def test_explorede_intensities_are_log2_with_zeros_missing(
    diann: DiannInput, tmp_path: Path
) -> None:
    report = pl.read_csv(diann.report, separator="\t")
    zero = (pl.col("Run") == "run_A1") & (pl.col("Protein.Group") == "P1")
    report = report.with_columns(
        pl.when(zero).then(0.0).otherwise(pl.col("PG.MaxLFQ")).alias("PG.MaxLFQ")
    )
    report.write_csv(diann.report, separator="\t")
    adata = _explorede(_parsed(diann), _sdrf(tmp_path), tmp_path)
    raw = np.asarray(adata.layers["intensities_raw"])
    intensities = np.asarray(adata.layers["intensities"])

    np.testing.assert_array_equal(np.asarray(adata.X), intensities)
    np.testing.assert_allclose(intensities, np.log2(raw))
    cell = (list(adata.obs_names).index("run_A1"), list(adata.var_names).index("P1"))
    assert np.isnan(raw[cell]) and np.isnan(intensities[cell]), "a vendor zero is missing"
    assert np.isfinite(intensities[0]).sum() > 1


def test_explorede_takes_the_prolfquapp_grouping_column_as_its_factor(
    diann: DiannInput, tmp_path: Path
) -> None:
    obs = _frame(_explorede(_parsed(diann), write_design(tmp_path), tmp_path).obs)

    assert [name for name in obs.columns if name.startswith("factor_")] == ["factor_group"]
    assert "BioReplicate" in obs.columns


def _features_by_protein(adata: ad.AnnData, kind: str, levels: list[str]) -> pd.DataFrame:
    """The var of a non-protein exploreDE export, checked for its kind and feature names."""
    var = _frame(adata.var)
    assert adata.uns["feature_id_kind"] == kind
    assert list(adata.uns["hierarchy"]["levels"]) == levels
    assert var.index.is_unique
    for feature, protein in zip(var.index, var["label_protein_id"], strict=True):
        assert feature.startswith(f"{protein}_")
    return var


def test_explorede_writes_a_file_for_each_level_the_result_holds(
    diann: DiannInput, tmp_path: Path
) -> None:
    (tmp_path / "maxquant").mkdir()
    maxquant = write_maxquant(tmp_path / "maxquant")
    levels = Exporter("explorede").levels
    compiler = ParseRuleCompiler.detect(maxquant.folder, maxquant.params, requested_levels=levels)
    files = _explorede_files(compiler.compile().parse(), tmp_path / "out")

    assert levels == ("protein", "peptide", "peptidoform_ion")
    assert set(files) == {"protein", "peptide", "precursor"}
    assert set(_explorede_files(_parsed(diann), tmp_path / "diann")) == {"protein", "precursor"}


def test_explorede_names_peptides_by_protein_and_sequence(tmp_path: Path) -> None:
    maxquant = write_maxquant(tmp_path)
    levels = Exporter("explorede").levels
    compiler = ParseRuleCompiler.detect(maxquant.folder, maxquant.params, requested_levels=levels)
    adata = _explorede_files(compiler.compile().parse(), tmp_path / "out")["peptide"]

    var = _features_by_protein(adata, "peptide", ["protein", "peptide"])
    assert {feature.rsplit("_", 1)[1] for feature in var.index} <= {p[0] for p in PEPTIDES}
    obs = _frame(adata.obs)
    assert not [name for name in obs.columns if name.startswith("factor_")], (
        "MaxQuant's own columns"
    )


def test_explorede_names_precursors_by_protein_and_ion(diann: DiannInput, tmp_path: Path) -> None:
    adata = _explorede(_parsed(diann), _sdrf(tmp_path), tmp_path, "precursor")

    _features_by_protein(adata, "precursor", ["protein", "peptide", "precursor"])
    assert "P1_PEPTIDEK/2" in adata.var_names
    assert _frame(adata.var).loc["P1_PEPTIDEK/2", "label_gene_symbol"] == "P1_PEPTIDEK/2"
    assert "factor_condition" in _frame(adata.obs).columns


def test_msstats_needs_the_condition_from_an_annotation(diann: DiannInput) -> None:
    with pytest.raises(ValueError, match="Condition"):
        Exporter("msstats").export(_parsed(diann))


def test_the_api_names_each_targets_levels_and_file() -> None:
    found = {
        target: (Exporter(target).levels, Exporter(target).extension)
        for target in Exporter.targets()
    }

    assert found == {
        "alphapepttools": (("peptidoform_ion", "peptide", "protein"), ".h5mu"),
        "explorede": (("protein", "peptide", "peptidoform_ion"), ".h5ads"),
        "msmu": (("peptidoform_ion",), ".h5mu"),
        "msstats": (("peptidoform_ion",), ".csv"),
        "perseus": (("protein",), ".txt"),
        "prolfqua": (("peptide", "peptidoform_ion", "protein"), ".h5ad"),
        "proteopy": (("protein",), ".h5ad"),
        "qfeatures": (("peptidoform_ion", "peptide", "protein"), ".h5mu"),
        "qpx": (("peptidoform_ion", "protein"), ".qpx"),
    }


def test_unknown_targets_and_misplaced_abundance_are_refused() -> None:
    with pytest.raises(ValueError, match="unknown export target 'nope'"):
        Exporter("nope")
    with pytest.raises(ValueError, match="abundance applies only to targets that write one level"):
        Exporter("alphapepttools", abundance="Intensity")
