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
from apb2.api import ParsedLevels, ParseRuleCompiler, read_parsed_levels, write_parsed_levels

from apb_export.api import Exporter
from conftest import (
    EXPERIMENTS,
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
        sequence: razor for sequence, _, razor in PEPTIDES
    }
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


def test_prolfquapp_file_is_written_again_unchanged(diann: DiannInput, tmp_path: Path) -> None:
    _written("prolfqua", _annotated(_parsed(diann)), tmp_path)

    assert_rewrite_identical(tmp_path / "prolfqua.h5ad", tmp_path / "again.h5ad")


def test_an_annotation_named_like_a_declared_column_is_refused(diann: DiannInput) -> None:
    clash = _annotated(_parsed(diann), column="sample_id")

    with pytest.raises(ValueError, match=r"writes columns \['sample_id'\] twice"):
        Exporter("proteopy").export(clash)


def test_the_api_names_each_targets_levels_and_file() -> None:
    found = {
        target: (Exporter(target).levels, Exporter(target).extension)
        for target in Exporter.targets()
    }

    assert found == {
        "alphapepttools": (("peptidoform_ion", "peptide", "protein"), ".h5mu"),
        "msmu": (("peptidoform_ion",), ".h5mu"),
        "prolfqua": (("peptide", "peptidoform_ion", "protein"), ".h5ad"),
        "proteopy": (("protein",), ".h5ad"),
        "qpx": (("peptidoform_ion", "protein"), ".qpx"),
    }


def test_unknown_targets_and_misplaced_abundance_are_refused() -> None:
    with pytest.raises(ValueError, match="unknown export target 'nope'"):
        Exporter("nope")
    with pytest.raises(ValueError, match="abundance applies only to targets that write one level"):
        Exporter("alphapepttools", abundance="Intensity")
