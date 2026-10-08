"""Each packaged target's side of its consumer's contract, checked without the consumer.

Whether the consumers really read these files is checked in the consumer image; here the
written files must carry what each consumer's reader looks for.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import polars as pl
import pytest
from apb2.api import ParsedLevels, ParseRuleCompiler

from apb_export.api import Exporter
from conftest import (
    EXPERIMENTS,
    GROUPS,
    PEPTIDES,
    PRECURSORS,
    RUNS,
    DiannInput,
    MaxQuantInput,
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
        ParseRuleCompiler(diann.report, diann.log, requested_levels=("ion", "protein"))
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


def _written(target: str, parsed: ParsedLevels, folder: Path) -> ad.AnnData | md.MuData:
    """Export, write the file the target reads, and read it back."""
    export = Exporter(target)
    result = export.export(parsed)
    path = folder / f"{target}{export.extension}"
    if isinstance(result, md.MuData):
        result.write_h5mu(path)
        return md.read_h5mu(path)
    result.write_h5ad(path)
    return ad.read_h5ad(path)


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

    again = _anndata(export.export(parsed))

    assert again.uns["prolfquapp"]["analysis_configuration"]["factors"] == {}


def _per_cell(adata: ad.AnnData, value: Callable[[int, int], float]) -> np.ndarray:
    """``value(run, precursor)`` at each cell of a prolfqua export; NaN where none was reported."""
    groups = [group for _, group in PRECURSORS]
    precursors = [groups.index(group) for group in _frame(adata.var)["protein_Id"]]
    runs = [RUNS.index(run) for run in adata.obs_names]
    return np.array([[value(r, p) if observed(r, p) else np.nan for p in precursors] for r in runs])


def test_prolfquapp_reads_protein_group_q_values_through_each_precursors_group(
    diann: DiannInput,
) -> None:
    adata = _anndata(Exporter("prolfqua").export(_parsed(diann)))

    run_wise = _per_cell(adata, lambda run, precursor: 0.0001 * (precursor + 1) * (run + 1))
    experiment = _per_cell(adata, lambda run, precursor: 0.002 * (precursor + 1) / 2)
    np.testing.assert_allclose(np.asarray(adata.layers["pg_qValue"]), run_wise)
    np.testing.assert_allclose(
        np.asarray(adata.layers["pg_qValue_experiment"]),
        experiment,
        err_msg="without match-between-runs the global q-value stands in for the library's",
    )


def test_prolfquapp_prefers_the_library_q_value_with_match_between_runs(tmp_path: Path) -> None:
    adata = _anndata(Exporter("prolfqua").export(_parsed(write_diann(tmp_path, library_q=0.003))))

    library = _per_cell(adata, lambda run, precursor: 0.003 * (precursor + 1) / 2)
    np.testing.assert_allclose(np.asarray(adata.layers["pg_qValue_experiment"]), library)


def test_prolfquapp_without_a_protein_level_has_no_protein_group_q_values(
    diann: DiannInput,
) -> None:
    ion_only = ParseRuleCompiler(diann.report, diann.log, requested_levels=("ion",))
    adata = _anndata(Exporter("prolfqua").export(ion_only.compile().parse()))

    assert {"pg_qValue", "pg_qValue_experiment"}.isdisjoint(adata.layers)
    assert "qValue" in adata.layers


def test_prolfquapp_reads_each_precursors_pep(diann: DiannInput) -> None:
    adata = _anndata(Exporter("prolfqua").export(_parsed(diann)))

    pep = _per_cell(adata, lambda _run, precursor: 0.01 * (precursor + 1))
    np.testing.assert_allclose(np.asarray(adata.layers["pep"]), pep)


def _maxquant(maxquant: MaxQuantInput) -> ad.AnnData:
    """The prolfqua export of a MaxQuant folder, parsed into the levels prolfqua reads."""
    levels = Exporter("prolfqua").levels
    compiler = ParseRuleCompiler(maxquant.folder, maxquant.params, requested_levels=levels)
    return _anndata(Exporter("prolfqua").export(compiler.compile().parse()))


def test_prolfquapp_reads_maxquant_from_its_peptide_level(tmp_path: Path) -> None:
    adata = _maxquant(write_maxquant(tmp_path))
    var = _frame(adata.var)
    config = adata.uns["prolfquapp"]["analysis_configuration"]

    assert level_part(adata)["export"]["provenance"]["source_level"] == "peptide", (
        "peptides.txt, as prolfquapp reads it"
    )
    assert list(var.columns) == [
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

    assert level_part(adata)["export"]["provenance"]["source_level"] == "ion"
    assert "precursor_Id" in _frame(adata.var).columns
    assert "pg_qValue_experiment" not in adata.layers


def test_proteopy_reads_precursors_as_peptides_without_a_protein_level(tmp_path: Path) -> None:
    """ProteoPy's peptide-level data, as its own DIA-NN reader writes precursors."""
    maxquant = write_maxquant(tmp_path, peptides=False)
    levels = Exporter("proteopy").levels
    compiler = ParseRuleCompiler(maxquant.folder, maxquant.params, requested_levels=levels)
    adata = _anndata(Exporter("proteopy").export(compiler.compile().parse()))
    var = _frame(adata.var)

    assert level_part(adata)["export"]["provenance"]["source_level"] == "ion"
    peptides = list(var["peptide_id"])
    assert peptides == list(adata.var_names)
    assert len(peptides) == len(set(peptides)), "peptide ids must be unique"
    proteins = list(var["protein_id"])
    assert all(isinstance(protein, str) and protein for protein in proteins), "one protein each"


@pytest.mark.parametrize(
    ("target", "protein"), [("proteopy", "protein_id"), ("prolfqua", "protein_Id")]
)
def test_precursors_without_a_protein_are_left_out(
    tmp_path: Path, target: str, protein: str
) -> None:
    """ProteoPy and prolfquapp need a protein per feature; PEAKS and DIA-NN leave some without."""
    maxquant = write_maxquant(tmp_path, peptides=False)
    levels = Exporter(target).levels
    compiler = ParseRuleCompiler(maxquant.folder, maxquant.params, requested_levels=levels)
    parsed = compiler.compile().parse()
    ion = parsed.levels["ion"]
    column = ion.var.roles["protein_assignment"]
    unassigned = pl.when(pl.int_range(pl.len()) == 0).then(None).otherwise(pl.col(column))
    var = replace(ion.var, frame=ion.var.frame.with_columns(unassigned.alias(column)))
    result = replace(parsed, levels={**parsed.levels, "ion": replace(ion, var=var)})

    adata = _anndata(Exporter(target).export(result))

    assert adata.n_vars == ion.var.frame.height - 1
    assert adata.X is not None and adata.X.shape == (adata.n_obs, adata.n_vars)
    sources = level_part(adata)["export"]["provenance"]["sources"]
    assert sources[f"var.{protein}"]["dropped"] == "1"


def test_proteopy_leaves_out_a_protein_group_without_a_name(diann: DiannInput) -> None:
    """DIA-NN's protein level holds an unnamed group when precursors lack Protein.Group."""
    parsed = _parsed(diann)
    protein = parsed.levels["protein"]
    column = protein.var.roles["protein_assignment"]
    unnamed = pl.when(pl.int_range(pl.len()) == 0).then(pl.lit("")).otherwise(pl.col(column))
    var = replace(protein.var, frame=protein.var.frame.with_columns(unnamed.alias(column)))
    result = replace(parsed, levels={**parsed.levels, "protein": replace(protein, var=var)})

    adata = _anndata(Exporter("proteopy").export(result))

    assert adata.n_vars == protein.var.frame.height - 1
    assert "" not in set(_frame(adata.var)["protein_id"])


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
        "alphapepttools": (("ion", "peptide", "protein"), ".h5mu"),
        "msmu": (("ion",), ".h5mu"),
        "prolfqua": (("peptide", "ion", "protein"), ".h5ad"),
        "proteopy": (("protein", "ion"), ".h5ad"),
    }


def test_unknown_targets_and_misplaced_abundance_are_refused() -> None:
    with pytest.raises(ValueError, match="unknown export target 'nope'"):
        Exporter("nope")
    with pytest.raises(ValueError, match="abundance applies only to targets that write one level"):
        Exporter("alphapepttools", abundance="Intensity")
