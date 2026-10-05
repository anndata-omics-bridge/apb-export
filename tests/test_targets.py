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
from conftest import PRECURSORS, RUNS, DiannInput, observed, write_diann

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
    assert set(config["hierarchy"]) == {"protein_Id", "peptide_Id", "precursor_Id"}
    assert set(config["hierarchy"]) <= set(var.columns)
    assert {config["sample_name"], config["file_name"]} <= set(obs.columns)
    assert list(obs[config["sample_name"]]) == list(adata.obs_names)
    assert dict(config["factors"]) == {"condition": "condition"}
    assert list(config["work_intensity"]) == [adata.uns["X_layer_name"]]
    assert list(meta["layer_names"]) == ["intensity", "qValue", "pg_qValue", "pg_qValue_experiment"]
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
        "prolfqua": (("ion", "protein"), ".h5ad"),
        "proteopy": (("protein",), ".h5ad"),
    }


def test_unknown_targets_and_misplaced_abundance_are_refused() -> None:
    with pytest.raises(ValueError, match="unknown export target 'nope'"):
        Exporter("nope")
    with pytest.raises(ValueError, match="abundance applies only to targets that write one level"):
        Exporter("alphapepttools", abundance="Intensity")
