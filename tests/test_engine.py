"""The rule engine: q-value candidates, absent sources, selected cells and misplaced sources."""

from __future__ import annotations

import json
import math
from importlib import resources
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
import polars as pl
import pytest
from apb2.api import ParsedLevels, ParseRuleCompiler

from apb_export.api import MsmuExporter
from apb_export.engine import CompiledExport
from apb_export.export_rules.loader import effective_rules
from apb_export.export_rules.schema import ExportRuleDocument
from apb_export.rows import select_cells
from conftest import DiannInput, write_diann

# identification_confidence catalogues no library q-value, so this candidate never resolves.
UNOFFERED = {
    "catalogue": "identification_confidence",
    "concept": "confidence",
    "kind": "library_q_value",
}
RUN_Q_VALUE = {"catalogue": "identification_confidence", "concept": "confidence", "kind": "q_value"}


def _parsed(diann: DiannInput) -> ParsedLevels:
    return ParseRuleCompiler(diann.report, diann.log, requested_levels=("ion",)).compile().parse()


def _msmu_document() -> dict[str, Any]:
    rules = resources.files("apb_export.export_rules").joinpath("documents/msmu/v0_4/rules.json")
    return json.loads(rules.read_text(encoding="utf-8"))


def _ion(document: dict[str, Any]) -> dict[str, Any]:
    return document["tables"][0]["levels"]["ion"]


def _compiled(document: dict[str, Any]) -> CompiledExport:
    return CompiledExport(effective_rules(ExportRuleDocument.model_validate(document), "test"))


def _with_var(replacements: dict[str, dict[str, Any]]) -> CompiledExport:
    """msmu's rule with var entries replaced by name."""
    document = _msmu_document()
    columns = _ion(document)["columns"]
    columns["var"] = [replacements.get(entry["name"], entry) for entry in columns["var"]]
    return _compiled(document)


def _psm(compiled: CompiledExport, diann: DiannInput) -> ad.AnnData:
    psm = compiled.export(_parsed(diann))["psm"]
    assert isinstance(psm, ad.AnnData)
    return psm


def _var(psm: object) -> pd.DataFrame:
    """The psm ``var`` frame; MuData and AnnData type modalities and frames as unions."""
    assert isinstance(psm, ad.AnnData)
    assert isinstance(psm.var, pd.DataFrame)
    return psm.var


def test_zero_global_q_values_stay_selected(tmp_path: Path) -> None:
    diann = write_diann(tmp_path, global_q=0.0)
    psm = MsmuExporter().export(_parsed(diann))["psm"]

    assert psm.uns["apb"]["sources"]["var.q_value"]["kind"] == "global_q_value"
    assert (_var(psm)["q_value"] == 0).all(), "zero is a value; only the library needs nonzero"


def test_a_candidate_the_catalogue_does_not_offer_falls_through(diann: DiannInput) -> None:
    q_value = {
        "name": "q_value",
        "source": [UNOFFERED, RUN_Q_VALUE],
        "type": "number",
        "width": 32,
        "required": False,
    }
    psm = _psm(_with_var({"q_value": q_value}), diann)

    assert psm.uns["apb"]["sources"]["var.q_value"]["kind"] == "q_value"
    assert _var(psm).loc["run_A1.PEPTIDEK/2", "q_value"] == pytest.approx(0.001), "per-run Q.Value"


def test_absent_optional_sources_are_filled_or_left_out(diann: DiannInput) -> None:
    pep = {
        "name": "PEP",
        "source": UNOFFERED,
        "type": "number",
        "width": 32,
        "required": False,
        "fill_value": "NaN",
    }
    q_value = {"name": "q_value", "source": UNOFFERED, "type": "number", "required": False}
    psm = _psm(_with_var({"PEP": pep, "q_value": q_value}), diann)

    assert "q_value" not in _var(psm).columns, "msmu's filter must not see an invented q-value"
    assert all(math.isnan(value) for value in _var(psm)["PEP"]), "msmu's to_peptide needs PEP"
    assert psm.uns["apb"]["sources"]["var.PEP"] == {"location": "absent", "filled": "NaN"}
    assert psm.uns["apb"]["sources"]["var.q_value"] == {"location": "absent"}


def test_an_absent_required_source_is_an_error(diann: DiannInput) -> None:
    pep = {"name": "PEP", "source": UNOFFERED, "type": "number"}

    with pytest.raises(ValueError, match="var entry 'PEP': the result carries no source for it"):
        _with_var({"PEP": pep}).export(_parsed(diann))


def test_a_layer_cannot_fill_a_uns_entry(diann: DiannInput) -> None:
    document = _msmu_document()
    _ion(document)["uns"].append({"name": "bad", "source": {"measurements": "primary_layer"}})

    with pytest.raises(ValueError, match="uns entry 'bad' cannot take layers"):
        _compiled(document).export(_parsed(diann))


def test_an_index_must_not_repeat(diann: DiannInput) -> None:
    runs_only = {"name": "feature_id", "source": {"axis": "obs_keys"}, "index_only": True}

    with pytest.raises(ValueError, match="var index 'feature_id' repeats values"):
        _with_var({"feature_id": runs_only}).export(_parsed(diann))


def test_long_output_writes_only_its_primary_layer() -> None:
    document = _msmu_document()
    layers = _ion(document)["measurements"]["layers"]
    layers.append({"name": "second", "source": {"measurements": "primary_layer"}})

    with pytest.raises(ValueError, match="long output writes only its primary layer"):
        _compiled(document)


def test_cells_exclude_zero_null_nan_and_infinity(diann: DiannInput) -> None:
    level = _parsed(diann).levels["ion"]
    rows, runs = level.var.frame.height, level.obs.frame.height
    values: list[list[float | None]] = [[1.0] * runs for _ in range(rows)]
    for position, value in enumerate((0.0, None, math.nan, math.inf)):
        values[position][position] = value
    x = pl.DataFrame(
        {f"run{run}": [values[row][run] for row in range(rows)] for run in range(runs)},
        schema={f"run{run}": pl.Float64 for run in range(runs)},
    )

    cells = select_cells(level, x, [0.0])

    kept = set(zip(cells.var_rows.tolist(), cells.obs.tolist(), strict=True))
    assert len(kept) == rows * runs - 4
    assert not kept & {(0, 0), (1, 1), (2, 2), (3, 3)}
    keys = level.var.frame.get_column("ProForma_ion").to_list()
    order = [(run, keys[row]) for row, run in zip(cells.var_rows, cells.obs, strict=True)]
    assert order == sorted(order), "rows run by observation, then by the variable key"
