"""The rule engine: q-value candidates, absent sources, layer minima and misplaced sources."""

from __future__ import annotations

import json
import math
from importlib import resources
from pathlib import Path
from typing import Any

import anndata as ad
import mudata as md
import pandas as pd
import pytest
from apb2.api import AnnotationCompiler, ParsedLevels, ParseRuleCompiler
from apb_catalog.api import UnresolvedField

from apb_export.api import Exporter
from apb_export.engine import CompiledExport
from apb_export.export_rules.loader import effective_rules
from apb_export.export_rules.schema import ExportRuleDocument
from conftest import RUNS, DiannInput, export_sources, read_back, write_diann

# The confidence set catalogues no per-precursor PEP for DIA-NN, so this candidate never resolves.
UNOFFERED = {"catalog": "confidence", "kind": "pep", "location": "var"}
RUN_Q_VALUE = {"catalog": "confidence", "kind": "q_value"}


def _parsed(diann: DiannInput) -> ParsedLevels:
    return ParseRuleCompiler.detect(diann.report, diann.log, ("peptidoform_ion",)).compile().parse()


def _msmu_document() -> dict[str, Any]:
    rules = resources.files("apb_export.export_rules").joinpath("documents/msmu/v0_4/rules.json")
    return json.loads(rules.read_text(encoding="utf-8"))


def _ion(document: dict[str, Any]) -> dict[str, Any]:
    return document["tables"][0]["levels"]["peptidoform_ion"]


def _compiled(document: dict[str, Any]) -> CompiledExport:
    return CompiledExport(effective_rules(ExportRuleDocument.model_validate(document), "test"))


def _with_var(replacements: dict[str, dict[str, Any]]) -> CompiledExport:
    """msmu's rule with var entries replaced by name."""
    document = _msmu_document()
    columns = _ion(document)["columns"]
    columns["var"] = [replacements.get(entry["name"], entry) for entry in columns["var"]]
    return _compiled(document)


def _psm(compiled: CompiledExport, diann: DiannInput) -> ad.AnnData:
    return _modality(compiled.export(_parsed(diann)), diann.report.with_name("psm.h5mu"))


def _modality(exported: ParsedLevels, path: Path) -> ad.AnnData:
    """msmu's psm modality, as msmu reads it from the written file."""
    container = read_back(exported, path)
    assert isinstance(container, md.MuData)
    psm = container["psm"]
    assert isinstance(psm, ad.AnnData)
    return psm


def _var(psm: object) -> pd.DataFrame:
    """The psm ``var`` frame; MuData and AnnData type modalities and frames as unions."""
    assert isinstance(psm, ad.AnnData)
    assert isinstance(psm.var, pd.DataFrame)
    return psm.var


def test_an_all_zero_global_q_value_leaves_the_run_q_value_to_stand_in(tmp_path: Path) -> None:
    """apb2 creates no column whose every value is the same; the next candidate answers."""
    diann = write_diann(tmp_path, global_q=0.0)
    parsed = _parsed(diann)
    psm = _modality(Exporter("msmu").export(parsed), tmp_path / "psm.h5mu")

    assert "Global_Q_Value" not in parsed.levels["peptidoform_ion"].var.frame.columns
    assert export_sources(psm)["var.q_value"]["kind"] == "q_value"


def test_a_candidate_the_catalogue_does_not_offer_falls_through(diann: DiannInput) -> None:
    q_value = {
        "name": "q_value",
        "how": "layer_minimum",
        "source": [UNOFFERED, RUN_Q_VALUE],
        "type": "number",
        "width": 32,
        "required": False,
    }
    psm = _psm(_with_var({"q_value": q_value}), diann)

    assert export_sources(psm)["var.q_value"]["kind"] == "q_value"
    assert _var(psm).loc["PEPTIDEK/2", "q_value"] == pytest.approx(0.001), "Q.Value, not global"


# A prolfquapp table with both a group and a bait column: apb-catalog's condition is ambiguous.
CONDITION = {"catalog": "apb", "kind": "condition", "location": "obs"}


def _with_condition(source: object) -> CompiledExport:
    """msmu's rule with a ``condition`` obs entry read from ``source``."""
    document = _msmu_document()
    _ion(document)["columns"]["obs"].append({"name": "condition", "source": source})
    return _compiled(document)


def _two_conditions(diann: DiannInput) -> ParsedLevels:
    """The DIA-NN result annotated with a table that names a condition twice."""
    design = diann.report.with_name("design.tsv")
    rows = [f"{run}\t{run[4]}\tbait_{run[4]}" for run in RUNS]
    design.write_text("\n".join(["raw_file\tgroup\tbait", *rows]) + "\n", encoding="utf-8")
    return AnnotationCompiler().compile(design).parse(_parsed(diann)).annotate().parsed


def test_a_candidate_the_catalogue_cannot_pick_passes_to_the_next(diann: DiannInput) -> None:
    compiled = _with_condition([CONDITION, {"axis": "obs_keys"}])
    psm = _modality(compiled.export(_two_conditions(diann)), diann.report.with_name("psm.h5mu"))

    assert export_sources(psm)["obs.condition"]["passed"] == "apb condition ambiguous"


def test_the_last_candidate_the_catalogue_cannot_pick_is_an_error(diann: DiannInput) -> None:
    with pytest.raises(UnresolvedField, match="ambiguous"):
        _with_condition(CONDITION).export(_two_conditions(diann))


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
    assert export_sources(psm)["var.PEP"] == {"location": "absent", "filled": "NaN"}
    assert export_sources(psm)["var.q_value"] == {"location": "absent"}


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
    same = {"name": "precursor", "how": "constant", "value": "same"}

    with pytest.raises(ValueError, match="var index 'precursor' repeats values"):
        _with_var({"precursor": same}).export(_parsed(diann))


def test_a_layer_minimum_leaves_out_missing_cells(diann: DiannInput) -> None:
    lowest = {
        "name": "lowest",
        "how": "layer_minimum",
        "source": {"measurements": "primary_layer"},
        "type": "number",
    }
    document = _msmu_document()
    _ion(document)["columns"]["var"].append(lowest)

    var = _var(_psm(_compiled(document), diann))

    # Precursor_Normalised is 900 x run x precursor; SHAREDK is absent from run_A1.
    assert var.loc["PEPTIDEK/2", "lowest"] == pytest.approx(900.0)
    assert var.loc["SHAREDK/2", "lowest"] == pytest.approx(900.0 * 2 * 5)
    assert export_sources(_psm(_compiled(document), diann))["var.lowest"]["how"] == "layer_minimum"
