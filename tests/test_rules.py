"""Export rule documents: published schemas, the packaged rule and refused declarations."""

from __future__ import annotations

import json
from collections.abc import Callable
from importlib import resources
from typing import Any

import pytest

from apb_export.engine import CompiledExport
from apb_export.export_rules.loader import effective_rules, packaged_documents, packaged_rules
from apb_export.export_rules.schema import ExportRuleDocument
from apb_export.export_rules.schema_artifact import SCHEMA_DIRECTORY, json_schemas

type Change = Callable[[dict[str, Any]], None]


def _msmu_document() -> dict[str, Any]:
    rules = resources.files("apb_export.export_rules").joinpath("documents/msmu/v0_4/rules.json")
    return json.loads(rules.read_text(encoding="utf-8"))


def _ion(document: dict[str, Any]) -> dict[str, Any]:
    return document["tables"][0]["levels"]["ion"]


def _var(document: dict[str, Any], name: str) -> dict[str, Any]:
    (entry,) = [entry for entry in _ion(document)["columns"]["var"] if entry["name"] == name]
    return entry


def test_the_published_artifacts_are_the_schemas_the_models_declare() -> None:
    for name, schema in json_schemas().items():
        committed = json.loads((SCHEMA_DIRECTORY / name).read_text(encoding="utf-8"))
        assert committed == schema, "rerun python -m apb_export.export_rules.schema_artifact"


def test_every_packaged_rule_validates_and_compiles() -> None:
    assert packaged_documents() == (
        "alphapepttools/v0_4/rules.json",
        "msmu/v0_4/rules.json",
        "prolfqua/v2_11/rules.json",
        "proteopy/v0_1/rules.json",
    )
    extensions = {
        name.split("/", 1)[0]: CompiledExport(packaged_rules(name.split("/", 1)[0])).extension
        for name in packaged_documents()
    }

    assert extensions == {
        "alphapepttools": ".h5mu",
        "msmu": ".h5mu",
        "prolfqua": ".h5ad",
        "proteopy": ".h5ad",
    }


def test_base_blocks_reach_every_level() -> None:
    (msmu,) = packaged_rules("msmu")
    precursors, peptides, proteins = packaged_rules("alphapepttools")

    assert "msmu" in msmu.rule.accession_syntax
    assert {rule.rule.measurements.primary_layer for rule in (precursors, peptides, proteins)} == {
        "intensity"
    }


def test_an_unknown_target_is_refused() -> None:
    with pytest.raises(ValueError, match="expected one packaged export rule for 'nope'"):
        packaged_rules("nope")


def _unknown_key(document: dict[str, Any]) -> None:
    _ion(document)["colour"] = "red"


def _undeclared_primary(document: dict[str, Any]) -> None:
    _ion(document)["measurements"]["primary_layer"] = "missing"


def _undeclared_axis(document: dict[str, Any]) -> None:
    _ion(document)["axis"]["var_keys"] = ["nope"]


def _undeclared_syntax(document: dict[str, Any]) -> None:
    document["tables"][0]["base"] = {}


def _width_on_text(document: dict[str, Any]) -> None:
    _var(document, "peptide")["width"] = 32


def _index_only_column(document: dict[str, Any]) -> None:
    _var(document, "peptide")["index_only"] = True


def _no_modality(document: dict[str, Any]) -> None:
    del _ion(document)["modality"]


def _repeated_name(document: dict[str, Any]) -> None:
    _ion(document)["columns"]["var"].append(_var(document, "peptide"))


def _map_in_columns(document: dict[str, Any]) -> None:
    _ion(document)["columns"]["var"].append({"name": "factors", "how": "annotation_map"})


def _mapping_constant_in_columns(document: dict[str, Any]) -> None:
    _ion(document)["columns"]["var"].append({"name": "x", "how": "constant", "value": {"a": 1}})


def _named_column(document: dict[str, Any]) -> None:
    _var(document, "peptide")["named_index"] = True


def _two_tables(document: dict[str, Any]) -> None:
    document["tables"].append(json.loads(json.dumps(document["tables"][0])))


def _h5ad_with_two_levels(document: dict[str, Any]) -> None:
    table = document["tables"][0]
    table["output"]["extensions"] = [".h5ad"]
    table["levels"]["protein"] = json.loads(json.dumps(_ion(document)))


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (_unknown_key, "Extra inputs are not permitted"),
        (_undeclared_primary, "primary_layer='missing' matches no layer"),
        (_undeclared_axis, r"axis.var_keys=\['nope'\] matches no columns.var entry"),
        (_undeclared_syntax, "names no accession_syntax 'msmu'"),
        (_width_on_text, "width applies only to type number"),
        (_index_only_column, "'peptide' is index_only but no axis key"),
        (_no_modality, "writes .h5mu but names no modality"),
        (_repeated_name, "columns.var entry names must be unique"),
        (_map_in_columns, "'factors' can only be written to uns"),
        (_mapping_constant_in_columns, "'x' can only be written to uns"),
        (_named_column, "'peptide' is named_index but no axis key"),
        (_two_tables, "an export rule writes one table, not 2"),
        (_h5ad_with_two_levels, ".h5ad holds one level"),
    ],
)
def test_inconsistent_declarations_are_refused(change: Change, message: str) -> None:
    document = _msmu_document()
    change(document)

    with pytest.raises(ValueError, match=message):
        effective_rules(ExportRuleDocument.model_validate(document), "test")
