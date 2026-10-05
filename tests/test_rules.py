"""Export rule documents: published schemas, the packaged rule and refused declarations."""

from __future__ import annotations

import json
from collections.abc import Callable
from importlib import resources
from typing import Any

import pytest

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


def test_the_packaged_msmu_rule_validates() -> None:
    assert packaged_documents() == ("msmu/v0_4/rules.json",)
    (rule,) = packaged_rules("msmu")

    assert (rule.target_name, rule.level, rule.rule.modality) == ("msmu", "ion", "psm")
    assert "msmu" in rule.rule.accession_syntax, "base blocks reach the level"


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
    ],
)
def test_inconsistent_declarations_are_refused(change: Change, message: str) -> None:
    document = _msmu_document()
    change(document)

    with pytest.raises(ValueError, match=message):
        effective_rules(ExportRuleDocument.model_validate(document), "test")
