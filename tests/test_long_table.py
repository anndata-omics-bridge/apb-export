"""The long-table output: one row per observed precursor and run, the layers side by side."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import polars as pl
import pytest
from apb2.api import ParsedLevel, ParsedLevels, ParseRuleCompiler

from apb_export.engine import CompiledExport
from apb_export.export_rules.loader import effective_rules
from apb_export.export_rules.schema import ExportRuleDocument
from apb_export.long_table import long_table
from conftest import FRAGPIPE_IONS, FRAGPIPE_SAMPLES, DiannInput, write_fragpipe

COLUMNS = ["Run", "precursor", "proteins", "Intensity", "q_value"]


def _document(var: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """A one-level ``.csv`` rule: the run key, two var columns, the primary layer and a q-value."""
    return {
        "schema_version": "0.4",
        "file_version": "1",
        "target_name": "Table",
        "target_version_pattern": "^1\\.",
        "tables": [
            {
                "output": {"extensions": [".csv"]},
                "levels": {
                    "peptidoform_ion": {
                        "axis": {"obs_keys": ["Run"], "var_keys": ["precursor"]},
                        "measurements": {
                            "primary_layer": "Intensity",
                            "layers": [
                                {"name": "Intensity", "source": {"measurements": "primary_layer"}},
                                {
                                    "name": "q_value",
                                    "source": {
                                        "catalog": "confidence",
                                        "kind": "q_value",
                                        "qualifiers": {"stage": "identification"},
                                    },
                                },
                            ],
                        },
                        "columns": {
                            "obs": [{"name": "Run", "source": {"axis": "obs_keys"}}],
                            "var": var
                            or [
                                {
                                    "name": "precursor",
                                    "source": {
                                        "catalog": "apb",
                                        "kind": "peptidoform_ion",
                                        "location": "var",
                                    },
                                },
                                {
                                    "name": "proteins",
                                    "source": {
                                        "catalog": "apb",
                                        "kind": "protein_assignment",
                                        "location": "var",
                                    },
                                },
                            ],
                        },
                    }
                },
            }
        ],
    }


def _compiled(document: dict[str, Any]) -> CompiledExport:
    return CompiledExport(effective_rules(ExportRuleDocument.model_validate(document), "test"))


def _exported(diann: DiannInput, document: dict[str, Any] | None = None) -> ParsedLevel:
    parsed = ParseRuleCompiler.detect(diann.report, diann.log, ("peptidoform_ion",)).compile()
    exported = _compiled(document or _document()).export(parsed.parse())
    (level,) = exported.levels.values()
    return level


def _without(level: ParsedLevel, layer: str, rows: pl.Expr) -> ParsedLevel:
    """The level with the cells of one layer that ``rows`` selects removed."""
    table = level.layers[layer]
    kept = replace(table, cells=table.cells.filter(~rows))
    return replace(level, layers={**level.layers, layer: kept})


def test_one_row_per_observed_cell_in_obs_var_layer_order(diann: DiannInput) -> None:
    level = _exported(diann)
    table = long_table(level)
    intensity = level.layers["Intensity"].cells.filter(pl.col("value").is_not_null())

    assert table.columns == COLUMNS
    assert table.height == intensity.height == diann.cells
    assert table.get_column("Intensity").null_count() == 0
    assert sorted(table.get_column("Intensity").to_list()) == sorted(
        intensity.get_column("value").to_list()
    )
    assert set(table.get_column("Run")) == set(level.obs.frame.get_column("Run"))


def test_a_cell_missing_in_every_layer_is_dropped_and_in_one_kept(diann: DiannInput) -> None:
    level = _exported(diann)
    first = (pl.col("var_pos") == 0) & (pl.col("obs_pos") == 0)
    second = (pl.col("var_pos") == 0) & (pl.col("obs_pos") == 1)
    level = _without(level, "Intensity", first | second)
    level = _without(level, "q_value", first)

    table = long_table(level)
    precursor = level.var.frame.get_column("precursor")[0]
    runs = level.obs.frame.get_column("Run")
    rows = table.filter(pl.col("precursor") == precursor)

    assert table.height == diann.cells - 1
    assert runs[0] not in rows.get_column("Run").to_list(), "no value in any layer"
    kept = rows.filter(pl.col("Run") == runs[1])
    assert kept.get_column("Intensity").to_list() == [None]
    assert kept.get_column("q_value").null_count() == 0


def test_the_csv_reads_back_with_na_for_missing(diann: DiannInput, tmp_path: Path) -> None:
    compiled = _compiled(_document())
    parsed = ParseRuleCompiler.detect(diann.report, diann.log, ("peptidoform_ion",)).compile()
    path = tmp_path / "table.csv"
    compiled.write(compiled.export(parsed.parse()), path)

    back = pl.read_csv(path, null_values="NA")
    assert compiled.extension == ".csv"
    assert compiled.writes_one_level
    assert back.columns == COLUMNS
    assert back.height == diann.cells
    assert not list(tmp_path.glob("*.apb.json")), "a table has no sidecar"


def test_a_name_used_twice_is_refused(diann: DiannInput) -> None:
    precursor = {"catalog": "apb", "kind": "peptidoform_ion", "location": "var"}
    clash = [{"name": "precursor", "source": precursor}, {"name": "Run", "source": precursor}]
    level = _exported(diann, _document(var=clash))

    with pytest.raises(ValueError, match=r"would name \['Run'\] twice"):
        long_table(level)


def test_a_table_of_several_levels_is_refused(diann: DiannInput, tmp_path: Path) -> None:
    compiled = _compiled(_document())
    parsed = ParseRuleCompiler.detect(diann.report, diann.log, ("peptidoform_ion",)).compile()
    exported = compiled.export(parsed.parse())
    (level,) = exported.levels.values()
    two: ParsedLevels = replace(exported, levels={"a": level, "b": level})

    with pytest.raises(ValueError, match="holds one level"):
        compiled.write(two, tmp_path / "table.csv")


def test_csv_levels_before_the_last_name_their_software() -> None:
    document = _document()
    levels = document["tables"][0]["levels"]
    levels["peptide"] = levels["peptidoform_ion"]

    with pytest.raises(ValueError, match=r"\.csv levels \['peptidoform_ion'\] precede another"):
        _compiled(document)


def test_a_categorical_layer_is_written_as_its_labels(tmp_path: Path) -> None:
    """FragPipe's match types are stored as codes; the table carries the rule's labels."""
    source = write_fragpipe(tmp_path)
    compiler = ParseRuleCompiler.from_software(
        source, software="fragpipe", requested_levels=("peptidoform_ion",)
    )
    level = compiler.compile().parse().levels["peptidoform_ion"]

    table = long_table(level)
    found = {
        (row["Peptide_Sequence"], row["sample"]): row["Match_Type"]
        for row in table.iter_rows(named=True)
    }

    assert found == {
        (sequence, sample): match
        for sequence, _, matches in FRAGPIPE_IONS
        for sample, match in zip(FRAGPIPE_SAMPLES, matches, strict=True)
    }
