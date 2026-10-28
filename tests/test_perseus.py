"""The Perseus table: proteinGroups.txt's names, its sample columns typed main, header rows."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import polars as pl
import pytest
from apb2.api import AnnotationCompiler, ParsedLevels, ParseRuleCompiler, sidecar_path

from apb_export.api import Exporter
from apb_export.perseus_table import perseus_table
from conftest import EXPERIMENTS, GROUPS, RUNS, DiannInput, write_design, write_maxquant

VAR_COLUMNS = ["Protein IDs", "Q-value", "Reverse", "Potential contaminant"]


def _maxquant(folder: Path) -> ParsedLevels:
    """The protein level of the synthetic MaxQuant folder, with LFQ intensities."""
    maxquant = write_maxquant(folder, lfq=True)
    compiler = ParseRuleCompiler.detect(maxquant.folder, maxquant.params, ("protein",))
    return compiler.compile().parse()


def _diann(diann: DiannInput) -> ParsedLevels:
    return ParseRuleCompiler.detect(diann.report, diann.log, ("protein",)).compile().parse()


def _write(parsed: ParsedLevels, path: Path, abundance: str | None = None) -> list[list[str]]:
    """Export and write the Perseus table; its lines split at tabs."""
    exporter = Exporter("perseus", abundance=abundance)
    exporter.write(exporter.export(parsed), path)
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]


def test_maxquant_proteins_get_proteingroups_names_and_types(tmp_path: Path) -> None:
    lines = _write(_maxquant(tmp_path), tmp_path / "perseus.txt")
    samples = [f"Intensity {name}" for name in EXPERIMENTS]

    assert lines[0] == [*VAR_COLUMNS, *samples]
    assert lines[1] == ["#!{Type}T", "N", "C", "C", "E", "E", "E", "E"]
    assert [line[0] for line in lines[2:]] == [group for group, _ in GROUPS], "no header rows"
    assert lines[2][1:8] == ["0.001", "", "", "5000.0", "10000.0", "15000.0", "20000.0"]
    assert not sidecar_path(tmp_path / "perseus.txt").exists(), "a table has no sidecar"


def test_the_sample_columns_are_named_by_the_abundance_kind(tmp_path: Path) -> None:
    lines = _write(_maxquant(tmp_path), tmp_path / "perseus.txt", abundance="LFQ_Intensity")

    assert lines[0][4:] == [f"LFQ intensity {name}" for name in EXPERIMENTS]
    assert lines[2][4:] == ["2500.0", "5000.0", "7500.0", "10000.0"]


def test_an_unlisted_kind_keeps_the_entry_name_and_missing_values_are_nan(
    diann: DiannInput, tmp_path: Path
) -> None:
    """DIA-NN 1.8.1's MaxLFQ is catalogued as an unknown kind."""
    lines = _write(_diann(diann), tmp_path / "perseus.txt")

    assert lines[0][4:] == [f"Intensity {run}" for run in RUNS]
    assert "NaN" in {value for line in lines[2:] for value in line[4:]}


def test_annotation_columns_become_categorical_rows(diann: DiannInput, tmp_path: Path) -> None:
    annotated = AnnotationCompiler().compile(write_design(tmp_path)).parse(_diann(diann))
    lines = _write(annotated.annotate().parsed, tmp_path / "perseus.txt")
    empty = [""] * (len(VAR_COLUMNS) - 1)

    assert lines[2] == ["#!{C:group}", *empty, *(run[4] for run in RUNS)]
    assert lines[3] == ["#!{C:BioReplicate}", *empty, *(run[4:] for run in RUNS)]


def test_decoys_and_contaminants_are_marked_plus(diann: DiannInput, tmp_path: Path) -> None:
    parsed = _diann(diann)
    level = parsed.levels["protein"]
    groups = level.var.frame.get_column("Protein_Group")
    level.var.frame = level.var.frame.with_columns(
        (groups == "P1").alias("apb_Decoy"), (groups == "P7").alias("apb_Contaminant")
    )
    lines = _write(parsed, tmp_path / "perseus.txt")
    marks = {line[0]: (line[2], line[3]) for line in lines[2:]}

    assert marks["P1"] == ("+", "")
    assert marks["P7"] == ("", "+")
    assert marks["P3"] == ("", "")


def test_a_perseus_table_of_several_levels_is_refused(diann: DiannInput, tmp_path: Path) -> None:
    exporter = Exporter("perseus")
    exported = exporter.export(_diann(diann))
    (level,) = exported.levels.values()
    two: ParsedLevels = replace(exported, levels={"a": level, "b": level})

    with pytest.raises(ValueError, match="holds one level"):
        exporter.write(two, tmp_path / "perseus.txt")


def test_the_table_reads_back_with_its_types(tmp_path: Path) -> None:
    path = tmp_path / "perseus.txt"
    _write(_maxquant(tmp_path), path)
    table = pl.read_csv(path, separator="\t", skip_rows_after_header=1)

    assert table.height == len(GROUPS)
    assert table.get_column("Intensity S1").dtype == pl.Float64


def test_a_name_used_twice_is_refused(diann: DiannInput) -> None:
    (level,) = Exporter("perseus").export(_diann(diann)).levels.values()
    level.var.frame = level.var.frame.with_columns(pl.lit(1.0).alias("Intensity run_A1"))

    with pytest.raises(ValueError, match=r"would name \['Intensity run_A1'\] twice"):
        perseus_table(level)
