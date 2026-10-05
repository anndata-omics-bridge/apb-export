"""Synthetic DIA-NN input: small, fully known, no real data."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

COLUMNS = (
    "Run",
    "Modified.Sequence",
    "Stripped.Sequence",
    "Precursor.Charge",
    "Precursor.Id",
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "Precursor.Quantity",
    "Precursor.Normalised",
    "Q.Value",
    "Global.Q.Value",
    "Lib.Q.Value",
    "PEP",
    "RT",
    "PG.MaxLFQ",
)
RUNS = ("run_A1", "run_A2", "run_B1", "run_B2")
# Modified sequence and protein group. Covers a fixed modification, a contaminant, a
# UniProt "db|ACC|NAME" member, a marker-prefixed contaminant and a two-member group.
PRECURSORS = (
    ("PEPTIDEK", "P1"),
    ("AAC(UniMod:4)LLK", "sp|P2|P2_HUMAN"),
    ("YEASTPEPK", "P3"),
    ("CONTPEPK", "Cont_P4"),
    ("SHAREDK", "P5;CON__P6"),
    ("LASTPEPK", "P7"),
)


@dataclass(frozen=True, slots=True)
class DiannInput:
    """A synthetic DIA-NN report, its log, and what was written into it."""

    report: Path
    log: Path
    cells: int


def observed(run: int, precursor: int) -> bool:
    """Whether a precursor is reported in a run; one in five cells is absent."""
    return (run + precursor) % 5 != 4


def write_diann(folder: Path, /, *, library_q: float = 0.0, global_q: float = 0.002) -> DiannInput:
    """Write a DIA-NN 1.8.1 report and log; ``library_q`` > 0 mimics match-between-runs."""
    rows = ["\t".join(COLUMNS)]
    cells = 0
    for run_index, run in enumerate(RUNS):
        for index, (modified, group) in enumerate(PRECURSORS):
            if not observed(run_index, index):
                continue
            cells += 1
            stripped = modified.replace("(UniMod:4)", "")
            scale = (run_index + 1) * (index + 1)
            values = (
                run,
                modified,
                stripped,
                2,
                f"{modified}2",
                group,
                group,
                group,
                f"G{index}",
                1000.0 * scale,
                900.0 * scale,
                0.001 * (index + 1),
                global_q * (index + 1),
                library_q * (index + 1),
                0.01 * (index + 1),
                10.0 + index,
                500.0 * (index + 1),
            )
            rows.append("\t".join(str(value) for value in values))
    report = folder / "report.tsv"
    report.write_text("\n".join(rows) + "\n", encoding="utf-8")
    log = folder / "report.log.txt"
    log.write_text(
        "DIA-NN 1.8.1 (Data-Independent Acquisition by Neural Networks)\ndiann --unimod4\n",
        encoding="utf-8",
    )
    return DiannInput(report=report, log=log, cells=cells)


@pytest.fixture
def diann(tmp_path: Path) -> DiannInput:
    """A synthetic DIA-NN run without match-between-runs."""
    return write_diann(tmp_path)
