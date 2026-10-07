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
    "PG.Q.Value",
    "Lib.PG.Q.Value",
    "Global.PG.Q.Value",
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
                0.0001 * (index + 1) * (run_index + 1),
                library_q * (index + 1) / 2,
                global_q * (index + 1) / 2,
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


EXPERIMENTS = ("S1", "S2", "S3", "S4")
# Sequence, proteins and razor protein of each peptide. The razor protein names a one-member
# group whole, leads a two-member group, or leads a group with a UniProt "db|ACC|NAME" member.
# Every razor protein is in a group, as in real MaxQuant output, which APB2 requires.
PEPTIDES = (
    ("PEPTIDEK", "P1", "P1"),
    ("SECONDK", "P2;P3", "P2"),
    ("THIRDPEPK", "P2;P3", "P2"),
    ("FOURTHK", "sp|P4|P4_HUMAN;P5", "sp|P4|P4_HUMAN"),
    ("FIFTHK", "P9", "P9"),
)
# Each protein group and the peptides it holds.
GROUPS = (("P1", (0,)), ("P2;P3", (1, 2)), ("sp|P4|P4_HUMAN;P5", (3,)), ("P9", (4,)))
# The protein-group id of each peptide, as evidence.txt links it.
GROUP_OF = {member: group for group, (_, members) in enumerate(GROUPS) for member in members}
MQPAR = """<?xml version="1.0" encoding="utf-8"?>
<MaxQuantParams>
  <maxQuantVersion>2.4.0.0</maxQuantVersion>
  <peptideFdr>0.01</peptideFdr>
  <proteinFdr>0.01</proteinFdr>
  <matchBetweenRuns>False</matchBetweenRuns>
  <minPepLen>7</minPepLen>
  <parameterGroups>
    <parameterGroup>
      <mainSearchTol>4.5</mainSearchTol>
      <enzymeMode>0</enzymeMode>
      <enzymes><string>Trypsin/P</string></enzymes>
      <maxMissedCleavages>2</maxMissedCleavages>
      <fixedModifications><string>Carbamidomethyl (C)</string></fixedModifications>
      <variableModifications><string>Oxidation (M)</string></variableModifications>
      <maxNmods>5</maxNmods>
      <maxCharge>7</maxCharge>
    </parameterGroup>
  </parameterGroups>
  <msmsParamsArray>
    <msmsParams Name="FTMS" MatchTolerance="20" MatchToleranceInPpm="True" />
    <msmsParams Name="ITMS" MatchTolerance="0.5" MatchToleranceInPpm="False" />
  </msmsParamsArray>
</MaxQuantParams>
"""


EVIDENCE_COLUMNS = (
    "Raw file",
    "Experiment",
    "Sequence",
    "Modified sequence",
    "Modifications",
    "Charge",
    "Mass",
    "Proteins",
    "Leading razor protein",
    "Intensity",
    "MS/MS count",
    "Retention time",
    "Score",
    "PEP",
    "id",
    "Protein group IDs",
)
# The per-experiment intensity columns sit between the first four columns and the IDs.
PEPTIDE_COLUMNS = (
    "Sequence",
    "Proteins",
    "Leading razor protein",
    "PEP",
    "id",
    "Mod. peptide IDs",
    "Evidence IDs",
)
GROUP_COLUMNS = ("Protein IDs", "Q-value", "id", "Peptide IDs", "Evidence IDs")


@dataclass(frozen=True, slots=True)
class MaxQuantInput:
    """A synthetic MaxQuant txt folder and its mqpar.xml."""

    folder: Path
    params: Path


def quantified(experiment: int, peptide: int) -> bool:
    """Whether a peptide has an intensity in an experiment; one in five cells is absent."""
    return (experiment + peptide) % 5 != 4


def _table(path: Path, header: list[str], rows: list[list[object]]) -> None:
    lines = ["\t".join(header), *("\t".join(str(value) for value in row) for row in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_maxquant(folder: Path, /, *, peptides: bool = True) -> MaxQuantInput:
    """Write MaxQuant 2.4 evidence, peptide and protein-group tables and an mqpar.xml.

    ``peptides=False`` writes evidence.txt alone, as ProteoBench submissions carry it. The
    q-value of group ``g`` is ``0.001 * (g + 1)``; the PEP of peptide ``p`` is ``0.01 * (p + 1)``.
    """
    evidence: list[list[object]] = []
    for experiment, name in enumerate(EXPERIMENTS):
        for index, (sequence, proteins, razor) in enumerate(PEPTIDES):
            if quantified(experiment, index):
                intensity = 1000.0 * (experiment + 1) * (index + 1)
                evidence.append(
                    [
                        *(f"raw_{name}", name, sequence, f"_{sequence}_", "Unmodified", 2),
                        *(1000.5, proteins, razor, intensity, 1, 10.0, 50.0, 0.01 * (index + 1)),
                        len(evidence),
                        GROUP_OF[index],
                    ]
                )
    _table(folder / "evidence.txt", list(EVIDENCE_COLUMNS), evidence)
    if peptides:
        intensities = [f"Intensity {name}" for name in EXPERIMENTS]
        peptide_rows: list[list[object]] = [
            [
                *(sequence, proteins, razor, 0.01 * (index + 1)),
                *(
                    1000.0 * (experiment + 1) * (index + 1) if quantified(experiment, index) else 0
                    for experiment in range(len(EXPERIMENTS))
                ),
                *(index, index, index),
            ]
            for index, (sequence, proteins, razor) in enumerate(PEPTIDES)
        ]
        _table(
            folder / "peptides.txt",
            [*PEPTIDE_COLUMNS[:4], *intensities, *PEPTIDE_COLUMNS[4:]],
            peptide_rows,
        )
        group_rows: list[list[object]] = [
            [
                *(group, 0.001 * (index + 1)),
                *(5000.0 * (experiment + 1) for experiment in range(len(EXPERIMENTS))),
                *(index, ";".join(map(str, members)), ";".join(map(str, members))),
            ]
            for index, (group, members) in enumerate(GROUPS)
        ]
        _table(
            folder / "proteinGroups.txt",
            [*GROUP_COLUMNS[:2], *intensities, *GROUP_COLUMNS[2:]],
            group_rows,
        )
    params = folder / "mqpar.xml"
    params.write_text(MQPAR, encoding="utf-8")
    return MaxQuantInput(folder=folder, params=params)
