"""Synthetic DIA-NN input: small, fully known, no real data."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import anndata as ad
import mudata as md
import pandas as pd
import pytest
from apb2.api import ParsedLevels, UnsJsonCodec, read_parsed_levels, write_parsed_levels
from polars.testing import assert_frame_equal

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


def write_design(folder: Path) -> Path:
    """Write a prolfquapp-style sample table: each run's ``group`` (A or B) and ``BioReplicate``,
    the columns apb-catalog reads as condition and biological replicate."""
    path = folder / "design.tsv"
    _table(path, ["raw_file", "group", "BioReplicate"], [[run, run[4], run[4:]] for run in RUNS])
    return path


# Each ion of the FragPipe example and how it was found in samples A and B: a match type is a
# categorical layer, stored as codes, whose labels the FragPipe parsing rule declares.
FRAGPIPE_IONS = (("PEPTIDEK", "P1", ("MS/MS", "MBR")), ("SECONDK", "P2", ("MS/MS", "unmatched")))
FRAGPIPE_SAMPLES = ("A", "B")
FRAGPIPE_PER_SAMPLE = (
    "Spectral Count",
    "Apex Retention Time",
    "Retention Time Start",
    "Retention Time End",
    "Retention Time FWHM",
    "Intensity",
    "Match Type",
)


def write_fragpipe(folder: Path) -> Path:
    """Write a FragPipe 22 ``combined_ion.tsv``; an unmatched ion's intensity is 0."""
    fixed = ["Peptide Sequence", "Modified Sequence", "Charge", "M/Z", "Protein"]
    fixed += ["Mapped Proteins", "Protein ID", "Gene", "Assigned Modifications"]
    header = fixed + [
        f"{sample} {field}" for field in FRAGPIPE_PER_SAMPLE for sample in FRAGPIPE_SAMPLES
    ]
    rows: list[list[object]] = []
    for index, (sequence, protein, matches) in enumerate(FRAGPIPE_IONS):
        row: list[object] = [sequence, sequence, 2, 465.7, f"sp|{protein}|{protein}_HUMAN", ""]
        row += [protein, f"G{index}", ""]
        for field in FRAGPIPE_PER_SAMPLE:
            for sample, match in enumerate(matches):
                if field == "Match Type":
                    row.append(match)
                elif field == "Intensity":
                    row.append(0 if match == "unmatched" else 1000 * (index + 1) + sample)
                else:
                    row.append(1 + index + sample)
        rows.append(row)
    path = folder / "combined_ion.tsv"
    _table(path, header, rows)
    return path


def apb_part(namespace: Mapping[str, Any]) -> dict[str, Any]:
    """One stored APB part, decoded."""
    codec = UnsJsonCodec()
    return cast(dict[str, Any], codec.decode(namespace, codec.storage(namespace)))


def level_part(adata: ad.AnnData | md.MuData) -> dict[str, Any]:
    """An exported AnnData's own APB part: ``uns[<name>]["apb"]`` standalone, else ``uns["apb"]``."""
    owners = [
        value["apb"] for value in adata.uns.values() if isinstance(value, dict) and "apb" in value
    ]
    return apb_part(owners[0] if owners else adata.uns["apb"])


def export_sources(adata: ad.AnnData | md.MuData) -> dict[str, Any]:
    """The APB source of every entry the export wrote into one AnnData."""
    return cast(dict[str, Any], level_part(adata)["export"]["provenance"]["sources"])


def column(frame: pd.DataFrame, name: str) -> pd.Series:
    """One column; pandas types a column lookup as Series or DataFrame."""
    values = frame[name]
    assert isinstance(values, pd.Series)
    return values


def read_back(exported: ParsedLevels, path: Path) -> ad.AnnData | md.MuData:
    """An export as its tool opens it: written as ``Exporter.write`` writes it, then read."""
    write_parsed_levels(exported, path, var_name_prefix=False)
    return md.read_h5mu(path) if path.suffix == ".h5mu" else ad.read_h5ad(path)


def assert_rewrite_identical(written: Path, again: Path) -> None:
    """An export read by APB2 and written again holds the same result, APB2's bar for its own."""
    first = read_parsed_levels(written)
    write_parsed_levels(first, again, var_name_prefix=False)
    second = read_parsed_levels(again)
    assert (second.uns, second.metadata, second.export_uns) == (
        first.uns,
        first.metadata,
        first.export_uns,
    )
    assert list(second.feature_relations) == list(first.feature_relations)
    for name, relation in first.feature_relations.items():
        assert second.feature_relations[name].metadata == relation.metadata
        assert_frame_equal(second.feature_relations[name].coordinates, relation.coordinates)
    assert list(second.levels) == list(first.levels)
    for name, level in first.levels.items():
        other = second.levels[name]
        assert (other.uns, other.metadata, other.export_uns) == (
            level.uns,
            level.metadata,
            level.export_uns,
        )
        assert_frame_equal(other.obs.frame, level.obs.frame)
        assert_frame_equal(other.var.frame, level.var.frame)
        assert list(other.layers) == list(level.layers)
        for layer, table in level.layers.items():
            assert_frame_equal(other.layers[layer].wide(), table.wide())
