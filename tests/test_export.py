"""msmu's experimental wide export: ions by runs, and the vendor's confidence by meaning."""

from __future__ import annotations

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
from conftest import (
    PRECURSORS,
    RUNS,
    DiannInput,
    apb_part,
    assert_rewrite_identical,
    column,
    export_sources,
    level_part,
    write_diann,
)

MSMU_DIANN_COLUMNS = [
    "precursor",
    "proteins",
    "peptide",
    "stripped_peptide",
    "charge",
    "peptide_length",
    "decoy",
    "contaminant",
    "PEP",
    "q_value",
]


def _parsed(diann: DiannInput) -> ParsedLevels:
    return ParseRuleCompiler.detect(diann.report, diann.log, ("peptidoform_ion",)).compile().parse()


def _psm(exporter: Exporter, parsed: ParsedLevels, folder: Path) -> ad.AnnData:
    """msmu's psm modality as msmu reads it from the file the exporter writes."""
    target = folder / "msmu.h5mu"
    exporter.write(exporter.export(parsed), target)
    psm = md.read_h5mu(target)["psm"]
    assert isinstance(psm, ad.AnnData)
    return psm


def _var(psm: object) -> pd.DataFrame:
    """The psm ``var`` frame; MuData and AnnData type modalities and frames as unions."""
    assert isinstance(psm, ad.AnnData)
    assert isinstance(psm.var, pd.DataFrame)
    return psm.var


def _x(psm: ad.AnnData, layer: str | None = None) -> np.ndarray:
    matrix = psm.X if layer is None else psm.layers[layer]
    assert isinstance(matrix, np.ndarray), "wide output is dense"
    return matrix


def test_one_feature_per_precursor(diann: DiannInput) -> None:
    psm = _psm(
        Exporter("msmu", abundance="Precursor_Quantity"), _parsed(diann), diann.report.parent
    )

    assert list(psm.obs_names) == sorted(RUNS)
    assert psm.n_vars == len(PRECURSORS)
    assert list(psm.var.columns) == MSMU_DIANN_COLUMNS
    x = _x(psm)
    assert x.dtype == np.float32
    assert np.isnan(x).sum() == len(RUNS) * len(PRECURSORS) - diann.cells, "absent cells are NaN"
    first = _var(psm).iloc[0]
    assert psm.var_names[0] == f"{first['peptide']}/{first['charge']}"


def test_quantities_and_confidence_come_from_the_report(diann: DiannInput) -> None:
    psm = _psm(
        Exporter("msmu", abundance="Precursor_Quantity"), _parsed(diann), diann.report.parent
    )
    var = _var(psm)
    column = list(psm.var_names).index("YEASTPEPK/2")
    run = list(psm.obs_names).index("run_A2")

    # run_A2 is the second run and YEASTPEPK the third precursor: scale 2 * 3.
    assert _x(psm)[run, column] == pytest.approx(6000.0)
    assert var.loc["YEASTPEPK/2", "q_value"] == pytest.approx(0.006), "Global.Q.Value"
    assert var.loc["YEASTPEPK/2", "PEP"] == pytest.approx(0.03)
    assert _x(psm, "q_value")[run, column] == pytest.approx(0.006), "the run's value as a layer"
    assert _x(psm, "PEP")[run, column] == pytest.approx(0.03)
    assert export_sources(psm)["var.q_value"]["kind"] == "global_q_value"
    assert export_sources(psm)["var.q_value"]["how"] == "layer_minimum"


def test_library_q_value_wins_when_match_between_runs_filled_it(tmp_path: Path) -> None:
    diann = write_diann(tmp_path, library_q=0.0005)
    psm = _psm(Exporter("msmu"), _parsed(diann), diann.report.parent)

    assert export_sources(psm)["var.q_value"]["kind"] == "library_q_value"
    assert _var(psm).loc["PEPTIDEK/2", "q_value"] == pytest.approx(0.0005)


def test_protein_groups_are_written_as_msmu_writes_them(diann: DiannInput) -> None:
    var = _var(_psm(Exporter("msmu"), _parsed(diann), diann.report.parent))

    assert var.loc["AAC[UNIMOD:4]LLK/2", "proteins"] == "P2"
    assert var.loc["CONTPEPK/2", "proteins"] == "Cont_P4"
    assert var.loc["SHAREDK/2", "proteins"] == "P5;Cont_P6"
    # DIA-NN marks nothing and nothing was FASTA-checked, so nothing is a contaminant or decoy.
    assert (column(var, "contaminant") == 0).all()
    assert (column(var, "decoy") == 0).all()


def test_a_feature_without_proteins_is_left_out_and_counted(diann: DiannInput) -> None:
    """PEAKS leaves some precursors without an accession; msmu would group them as ""."""
    parsed = _parsed(diann)
    level = parsed.levels["peptidoform_ion"]
    accessions = level.var.roles["fasta_accessions"]
    unassigned = level.var.frame.get_column("ProForma_ion") == "PEPTIDEK/2"
    level.var.frame = level.var.frame.with_columns(
        pl.when(unassigned).then(pl.lit("")).otherwise(pl.col(accessions)).alias(accessions)
    )
    full = _psm(Exporter("msmu"), _parsed(diann), diann.report.parent)

    psm = _psm(Exporter("msmu"), parsed, diann.report.parent)

    kept = [index for index, name in enumerate(full.var_names) if name != "PEPTIDEK/2"]
    assert list(psm.var_names) == [full.var_names[index] for index in kept]
    assert (column(_var(psm), "proteins") != "").all()
    np.testing.assert_array_equal(_x(psm), _x(full)[:, kept])
    assert export_sources(psm)["var.proteins"]["dropped"] == "1"


def test_contaminants_merge_the_vendor_marking_and_the_fasta_match(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    level = parsed.levels["peptidoform_ion"]
    ions = level.var.frame.get_column("ProForma_ion")
    level.var.frame = level.var.frame.with_columns(
        (ions == "SHAREDK/2").alias("apb_Contaminant"),
        (ions == "LASTPEPK/2").alias("apb_Decoy"),
    )
    level.varm["fasta_validation"] = pl.DataFrame(
        {"fasta_matches_contaminant": (ions == "CONTPEPK/2").to_list()}
    )

    var = _var(_psm(Exporter("msmu"), parsed, diann.report.parent))

    assert var.loc["CONTPEPK/2", "contaminant"] == 1
    assert var.loc["SHAREDK/2", "contaminant"] == 1
    assert var.loc["PEPTIDEK/2", "contaminant"] == 0
    assert var.loc["LASTPEPK/2", "decoy"] == 1
    assert var.loc["PEPTIDEK/2", "decoy"] == 0


def test_written_file_reads_back(diann: DiannInput, tmp_path: Path) -> None:
    exporter = Exporter("msmu")
    target = tmp_path / "out.h5mu"
    exporter.write(exporter.export(_parsed(diann)), target)
    back = md.read_h5mu(target)

    psm = back["psm"]
    assert isinstance(psm, ad.AnnData)
    assert psm.shape == (len(RUNS), len(PRECURSORS))
    root = apb_part(back.uns["apb"])["export"]
    assert root["schema_version"] == "1"
    assert root["provenance"]["package"] == "apb-export"
    assert root["provenance"]["export_rule"] == "msmu/v0_4/rules.json"
    record = level_part(psm)["export"]
    assert [entry["name"] for entry in record["summary"]] == ["exported_matrices", "absent_entries"]
    assert record["provenance"]["source_level"] == "peptidoform_ion"
    assert "parse" in level_part(psm) and "parse" in apb_part(back.uns["apb"])
    assert psm.uns["search_engine"] == "dia-nn"
    assert psm.uns["apb_export_status"] == "experimental", "msmu's wide file is experimental"
    assert list(psm.obs.columns) == ["filename"], "the run key stays a column and names the index"
    assert psm.obs.index.name == "filename"
    assert psm.var.index.name == "precursor"
    assert _x(psm).dtype == np.float32
    assert {"PEP", "q_value"} <= set(psm.layers.keys()), "each run's value as a layer"
    assert _var(psm)["PEP"].dtype == np.float64, "APB2 writes number columns as Float64"
    assert_rewrite_identical(target, tmp_path / "again.h5mu")


def test_refusals_name_the_problem(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    with pytest.raises(ValueError, match="level has no layer 'nope'"):
        Exporter("msmu", abundance="nope").export(parsed)
    with pytest.raises(ValueError, match="'Q_Value' does not carry the abundance role"):
        Exporter("msmu", abundance="Q_Value").export(parsed)

    without_ion = replace(parsed, levels={})
    with pytest.raises(ValueError, match="needs APB2's peptidoform_ion level"):
        Exporter("msmu").export(without_ion)

    level = parsed.levels["peptidoform_ion"]
    plexed = replace(level, obs=replace(level.obs, key_columns=("Run", "Channel")))
    with pytest.raises(ValueError, match="one observation per run"):
        Exporter("msmu").export(replace(parsed, levels={"peptidoform_ion": plexed}))
