"""The export has msmu's layout and carries the vendor's confidence by meaning."""

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
from scipy import sparse

from apb_export.api import Exporter
from conftest import RUNS, DiannInput, apb_part, export_sources, level_part, write_diann

MSMU_DIANN_COLUMNS = [
    "proteins",
    "peptide",
    "stripped_peptide",
    "filename",
    "charge",
    "peptide_length",
    "decoy",
    "contaminant",
    "PEP",
    "q_value",
]


def _parsed(diann: DiannInput) -> ParsedLevels:
    return ParseRuleCompiler(diann.report, diann.log, requested_levels=("ion",)).compile().parse()


def _var(psm: object) -> pd.DataFrame:
    """The psm ``var`` frame; MuData and AnnData type modalities and frames as unions."""
    assert isinstance(psm, ad.AnnData)
    assert isinstance(psm.var, pd.DataFrame)
    return psm.var


def _search_result(psm: object) -> pd.DataFrame:
    assert isinstance(psm, ad.AnnData)
    table = psm.varm["search_result"]
    assert isinstance(table, pd.DataFrame)
    return table


def _matrix(psm: object) -> sparse.csr_matrix:
    assert isinstance(psm, ad.AnnData)
    assert sparse.issparse(psm.X)
    return sparse.csr_matrix(psm.X)


def test_one_feature_per_observed_cell(diann: DiannInput) -> None:
    psm = Exporter("msmu", abundance="Precursor_Quantity").export(_parsed(diann))["psm"]

    assert list(psm.obs_names) == sorted(RUNS)
    assert psm.n_vars == diann.cells
    assert list(psm.var.columns) == MSMU_DIANN_COLUMNS
    matrix = _matrix(psm)
    assert matrix.dtype == np.float32
    assert (np.diff(sparse.csc_matrix(matrix).indptr) == 1).all(), "one run per feature"
    first = _var(psm).iloc[0]
    assert psm.var_names[0] == f"{first['filename']}.{first['peptide']}/{first['charge']}"


def test_quantities_and_confidence_come_from_the_report(diann: DiannInput) -> None:
    psm = Exporter("msmu", abundance="Precursor_Quantity").export(_parsed(diann))["psm"]
    var = _var(psm)
    row = var.index.get_loc("run_A2.YEASTPEPK/2")
    run = list(psm.obs_names).index("run_A2")

    # run_A2 is the second run and YEASTPEPK the third precursor: scale 2 * 3.
    assert _matrix(psm)[run, row] == pytest.approx(6000.0)
    assert var.iloc[row]["q_value"] == pytest.approx(0.006), "Global.Q.Value, not Q.Value"
    assert var.iloc[row]["PEP"] == pytest.approx(0.03)
    assert export_sources(psm)["var.q_value"]["kind"] == "global_q_value"


def test_library_q_value_wins_when_match_between_runs_filled_it(tmp_path: Path) -> None:
    diann = write_diann(tmp_path, library_q=0.0005)
    psm = Exporter("msmu").export(_parsed(diann))["psm"]

    assert export_sources(psm)["var.q_value"]["kind"] == "library_q_value"
    assert _var(psm).loc["run_A1.PEPTIDEK/2", "q_value"] == pytest.approx(0.0005)


def test_protein_groups_are_written_as_msmu_writes_them(diann: DiannInput) -> None:
    var = _var(Exporter("msmu").export(_parsed(diann))["psm"])

    assert var.loc["run_A1.AAC[UNIMOD:4]LLK/2", "proteins"] == "P2"
    assert var.loc["run_A1.CONTPEPK/2", "proteins"] == "Cont_P4"
    # SHAREDK is absent from run_A1 by design (see conftest.observed).
    assert var.loc["run_A2.SHAREDK/2", "proteins"] == "P5;Cont_P6"
    # DIA-NN marks nothing and nothing was FASTA-checked, so nothing is a contaminant or decoy.
    assert (var["contaminant"] == 0).all()
    assert (var["decoy"] == 0).all()


def test_contaminants_merge_the_vendor_marking_and_the_fasta_match(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    level = parsed.levels["ion"]
    ions = level.var.frame.get_column("ProForma_ion")
    level.var.frame = level.var.frame.with_columns(
        (ions == "SHAREDK/2").alias("apb_Contaminant"),
        (ions == "LASTPEPK/2").alias("apb_Decoy"),
    )
    level.varm["fasta_validation"] = pl.DataFrame(
        {"fasta_matches_contaminant": (ions == "CONTPEPK/2").to_list()}
    )

    var = _var(Exporter("msmu").export(parsed)["psm"])

    assert var.loc["run_A1.CONTPEPK/2", "contaminant"] == 1
    assert var.loc["run_A2.SHAREDK/2", "contaminant"] == 1
    assert var.loc["run_A1.PEPTIDEK/2", "contaminant"] == 0
    assert var.loc["run_A1.LASTPEPK/2", "decoy"] == 1
    assert var.loc["run_A1.PEPTIDEK/2", "decoy"] == 0


def test_every_source_value_is_kept_in_search_result(diann: DiannInput) -> None:
    psm = Exporter("msmu").export(_parsed(diann))["psm"]
    source = _search_result(psm)

    assert list(source.index) == list(psm.var_names)
    for name in ("Precursor_Quantity", "Q_Value", "Global_Q_Value", "RT", "Protein_Group"):
        assert name in source.columns
    assert source.loc["run_A1.PEPTIDEK/2", "RT"] == pytest.approx(10.0)


def test_written_file_reads_back(diann: DiannInput, tmp_path: Path) -> None:
    mdata = Exporter("msmu").export(_parsed(diann))
    assert isinstance(mdata, md.MuData)
    target = tmp_path / "out.h5mu"
    mdata.write_h5mu(target)
    back = md.read_h5mu(target)

    psm = back["psm"]
    assert psm.shape == mdata["psm"].shape
    root = apb_part(back.uns["apb"])["export"]
    assert root["schema_version"] == "1"
    assert root["provenance"]["package"] == "apb-export"
    assert root["provenance"]["export_rule"] == "msmu/v0_4/rules.json"
    record = level_part(psm)["export"]
    assert [entry["name"] for entry in record["summary"]] == ["exported_matrices", "absent_entries"]
    assert record["provenance"]["source_level"] == "ion"
    assert "parse" in level_part(psm) and "parse" in apb_part(back.uns["apb"])
    assert psm.uns["search_engine"] == "dia-nn"
    assert list(psm.obs.columns) == [], "runs are the index only"
    assert "feature_id" not in _var(psm).columns, "feature ids are the index only"
    assert _matrix(psm).dtype == np.float32
    assert _var(psm)["PEP"].dtype == np.float32
    assert _var(psm)["q_value"].dtype == np.float32


def test_refusals_name_the_problem(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    with pytest.raises(ValueError, match="level has no layer 'nope'"):
        Exporter("msmu", abundance="nope").export(parsed)
    with pytest.raises(ValueError, match="'Q_Value' does not carry the abundance role"):
        Exporter("msmu", abundance="Q_Value").export(parsed)

    without_ion = replace(parsed, levels={})
    with pytest.raises(ValueError, match="needs APB2's ion level"):
        Exporter("msmu").export(without_ion)

    level = parsed.levels["ion"]
    plexed = replace(level, obs=replace(level.obs, key_columns=("Run", "Channel")))
    with pytest.raises(ValueError, match="one observation per run"):
        Exporter("msmu").export(replace(parsed, levels={"ion": plexed}))
