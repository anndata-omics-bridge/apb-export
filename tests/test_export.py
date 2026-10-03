"""The export has msmu's layout and carries the vendor's confidence by meaning."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from apb2.api import ParsedLevels, ParseRuleCompiler
from scipy import sparse

from apb_msmu.api import MsmuExporter
from apb_msmu.confidence import Confidence, resolve_confidence
from apb_msmu.psm import psm_rows
from conftest import RUNS, DiannInput, write_diann

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
    psm = MsmuExporter(abundance="Precursor_Quantity").export(_parsed(diann))["psm"]

    assert list(psm.obs_names) == sorted(RUNS)
    assert psm.n_vars == diann.cells
    assert list(psm.var.columns) == MSMU_DIANN_COLUMNS
    matrix = _matrix(psm)
    assert matrix.dtype == np.float32
    assert (np.diff(sparse.csc_matrix(matrix).indptr) == 1).all(), "one run per feature"
    first = _var(psm).iloc[0]
    assert psm.var_names[0] == f"{first['filename']}.{first['peptide']}/{first['charge']}"


def test_quantities_and_confidence_come_from_the_report(diann: DiannInput) -> None:
    psm = MsmuExporter(abundance="Precursor_Quantity").export(_parsed(diann))["psm"]
    var = _var(psm)
    row = var.index.get_loc("run_A2.YEASTPEPK/2")
    run = list(psm.obs_names).index("run_A2")

    # run_A2 is the second run and YEASTPEPK the third precursor: scale 2 * 3.
    assert _matrix(psm)[run, row] == pytest.approx(6000.0)
    assert var.iloc[row]["q_value"] == pytest.approx(0.006), "Global.Q.Value, not Q.Value"
    assert var.iloc[row]["PEP"] == pytest.approx(0.03)
    assert psm.uns["apb"]["q_value_kind"] == "global_q_value"


def test_library_q_value_wins_when_match_between_runs_filled_it(tmp_path: Path) -> None:
    diann = write_diann(tmp_path, library_q=0.0005)
    psm = MsmuExporter().export(_parsed(diann))["psm"]

    assert psm.uns["apb"]["q_value_kind"] == "library_q_value"
    assert _var(psm).loc["run_A1.PEPTIDEK/2", "q_value"] == pytest.approx(0.0005)


def test_per_run_q_value_when_nothing_better_is_reported(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    confidence = resolve_confidence(parsed)
    run_only = replace(
        confidence, q_value=parsed.levels["ion"].layers["Q_Value"], q_value_kind="q_value"
    )
    rows = psm_rows(parsed.levels["ion"], abundance="Precursor_Quantity", confidence=run_only)

    assert rows.var.filter(rows.var["feature_id"] == "run_A1.PEPTIDEK/2")["q_value"].item() == (
        pytest.approx(0.001)
    )


def test_unreported_pep_is_nan_and_unreported_q_value_absent(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    nothing = Confidence(q_value=None, q_value_kind=None, pep=None)
    rows = psm_rows(parsed.levels["ion"], abundance="Precursor_Quantity", confidence=nothing)

    assert "q_value" not in rows.var.columns, "msmu's filter must not see an invented q-value"
    assert all(math.isnan(value) for value in rows.var["PEP"].to_list()), (
        "msmu's to_peptide needs the PEP column; NaN says not reported"
    )


def test_protein_groups_are_written_as_msmu_writes_them(diann: DiannInput) -> None:
    var = _var(MsmuExporter().export(_parsed(diann))["psm"])

    assert var.loc["run_A1.AAC[UNIMOD:4]LLK/2", "proteins"] == "P2"
    assert var.loc["run_A1.CONTPEPK/2", "proteins"] == "Cont_P4"
    assert var.loc["run_A1.CONTPEPK/2", "contaminant"] == 1
    # SHAREDK is absent from run_A1 by design (see conftest.observed).
    assert var.loc["run_A2.SHAREDK/2", "proteins"] == "P5;Cont_P6"
    assert var.loc["run_A2.SHAREDK/2", "contaminant"] == 1
    assert var.loc["run_A1.PEPTIDEK/2", "contaminant"] == 0


def test_every_source_value_is_kept_in_search_result(diann: DiannInput) -> None:
    psm = MsmuExporter().export(_parsed(diann))["psm"]
    source = _search_result(psm)

    assert list(source.index) == list(psm.var_names)
    for name in ("Precursor_Quantity", "Q_Value", "Global_Q_Value", "RT", "Protein_Group"):
        assert name in source.columns
    assert source.loc["run_A1.PEPTIDEK/2", "RT"] == pytest.approx(10.0)


def test_written_file_reads_back(diann: DiannInput, tmp_path: Path) -> None:
    mdata = MsmuExporter().export(_parsed(diann))
    target = tmp_path / "out.h5mu"
    mdata.write_h5mu(target)
    back = md.read_h5mu(target)

    assert back["psm"].shape == mdata["psm"].shape
    assert back["psm"].uns["apb"]["package"] == "apb-msmu"
    assert back["psm"].uns["search_engine"] == "dia-nn"


def test_refusals_name_the_problem(diann: DiannInput) -> None:
    parsed = _parsed(diann)
    with pytest.raises(ValueError, match="level has no layer 'nope'"):
        MsmuExporter(abundance="nope").export(parsed)
    with pytest.raises(ValueError, match="'Q_Value' does not carry the abundance role"):
        MsmuExporter(abundance="Q_Value").export(parsed)

    without_ion = replace(parsed, levels={})
    with pytest.raises(ValueError, match="needs APB2's ion level"):
        MsmuExporter().export(without_ion)

    level = parsed.levels["ion"]
    plexed = replace(level, obs=replace(level.obs, key_columns=("Run", "Channel")))
    with pytest.raises(ValueError, match="one observation per run"):
        MsmuExporter().export(replace(parsed, levels={"ion": plexed}))
