"""msmu's own DIA-NN reader and an apb-msmu export give the same msmu results.

Runs only with the ``parity`` dependency group, which installs msmu.
"""

from __future__ import annotations

from pathlib import Path

import mudata as md
import numpy as np
import pandas as pd
import pytest
from apb2.api import ParseRuleCompiler
from scipy import sparse

from apb_msmu.api import MsmuExporter
from conftest import DiannInput

mm = pytest.importorskip("msmu")


def _pipeline(mdata: md.MuData) -> md.MuData:
    """msmu's DIA label-free tutorial, without decoys."""
    mdata = mm.pp.add_filter(mdata, modality="psm", column="q_value", keep="lt", value=0.01)
    mdata = mm.pp.apply_filter(mdata, modality="psm", on="var")
    mdata = mm.pp.log2_transform(mdata, modality="psm")
    mdata = mm.pp.normalize(mdata, modality="psm", method="median")
    mdata = mm.pp.to_peptide(mdata, calculate_q=False)
    mdata = mm.pp.infer_protein(mdata)
    return mm.pp.to_protein(mdata, top_n=3, rank_method="total_intensity", calculate_q=False)


def _dense(matrix: object) -> np.ndarray:
    if isinstance(matrix, sparse.spmatrix | sparse.sparray):
        return np.asarray(sparse.csr_matrix(matrix).toarray(), dtype=float)
    return np.asarray(matrix, dtype=float)


@pytest.mark.parametrize("level", ["peptide", "protein"])
def test_same_features_and_values_as_read_diann(
    diann: DiannInput, tmp_path: Path, level: str
) -> None:
    export = tmp_path / "export.h5mu"
    parsed = ParseRuleCompiler(diann.report, diann.log, requested_levels=("ion",)).compile().parse()
    exported = MsmuExporter(abundance="Precursor_Quantity").export(parsed)
    exported.write_h5mu(export)
    reference = _pipeline(mm.read_diann(str(diann.report)))[level]
    ours = _pipeline(mm.read_h5mu(export))[level]
    if level == "peptide":
        # The export writes every vendor's peptides in ProForma; msmu's reader keeps DIA-NN's
        # notation. Same peptides, so compare under one notation.
        psm = exported["psm"]
        source = psm.varm["search_result"]
        assert isinstance(psm.var, pd.DataFrame)
        assert isinstance(source, pd.DataFrame)
        notation = dict(zip(psm.var["peptide"], source["Modified_Sequence"], strict=True))
        ours.var_names = [notation[name] for name in ours.var_names]

    assert set(ours.obs_names) == set(reference.obs_names)
    assert set(ours.var_names) == set(reference.var_names)
    aligned = ours[reference.obs_names, reference.var_names]
    np.testing.assert_array_equal(_dense(aligned.X), _dense(reference.X))
