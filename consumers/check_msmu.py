"""msmu reads apb-export's file and gets what its own DIA-NN reader gets from the report.

Inside the consumer image: ``python check_msmu.py EXAMPLES``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import msmu as mm
import mudata as md
import numpy as np
from scipy import sparse


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


def check(folder: Path, level: str) -> None:
    """Same samples, features and values at one level of msmu's pipeline."""
    export = folder / "msmu.h5mu"
    reference = _pipeline(mm.read_diann(str(folder / "report.tsv")))[level]
    ours = _pipeline(mm.read_h5mu(str(export)))[level]
    if level == "peptide":
        # The export writes peptides in ProForma; msmu's reader keeps DIA-NN's notation.
        psm = md.read_h5mu(export)["psm"]
        source = psm.varm["search_result"]
        notation = dict(zip(psm.var["peptide"], source["Modified_Sequence"], strict=True))
        ours.var_names = [notation[name] for name in ours.var_names]
    assert set(ours.obs_names) == set(reference.obs_names), f"{level}: samples differ"
    assert set(ours.var_names) == set(reference.var_names), f"{level}: features differ"
    aligned = ours[reference.obs_names, reference.var_names]
    np.testing.assert_array_equal(_dense(aligned.X), _dense(reference.X))
    print(f"msmu {level}: {reference.n_vars} features, identical to read_diann")


if __name__ == "__main__":
    for name in ("peptide", "protein"):
        check(Path(sys.argv[1]), name)
