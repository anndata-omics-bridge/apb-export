"""msmu's MuData, built from psm rows.

The layout is what ``msmu.read_diann`` writes, so msmu takes the file with
``msmu.read_h5mu`` and runs unchanged: a single ``psm`` modality with runs as observations,
one feature per observed (run, precursor) cell, a block-diagonal sparse float32 ``X``, the
curated columns in ``var`` and the complete source values in ``varm["search_result"]``.
"""

from __future__ import annotations

from collections.abc import Mapping

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse

from apb_msmu.psm import PsmRows

MODALITY = "psm"
_HIDDEN = ("feature_id", "run_index", "abundance")


def _search_result(frame: pl.DataFrame, index: pd.Index) -> pd.DataFrame:
    """The source values as pandas, made writable to HDF5.

    h5py rejects "/" in dataset keys, and AnnData stores object columns only as strings.
    """
    table = frame.to_pandas()
    table.columns = [str(name).replace("/", "_") for name in table.columns]
    for name in table.columns:
        if table[name].dtype == object:
            table[name] = table[name].astype(str)
    table.index = index
    return table


def build_mudata(rows: PsmRows, *, uns: Mapping[str, object]) -> md.MuData:
    """Assemble the msmu container.

    Args:
        rows: The psm rows and their source values.
        uns: Entries for the psm modality's ``uns``, msmu's reader settings among them.

    Returns:
        A MuData with one ``psm`` modality.
    """
    index = pd.Index(rows.var["feature_id"].to_list())
    matrix = sparse.csr_matrix(
        (
            rows.var["abundance"].to_numpy().astype(np.float32),
            (rows.var["run_index"].to_numpy(), np.arange(rows.var.height)),
        ),
        shape=(len(rows.runs), rows.var.height),
        dtype=np.float32,
    )
    var = rows.var.drop(_HIDDEN).to_pandas()
    var.index = index
    psm = ad.AnnData(X=matrix, obs=pd.DataFrame(index=pd.Index(rows.runs)), var=var)
    psm.varm["search_result"] = _search_result(rows.search_result, index)
    psm.uns.update(dict(uns))
    return md.MuData({MODALITY: psm})
