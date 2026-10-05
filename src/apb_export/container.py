"""AnnData and MuData assembly: the only module that builds them.

Wide output keeps APB2's samples-by-features matrices. Long output makes each row a feature
with a block-diagonal sparse ``X``, the layout msmu's own readers write. Row-aligned tables go
to ``varm``, made writable to HDF5. Several levels become one MuData, one modality each.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import anndata as ad
import mudata as md
import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
from scipy import sparse

type Matrix = npt.NDArray[np.floating]


@dataclass(frozen=True, slots=True)
class Index:
    """An axis index: one name per row, and the index's own name, if any."""

    values: Sequence[str]
    name: str | None = None

    def pandas(self) -> pd.Index:
        """The index as pandas builds it."""
        return pd.Index(list(self.values), name=self.name)


def _frame(columns: pl.DataFrame, index: pd.Index) -> pd.DataFrame:
    frame = columns.to_pandas() if columns.width else pd.DataFrame(index=index)
    frame.index = index
    return frame


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


def long_anndata(
    *,
    obs_index: Index,
    obs: pl.DataFrame,
    var_index: Index,
    var: pl.DataFrame,
    x: Matrix,
    x_obs: npt.NDArray[np.int64],
    varm: Mapping[str, pl.DataFrame],
    uns: Mapping[str, object],
) -> ad.AnnData:
    """Assemble long rows: each row a feature, X block-diagonal.

    Args:
        obs_index: One name per observation.
        obs: Observation columns, possibly none.
        var_index: One name per row.
        var: The rows' columns.
        x: Each row's X value, already in its written dtype.
        x_obs: Each row's observation position.
        varm: Row-aligned tables.
        uns: The ``uns`` entries.

    Returns:
        One AnnData.
    """
    observations, rows = obs_index.pandas(), var_index.pandas()
    matrix = sparse.csr_matrix(
        (x, (x_obs, np.arange(len(rows)))), shape=(len(observations), len(rows)), dtype=x.dtype
    )
    adata = ad.AnnData(X=matrix, obs=_frame(obs, observations), var=_frame(var, rows))
    for name, table in varm.items():
        adata.varm[name] = _search_result(table, rows)
    adata.uns.update(dict(uns))
    return adata


def wide_anndata(
    *,
    obs_index: Index,
    obs: pl.DataFrame,
    var_index: Index,
    var: pl.DataFrame,
    x: Matrix,
    layers: Mapping[str, Matrix],
    uns: Mapping[str, object],
) -> ad.AnnData:
    """Assemble samples-by-features matrices as APB2 keeps them.

    Args:
        obs_index: One name per observation.
        obs: Observation columns, possibly none.
        var_index: One name per feature.
        var: The features' columns.
        x: The primary matrix, observations by features.
        layers: Further matrices of the same shape.
        uns: The ``uns`` entries.

    Returns:
        One AnnData.
    """
    observations, features = obs_index.pandas(), var_index.pandas()
    adata = ad.AnnData(
        X=x, obs=_frame(obs, observations), var=_frame(var, features), layers=dict(layers)
    )
    adata.uns.update(dict(uns))
    return adata


def mudata(modalities: Mapping[str, ad.AnnData]) -> md.MuData:
    """One MuData holding each level's AnnData as its own modality."""
    return md.MuData(dict(modalities))
