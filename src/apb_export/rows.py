"""Rows an export writes, and resolved sources aligned to them.

Long output writes one row per selected cell: a cell is kept when the X layer's value is
present, finite and not declared missing. Rows run by observation, then by the level's
variable key. Every source is aligned to those rows: variable columns repeat per cell,
observation columns per observation, and layers give the cell's own value.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import polars as pl
from apb2.api import ParsedLevel

from apb_export.sources import LayerValues, ObsValues, Resolved, ScalarValue, VarValues

type IndexArray = npt.NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class Cells:
    """The selected cells, in row order.

    Attributes:
        var_rows: Each cell's variable position in the level.
        obs: Each cell's observation position in the level.
        values: Each cell's X value.
    """

    var_rows: IndexArray
    obs: IndexArray
    values: npt.NDArray[np.float64]


def _key_rank(level: ParsedLevel) -> IndexArray:
    """Each variable's position when the level is sorted by its key columns."""
    keys = list(level.var.key_columns)
    order = level.var.frame.select(keys).with_row_index("row").sort(keys).get_column("row")
    rank = np.empty(level.var.frame.height, dtype=np.int64)
    rank[order.to_numpy()] = np.arange(level.var.frame.height, dtype=np.int64)
    return rank


def select_cells(level: ParsedLevel, x: pl.DataFrame, missing_values: Sequence[float]) -> Cells:
    """Keep the cells whose X value is present, finite and not declared missing."""
    matrix = x.select(pl.all().cast(pl.Float64)).to_numpy()
    kept = np.isfinite(matrix)
    if missing_values:
        kept &= ~np.isin(matrix, np.asarray(missing_values, dtype=np.float64))
    var_rows, obs = np.nonzero(kept)
    order = np.lexsort((_key_rank(level)[var_rows], obs))
    var_rows = var_rows[order].astype(np.int64)
    obs = obs[order].astype(np.int64)
    return Cells(var_rows=var_rows, obs=obs, values=matrix[var_rows, obs])


def layer_at(values: pl.DataFrame, cells: Cells, name: str) -> pl.Series:
    """A layer's value at each cell, its type kept."""
    positional = values.rename({column: str(index) for index, column in enumerate(values.columns)})
    long = (
        positional.with_row_index("var_row")
        .unpivot(index="var_row", variable_name="obs", value_name=name)
        .with_columns(pl.col("var_row").cast(pl.Int64), pl.col("obs").cast(pl.Int64))
    )
    wanted = pl.DataFrame({"var_row": cells.var_rows, "obs": cells.obs})
    joined = wanted.join(long, on=["var_row", "obs"], how="left", maintain_order="left")
    return joined.get_column(name)


def all_values(level: ParsedLevel, cells: Cells) -> pl.DataFrame:
    """Every variable column, the key first, then every layer's value at each cell."""
    keys = list(level.var.key_columns)
    ordered = [*keys, *(name for name in level.var.frame.columns if name not in keys)]
    table = level.var.frame.select(pl.col(ordered).gather(cells.var_rows))
    for name, layer in level.layers.items():
        table = table.with_columns(layer_at(layer.decoded_values(), cells, name))
    return table


class CellRows:
    """One row per selected cell."""

    def __init__(self, cells: Cells) -> None:
        """Align to these cells."""
        self._cells = cells

    def count(self) -> int:
        """The number of rows."""
        return len(self._cells.var_rows)

    def align(self, resolved: Resolved) -> pl.Series:
        """The resolved values at each cell."""
        match resolved:
            case VarValues():
                return resolved.values.gather(self._cells.var_rows)
            case ObsValues():
                return resolved.values.gather(self._cells.obs)
            case LayerValues():
                return layer_at(resolved.values, self._cells, "value")
            case ScalarValue():
                return pl.Series([resolved.value] * self.count())


class FeatureRows:
    """One row per variable of the level, as wide output writes them."""

    def __init__(self, count: int) -> None:
        """Align to this many variables."""
        self._count = count

    def count(self) -> int:
        """The number of rows."""
        return self._count

    def align(self, resolved: Resolved) -> pl.Series:
        """Variable values as they are; a scalar repeated."""
        match resolved:
            case VarValues():
                return resolved.values
            case ScalarValue():
                return pl.Series([resolved.value] * self._count)
            case _:
                raise ValueError(f"a variable row cannot take {resolved.provenance}")


class ObservationRows:
    """One row per observation of the level."""

    def __init__(self, count: int) -> None:
        """Align to this many observations."""
        self._count = count

    def count(self) -> int:
        """The number of rows."""
        return self._count

    def align(self, resolved: Resolved) -> pl.Series:
        """Observation values as they are; a scalar repeated."""
        match resolved:
            case ObsValues():
                return resolved.values
            case ScalarValue():
                return pl.Series([resolved.value] * self._count)
            case _:
                raise ValueError(f"an observation row cannot take {resolved.provenance}")


class ScalarRow:
    """The single row of ``uns`` entries."""

    def count(self) -> int:
        """Always one."""
        return 1

    def align(self, resolved: Resolved) -> pl.Series:
        """A scalar as a one-value series."""
        if not isinstance(resolved, ScalarValue):
            raise ValueError(f"a uns entry cannot take {resolved.provenance}")
        return pl.Series([resolved.value])


type Rows = CellRows | FeatureRows | ObservationRows | ScalarRow
