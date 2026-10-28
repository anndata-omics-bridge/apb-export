"""The long table: one row per observation and feature, the layers side by side."""

from __future__ import annotations

from pathlib import Path

import polars as pl
from apb2.api import ParsedLevel, ParsedLevels

_POSITIONS = ["var_pos", "obs_pos"]


def _layer(level: ParsedLevel, name: str) -> pl.DataFrame:
    """One layer's stored cells, its value named as the layer: category codes as the labels
    the parsing rule declared, NaN read as missing."""
    layer = level.layers[name]
    values = layer.semantics.decode(layer.cells.select(pl.col("value").alias(name)))
    cells = layer.cells.select(_POSITIONS).hstack(values)
    if cells.schema[name].is_float():
        cells = cells.with_columns(pl.col(name).fill_nan(None))
    return cells.filter(pl.col(name).is_not_null())


def _cells(level: ParsedLevel) -> pl.DataFrame:
    """Every position some layer holds a value at, with each layer's value or null."""
    frames = [_layer(level, name) for name in level.layers]
    cells = frames[0]
    for frame in frames[1:]:
        cells = cells.join(frame, on=_POSITIONS, how="full", coalesce=True)
    return cells


def long_table(level: ParsedLevel) -> pl.DataFrame:
    """The level's obs columns, var columns and layers, one row per observed cell.

    Raises:
        ValueError: A column name repeats across obs, var and layers.
    """
    obs, var, layers = level.obs.frame, level.var.frame, list(level.layers)
    names = [*obs.columns, *var.columns, *layers]
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise ValueError(f"the long table would name {repeated} twice")
    return (
        _cells(level)
        .sort(_POSITIONS)
        .join(var.with_row_index("var_pos"), on="var_pos", how="left", maintain_order="left")
        .join(obs.with_row_index("obs_pos"), on="obs_pos", how="left", maintain_order="left")
        .select(*obs.columns, *var.columns, *layers)
    )


class LongTableOutput:
    """A CSV of one level, as MSstats reads one row per precursor and run."""

    @property
    def one_level(self) -> bool:
        """The table holds the first listed level the result holds."""
        return True

    def align(self, parsed: ParsedLevels, first_level: str) -> ParsedLevels:
        """Keep the observation group holding the exported level."""
        return next(group for group in parsed.observation_groups() if first_level in group.levels)

    def write(self, exported: ParsedLevels, output: Path) -> None:
        """Write the long table, missing values as ``NA``; no sidecar.

        Raises:
            ValueError: The export holds more than one level.
        """
        if len(exported.levels) != 1:
            raise ValueError(f"a long table holds one level, not {sorted(exported.levels)}")
        (level,) = exported.levels.values()
        output.parent.mkdir(parents=True, exist_ok=True)
        long_table(level).write_csv(output, null_value="NA")
