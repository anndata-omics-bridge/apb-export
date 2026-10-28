"""The Perseus table: one row per feature, ``#!{Type}`` and ``#!{C:<name>}`` header rows."""

from __future__ import annotations

from pathlib import Path

import polars as pl
from apb2.api import ParsedLevel, ParsedLevels

_EMPTY = ""


def _var_column(column: pl.Series) -> tuple[pl.Series, str]:
    """A var column as Perseus reads it, and its type: a flag as ``+`` or empty, categorical;
    a number numerical; anything else text."""
    if column.dtype == pl.Boolean:
        marked = pl.select(pl.when(column).then(pl.lit("+")).otherwise(pl.lit(_EMPTY)))
        return marked.to_series().alias(column.name), "C"
    if column.dtype.is_numeric():
        return column.cast(pl.Float64).fill_nan(None), "N"
    return column.cast(pl.Utf8), "T"


def _sample_columns(level: ParsedLevel, name: str, samples: list[str]) -> pl.DataFrame:
    """One layer, a ``<layer> <sample>`` column per sample, NaN read as missing."""
    values = level.layers[name].quantitative_values().cast(pl.Float64).fill_nan(None)
    return values.rename(
        dict(zip(values.columns, [f"{name} {sample}" for sample in samples], strict=True))
    )


def perseus_table(level: ParsedLevel) -> tuple[pl.DataFrame, list[str], dict[str, list[str]]]:
    """The level's var columns, then each layer as ``<layer> <sample>`` columns.

    Returns:
        The table, one type letter per column, and per obs annotation column its value under
        each sample column, empty under the others.

    Raises:
        ValueError: A column name repeats.
    """
    (key,) = level.obs.key_columns
    obs = level.obs.frame.with_columns(pl.all().cast(pl.Utf8))
    samples = obs.get_column(key).to_list()
    columns, types = zip(
        *(_var_column(column) for column in level.var.frame.iter_columns()), strict=True
    )
    layers = [_sample_columns(level, name, samples) for name in level.layers]
    names = [
        *(column.name for column in columns),
        *(name for layer in layers for name in layer.columns),
    ]
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise ValueError(f"the Perseus table would name {repeated} twice")
    table = pl.DataFrame(list(columns)).hstack([column for layer in layers for column in layer])
    categories = {
        name: [_EMPTY] * len(columns)
        + obs.get_column(name).fill_null(_EMPTY).to_list() * len(layers)
        for name in obs.columns
        if name != key
    }
    return table, [*types, *["E"] * (table.width - len(columns))], categories


class PerseusOutput:
    """A tab-separated table of one level, as Perseus loads MaxQuant's proteinGroups.txt."""

    @property
    def one_level(self) -> bool:
        """The table holds the first listed level the result holds."""
        return True

    def align(self, parsed: ParsedLevels, first_level: str) -> ParsedLevels:
        """Keep the observation group holding the exported level."""
        return next(group for group in parsed.observation_groups() if first_level in group.levels)

    def write(self, exported: ParsedLevels, output: Path) -> None:
        """Write the Perseus table, missing numbers as ``NaN``; no sidecar.

        Raises:
            ValueError: The export holds more than one level.
        """
        if len(exported.levels) != 1:
            raise ValueError(f"a Perseus table holds one level, not {sorted(exported.levels)}")
        (level,) = exported.levels.values()
        table, types, categories = perseus_table(level)
        header = [
            "\t".join(table.columns),
            "#!{Type}" + "\t".join(types),
            *("#!{C:" + name + "}" + "\t".join(values) for name, values in categories.items()),
        ]
        numbers = [
            name for name, kind in zip(table.columns, types, strict=True) if kind in {"E", "N"}
        ]
        body = table.with_columns(pl.col(numbers).fill_null(float("nan"))).write_csv(
            separator="\t", include_header=False, null_value=_EMPTY, quote_style="never"
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(header) + "\n" + body, encoding="utf-8")
