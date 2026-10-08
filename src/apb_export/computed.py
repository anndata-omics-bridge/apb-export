"""Computed entries: what an export rule's ``how`` names, applied to inputs aligned to the rows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import polars as pl

from apb_export.export_rules.schema import (
    Accessions,
    AccessionSyntax,
    ColumnType,
    Contaminant,
    FloatWidth,
    JoinNonempty,
    Lowercase,
    ProformaCharge,
    SequenceLength,
)

type Computed = (
    JoinNonempty | ProformaCharge | SequenceLength | Accessions | Contaminant | Lowercase
)


def dtype(column_type: ColumnType, width: FloatWidth | None) -> pl.DataType:
    """The Polars type an entry is written as; numbers are float64 unless width says 32."""
    match column_type:
        case "string":
            return pl.Utf8()
        case "integer":
            return pl.Int64()
        case "boolean":
            return pl.Boolean()
        case "number":
            return pl.Float32() if width == 32 else pl.Float64()


def _join_nonempty(inputs: Sequence[pl.Series], separator: str) -> pl.Series:
    frame = pl.DataFrame(
        {f"input{index}": series.cast(pl.Utf8) for index, series in enumerate(inputs)}
    )
    parts = [pl.when(pl.col(name) != "").then(pl.col(name)) for name in frame.columns]
    return frame.select(pl.concat_str(parts, separator=separator, ignore_nulls=True)).to_series()


def _is_contaminant(entry: pl.Expr, markers: Sequence[str]) -> pl.Expr:
    if not markers:
        return pl.lit(value=False)
    return pl.any_horizontal(*(entry.str.contains(marker, literal=True) for marker in markers))


def _accession(entry: pl.Expr, syntax: AccessionSyntax) -> pl.Expr:
    """One member as the syntax writes it: ``ACC`` from ``db|ACC|NAME``, contaminants prefixed."""
    fields = entry.str.split("|")
    middle = fields.list.get(1, null_on_oob=True)
    accession = pl.when(fields.list.len() == 3).then(middle).otherwise(entry)
    for marker in syntax.contaminant_markers:
        accession = accession.str.replace(marker, "", literal=True)
    contaminant = _is_contaminant(entry, syntax.contaminant_markers)
    return pl.when(contaminant).then(syntax.contaminant_prefix + accession).otherwise(accession)


def _members(series: pl.Series, syntax: AccessionSyntax) -> tuple[pl.DataFrame, pl.Expr]:
    frame = pl.DataFrame({"members": series.cast(pl.Utf8)})
    return frame, pl.col("members").fill_null("").str.split(syntax.separator)


def _accessions(series: pl.Series, syntax: AccessionSyntax) -> pl.Series:
    frame, members = _members(series, syntax)
    rewritten = members.list.eval(_accession(pl.element(), syntax)).list.join(syntax.separator)
    return frame.select(rewritten).to_series()


def compute(
    entry: Computed, inputs: Sequence[pl.Series], syntax: Mapping[str, AccessionSyntax]
) -> pl.Series:
    """Compute the entry from its aligned inputs."""
    match entry:
        case JoinNonempty():
            return _join_nonempty(inputs, entry.separator)
        case ProformaCharge():
            return inputs[0].cast(pl.Utf8).str.split("/").list.last().cast(pl.Int64)
        case SequenceLength():
            return inputs[0].cast(pl.Utf8).str.len_chars().cast(pl.Int64)
        case Accessions():
            return _accessions(inputs[0], syntax[entry.syntax])
        case Contaminant():
            return inputs[0]
        case Lowercase():
            return inputs[0].cast(pl.Utf8).str.to_lowercase()
