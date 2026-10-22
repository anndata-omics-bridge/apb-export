"""Rows an export writes, and resolved sources aligned to them.

Variable rows take variable values, observation rows observation values, and the single
``uns`` row a scalar; a scalar repeats on any row.
"""

from __future__ import annotations

import polars as pl

from apb_export.sources import ObsValues, Resolved, ScalarValue, VarValues


class FeatureRows:
    """One row per variable of the level."""

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


type Rows = FeatureRows | ObservationRows | ScalarRow
