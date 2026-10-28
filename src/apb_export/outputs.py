"""Select the target serializer once; APB2 continues to own HDF5 persistence."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from apb2.api import ParsedLevels, write_parsed_levels

from apb_export.long_table import LongTableOutput
from apb_export.perseus_table import PerseusOutput
from apb_export.qpx.output import QpxOutput


class OutputFormat(Protocol):
    """The operations that differ between matrix containers and tables."""

    @property
    def one_level(self) -> bool:
        """Whether the format holds one level: the first listed that the result holds."""
        ...

    def align(self, parsed: ParsedLevels, first_level: str) -> ParsedLevels:
        """Resolve the observation axes the format can represent."""
        ...

    def write(self, exported: ParsedLevels, output: Path) -> None:
        """Serialize the exported quantities."""
        ...


@dataclass(frozen=True)
class MatrixOutput:
    """AnnData/MuData require one shared observation axis, as before."""

    one_level: bool

    def align(self, parsed: ParsedLevels, first_level: str) -> ParsedLevels:
        """Keep the observation group holding the first exported level."""
        return next(group for group in parsed.observation_groups() if first_level in group.levels)

    def write(self, exported: ParsedLevels, output: Path) -> None:
        """Delegate HDF5 and its sidecar to APB2."""
        write_parsed_levels(exported, output, var_name_prefix=False)


_OUTPUTS: dict[str, OutputFormat] = {
    ".h5ad": MatrixOutput(one_level=True),
    ".h5mu": MatrixOutput(one_level=False),
    ".csv": LongTableOutput(),
    ".txt": PerseusOutput(),
    ".qpx": QpxOutput(),
}


def output_for(extension: str) -> OutputFormat:
    """Select the stateless serializer for one declared target output."""
    return _OUTPUTS[extension]
