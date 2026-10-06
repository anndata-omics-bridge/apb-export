"""Exports of APB2 results in the formats downstream tools read.

APB2 owns conversion and persistence: ``ParseRuleCompiler`` reads vendor files into
``ParsedLevels``, ``read_parsed_levels`` loads a stored result. Each target's packaged export
rule declares what that tool reads; ``Exporter`` applies one and returns the AnnData or MuData
the tool opens. The ``apb-export`` command composes conversion, annotation and export.
"""

from __future__ import annotations

import anndata as ad
import mudata as md
from apb2.api import ParsedLevels, QuantificationLevel

from apb_export.engine import CompiledExport
from apb_export.export_rules.loader import packaged_documents, packaged_rules


class Exporter:
    """Bind one packaged target's export rule and the layer that becomes X."""

    __slots__ = ("_abundance", "_compiled")

    def __init__(self, target: str, abundance: str | None = None) -> None:
        """Choose the target and, for a one-level target, the layer that becomes X.

        Args:
            target: One of :meth:`targets`.
            abundance: A layer of the target's level; APB2's primary layer when ``None``.

        Raises:
            ValueError: The target is not packaged, or ``abundance`` is given for a target that
                writes several levels.
        """
        if target not in self.targets():
            raise ValueError(f"unknown export target {target!r}; packaged: {list(self.targets())}")
        self._compiled = CompiledExport(packaged_rules(target))
        if abundance is not None and not self._compiled.writes_one_level:
            raise ValueError(
                f"{target} writes {list(self._compiled.written)}; "
                "abundance applies only to targets that write one level"
            )
        self._abundance = abundance

    @staticmethod
    def targets() -> tuple[str, ...]:
        """The packaged targets, by name."""
        return tuple(dict.fromkeys(name.split("/", 1)[0] for name in packaged_documents()))

    @property
    def levels(self) -> tuple[QuantificationLevel, ...]:
        """The APB2 levels the target reads; the command asks APB2 for exactly these."""
        return self._compiled.levels

    @property
    def extension(self) -> str:
        """The file the target opens: ``.h5ad`` or ``.h5mu``."""
        return self._compiled.extension

    def export(self, parsed: ParsedLevels) -> ad.AnnData | md.MuData:
        """Build the target's AnnData or MuData from an APB2 result.

        Args:
            parsed: An APB2 result holding the levels the target requires.

        Returns:
            One AnnData for an ``.h5ad`` target; a MuData with one modality per level otherwise.

        Raises:
            ValueError: The result lacks a level or field the target requires.
            UnresolvedField: apb-catalog finds an answer ambiguous or the vendor rule unreviewed.
        """
        written = self._compiled.written
        chosen = None if self._abundance is None else dict.fromkeys(written, self._abundance)
        return self._compiled.export(parsed, abundance=chosen)


__all__ = ["Exporter"]
