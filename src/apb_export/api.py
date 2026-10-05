"""Exports of APB2 results in the formats downstream tools read.

APB2 owns conversion and persistence: ``ParseRuleCompiler`` reads vendor files into
``ParsedLevels``, ``read_parsed_levels`` loads a stored result. An export rule declares what a
target reads; ``MsmuExporter`` applies msmu's rule and builds the MuData msmu's own readers
would. The ``apb-export`` command composes conversion and export in one step.
"""

from __future__ import annotations

from typing import ClassVar

import mudata as md
from apb2.api import ParsedLevels, QuantificationLevel

from apb_export.engine import CompiledExport
from apb_export.export_rules.loader import packaged_rules

_MSMU = CompiledExport(packaged_rules("msmu"))


class MsmuExporter:
    """Bind the abundance choice for exports to msmu's MuData."""

    __slots__ = ("_abundance",)

    level: ClassVar[QuantificationLevel] = _MSMU.levels[0]
    """The quantification level ``export`` reads; the CLI converts only this level."""

    def __init__(self, abundance: str | None = None) -> None:
        """Choose the layer that becomes msmu's ``X``; APB2's primary layer when ``None``."""
        self._abundance = abundance

    def export(self, parsed: ParsedLevels) -> md.MuData:
        """Build msmu's MuData from an APB2 result's ion level, by msmu's export rule.

        Args:
            parsed: An APB2 result holding an ``ion`` level.

        Returns:
            A MuData with one ``psm`` modality, ready for ``msmu.read_h5mu`` once written.

        Raises:
            ValueError: The result has no ion level, or lacks what msmu's psm modality needs.
            UnresolvedField: apb-catalog finds an answer ambiguous or the vendor rule unreviewed.
        """
        chosen = None if self._abundance is None else {self.level: self._abundance}
        exported = _MSMU.export(parsed, abundance=chosen)
        if not isinstance(exported, md.MuData):
            raise TypeError("msmu's export rule must write a MuData")
        return exported


__all__ = ["MsmuExporter"]
