"""msmu's MuData from an APB2 result.

APB2 owns conversion and persistence: ``ParseRuleCompiler`` reads vendor files into
``ParsedLevels``, ``read_parsed_levels`` loads a stored result. ``MsmuExporter`` takes that
result and builds the MuData msmu's own readers would. The ``apb-msmu`` command composes the
two in one step.
"""

from __future__ import annotations

import json
from importlib.metadata import version

import mudata as md
from apb2.api import ParsedLevels

from apb_msmu.confidence import LEVEL, Confidence, resolve_confidence
from apb_msmu.container import build_mudata
from apb_msmu.psm import psm_rows

NAMESPACE = "apb"


def _software(rule_json: object) -> str:
    rule = json.loads(rule_json) if isinstance(rule_json, str) else {}
    return (
        str(rule.get("software_name", "unknown")).lower() if isinstance(rule, dict) else "unknown"
    )


def _settings(parsed: ParsedLevels, *, abundance: str, confidence: Confidence) -> dict[str, object]:
    """msmu's reader settings, and what the export chose, for the psm modality's ``uns``.

    ``identification_file`` stays empty: an APB2 result does not record the files it was read
    from, so whoever converted them fills it in.
    """
    level = parsed.levels[LEVEL]
    software = _software(level.uns.get("rule_json"))
    return {
        "level": "precursor",
        "search_engine": software,
        "quantification": software,
        "label": "label_free",
        "acquisition": None,
        "identification_file": None,
        "quantification_file": None,
        NAMESPACE: {
            "package": "apb-msmu",
            "package_version": version("apb-msmu"),
            "apb2_version": version("apb2"),
            "source_level": LEVEL,
            "abundance_layer": abundance,
            "q_value_layer": None if confidence.q_value is None else confidence.q_value.layer_name,
            "q_value_kind": confidence.q_value_kind,
            "pep_layer": None if confidence.pep is None else confidence.pep.layer_name,
            "rule_json": level.uns.get("rule_json"),
            "plan_json": level.uns.get("plan_json"),
        },
    }


class MsmuExporter:
    """Bind the abundance choice for exports to msmu's MuData."""

    __slots__ = ("_abundance",)

    def __init__(self, *, abundance: str | None = None) -> None:
        """Choose the layer that becomes msmu's ``X``; APB2's primary layer when ``None``."""
        self._abundance = abundance

    def export(self, parsed: ParsedLevels, /) -> md.MuData:
        """Build msmu's MuData from an APB2 result's ion level.

        Args:
            parsed: An APB2 result holding an ``ion`` level.

        Returns:
            A MuData with one ``psm`` modality, ready for ``msmu.read_h5mu`` once written.

        Raises:
            ValueError: The result has no ion level, or lacks what msmu's psm modality needs.
        """
        if LEVEL not in parsed.levels:
            raise ValueError(
                f"msmu export needs APB2's ion level; this result has {sorted(parsed.levels)}"
            )
        level = parsed.levels[LEVEL]
        abundance = self._abundance or level.primary_layer_name
        confidence = resolve_confidence(parsed)
        rows = psm_rows(level, abundance=abundance, confidence=confidence)
        return build_mudata(rows, uns=_settings(parsed, abundance=abundance, confidence=confidence))


__all__ = ["MsmuExporter"]
