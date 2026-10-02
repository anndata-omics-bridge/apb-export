"""Vendor output to msmu's MuData in one call.

``convert`` reads the vendor files through APB2, finds the confidence layers through
apb-catalog, and returns the MuData msmu's own readers would build. ``to_msmu`` does the same
from an APB2 result already in memory.
"""

from __future__ import annotations

import json
from importlib.metadata import version
from pathlib import Path

import mudata as md
from apb2.api import ParsedLevels, ParseRuleCompiler

from apb_msmu.confidence import Confidence, resolve_confidence
from apb_msmu.container import build_mudata
from apb_msmu.psm import psm_rows

LEVEL = "ion"
NAMESPACE = "apb"


def _software(rule_json: object) -> str:
    rule = json.loads(rule_json) if isinstance(rule_json, str) else {}
    return (
        str(rule.get("software_name", "unknown")).lower() if isinstance(rule, dict) else "unknown"
    )


def _provenance(
    parsed: ParsedLevels, *, abundance: str, confidence: Confidence, source: str
) -> dict[str, object]:
    level = parsed.levels[LEVEL]
    return {
        "level": "precursor",
        "search_engine": _software(level.uns.get("rule_json")),
        "quantification": _software(level.uns.get("rule_json")),
        "label": "label_free",
        "acquisition": None,
        "identification_file": source,
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


def to_msmu(parsed: ParsedLevels, *, abundance: str | None = None, source: str = "") -> md.MuData:
    """Build msmu's MuData from an APB2 result.

    Args:
        parsed: An APB2 result with an ``ion`` level.
        abundance: The layer to become msmu's ``X``; APB2's primary layer when omitted.
        source: The vendor file named in msmu's ``identification_file``.

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
    chosen = abundance or level.primary_layer_name
    confidence = resolve_confidence(parsed)
    rows = psm_rows(level, abundance=chosen, confidence=confidence)
    return build_mudata(
        rows, uns=_provenance(parsed, abundance=chosen, confidence=confidence, source=source)
    )


def convert(
    data: Path,
    params: Path | None,
    /,
    *,
    software: str | None = None,
    abundance: str | None = None,
    strict: bool = False,
) -> md.MuData:
    """Convert vendor output to msmu's MuData.

    Args:
        data: One vendor table or a vendor-result directory.
        params: The vendor's parameter file; ``None`` only for rules that need none.
        software: Restricts recognition to one vendor, as ``apb2 convert --software`` does.
        abundance: The layer to become msmu's ``X``; APB2's primary layer when omitted.
        strict: Promote APB2 layer-contract warnings to errors.

    Returns:
        A MuData with one ``psm`` modality.
    """
    compiler = ParseRuleCompiler(
        data,
        params,
        requested_levels=(LEVEL,),
        checks="strict" if strict else "standard",
        software=software,
    )
    parsed = compiler.compile().parse()
    return to_msmu(parsed, abundance=abundance, source=str(data))
