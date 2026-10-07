"""Plot-ready QC of one quantitative matrix: counts, intensity densities, CV distributions.

Values are linear abundances; anything not finite or not positive counts as missing. The page
receives binned summaries only, never the matrix.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

type Floats = npt.NDArray[np.float64]

_DENSITY_BINS = 80
_CV_BINS = 50
_CV_MAX_PERCENT = 200.0
ALL_SAMPLES = "all samples"


@dataclass(frozen=True, slots=True)
class Matrix:
    """One layer: one row per feature, one column per sample."""

    level: str
    layer: str
    samples: list[str]
    values: Floats
    # One label per sample from the annotation column chosen for grouping; None without one.
    groups: list[str] | None = None
    grouping: str | None = None


def _round(values: Floats) -> list[float]:
    return [float(value) for value in np.round(values, 6)]


def _present(values: Floats) -> npt.NDArray[np.bool_]:
    return np.isfinite(values) & (values > 0)


def _density(
    log2: Floats, present: npt.NDArray[np.bool_], samples: Sequence[str]
) -> dict[str, object]:
    observed = log2[present]
    if observed.size == 0:
        return {"x": [], "series": []}
    edges = np.linspace(float(observed.min()), float(observed.max()) + 1e-9, _DENSITY_BINS + 1)
    series = [
        {
            "sample": sample,
            "y": _round(np.histogram(log2[present[:, column], column], bins=edges, density=True)[0])
            if present[:, column].any()
            else [0.0] * _DENSITY_BINS,
        }
        for column, sample in enumerate(samples)
    ]
    return {"x": _round((edges[:-1] + edges[1:]) / 2), "series": series}


def _cv_percent(values: Floats) -> Floats:
    """Per-feature CV over the samples given, for features seen in at least two of them."""
    present = _present(values)
    counts = present.sum(axis=1)
    kept = values[counts >= 2]
    masked = np.where(_present(kept), kept, np.nan)
    mean = np.nanmean(masked, axis=1)
    sd = np.nanstd(masked, axis=1, ddof=1)
    return 100.0 * sd / mean


def _groups(matrix: Matrix) -> dict[str, list[int]]:
    if matrix.groups is None:
        return {ALL_SAMPLES: list(range(len(matrix.samples)))}
    members: dict[str, list[int]] = {}
    for column, group in enumerate(matrix.groups):
        members.setdefault(group, []).append(column)
    return members


def _cv(matrix: Matrix) -> dict[str, object]:
    edges = np.linspace(0.0, _CV_MAX_PERCENT, _CV_BINS + 1)
    series: list[dict[str, object]] = []
    for group, columns in _groups(matrix).items():
        cv = _cv_percent(matrix.values[:, columns]) if len(columns) >= 2 else np.empty(0)
        clipped = np.minimum(cv, _CV_MAX_PERCENT - 1e-9)
        series.append(
            {
                "group": group,
                "samples": len(columns),
                "features": int(cv.size),
                "median": None if cv.size == 0 else round(float(np.median(cv)), 3),
                "counts": [int(count) for count in np.histogram(clipped, bins=edges)[0]],
            }
        )
    return {
        "grouping": matrix.grouping or ALL_SAMPLES,
        "x": _round((edges[:-1] + edges[1:]) / 2),
        "series": series,
    }


def summarize(matrix: Matrix) -> dict[str, object]:
    """Counts per sample, log2 intensity densities and CV histograms for one matrix."""
    present = _present(matrix.values)
    with np.errstate(divide="ignore", invalid="ignore"):
        log2 = np.log2(np.where(present, matrix.values, np.nan))
    features = matrix.values.shape[0]
    detected = present.sum(axis=0)
    return {
        "level": matrix.level,
        "layer": matrix.layer,
        "features": features,
        "samples": [
            {
                "name": sample,
                "group": None if matrix.groups is None else matrix.groups[column],
                "detected": int(detected[column]),
                "missing_fraction": round(1.0 - float(detected[column]) / features, 4)
                if features
                else None,
            }
            for column, sample in enumerate(matrix.samples)
        ],
        "density": _density(log2, present, matrix.samples),
        "cv": _cv(matrix),
    }


def _key(name: str) -> str:
    """A column name compared loosely: apb2 writes "factor value[x]" as "factor_value_x"."""
    return re.sub(r"[^0-9a-z]+", "_", name.lower()).strip("_")


def grouping_column(
    columns: Sequence[str], sample_table: dict[str, list[str]]
) -> tuple[str, list[str]] | None:
    """The first annotation column in the sample table that splits samples into groups.

    SDRF ``factor value[...]`` columns come first, since they name the design. A run key
    names every sample once and a constant such as "not available" names none, so neither
    groups. Names match as apb2 writes them into the sample table.
    """
    table = {_key(name): (name, values) for name, values in sample_table.items()}
    ordered = sorted(columns, key=lambda column: not _key(column).startswith("factor_value"))
    for column in ordered:
        found = table.get(_key(column))
        if found is not None and 1 < len(set(found[1])) < len(found[1]):
            return found
    return None
