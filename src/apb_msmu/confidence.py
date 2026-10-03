"""The identification confidence msmu filters and scores on, found by meaning.

msmu keeps two confidence columns on its ``psm`` modality: ``q_value``, which its tutorials
filter on, and ``PEP``, which ``to_peptide`` and ``to_protein`` aggregate and use for their
target-decoy q-values. Both are looked up in apb-catalog rather than by vendor column name,
so a new vendor rule needs a catalogue entry, not a change here.

The q-value follows msmu's own DIA-NN reader: the library q-value when match-between-runs
filled it, otherwise the experiment-wide q-value, and the per-run q-value only when the
vendor reports neither. That keeps an export filtered exactly as ``msmu.read_diann`` filters.
"""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl
from apb2.api import FinalLayerTable, ParsedLevels
from apb_catalog.catalog import Catalog

LEVEL = "ion"
CONCEPT = "confidence"
# Library and experiment-wide q-values are catalogued with ProteoBench's entrapment module,
# which ranks DIA-NN uploads by them; the per-run q-value and PEP with identification
# confidence. The kinds share one vocabulary across both catalogues.
_RUN_CATALOGUE = "identification_confidence"
_EXPERIMENT_CATALOGUE = "proteobench_entrapment"


@dataclass(frozen=True, slots=True)
class Confidence:
    """The vendor layers that fill msmu's ``q_value`` and ``PEP`` columns.

    Attributes:
        q_value: The layer msmu's psm filter reads, or ``None`` when the vendor reports none.
        q_value_kind: The catalogue kind it was found as, such as ``global_q_value``.
        pep: The posterior error probability layer, or ``None`` when the vendor reports none.
    """

    q_value: FinalLayerTable | None
    q_value_kind: str | None
    pep: FinalLayerTable | None


def _has_values(layer: FinalLayerTable) -> bool:
    """Say whether a layer holds any finite, nonzero value.

    msmu reads DIA-NN's library q-value only when its column sums to nonzero, because DIA-NN
    writes the column without match-between-runs and leaves it empty.
    """
    values = layer.quantitative_values()
    total = values.select(pl.sum_horizontal(pl.all().fill_nan(0).fill_null(0)).sum()).item()
    return bool(total)


def resolve_confidence(parsed: ParsedLevels) -> Confidence:
    """Find the q-value and PEP layers of the result's ion level.

    Args:
        parsed: An APB2 result holding an ``ion`` level.

    Returns:
        The layers to export, each ``None`` when the vendor does not report it.
    """
    run = Catalog(parsed, _RUN_CATALOGUE)
    experiment = Catalog(parsed, _EXPERIMENT_CATALOGUE)
    library = experiment.layer(LEVEL, concept=CONCEPT, kind="library_q_value")
    candidates = (
        ("library_q_value", library if library is not None and _has_values(library) else None),
        ("global_q_value", experiment.layer(LEVEL, concept=CONCEPT, kind="global_q_value")),
        ("q_value", run.layer(LEVEL, concept=CONCEPT, kind="q_value")),
    )
    kind, q_value = next(
        ((kind, layer) for kind, layer in candidates if layer is not None), (None, None)
    )
    return Confidence(
        q_value=q_value,
        q_value_kind=kind,
        pep=run.layer(LEVEL, concept=CONCEPT, kind="pep"),
    )
