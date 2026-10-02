"""APB2's ion level as the rows of msmu's ``psm`` modality.

msmu's ``psm`` modality holds one feature per observed (run, precursor) cell: its readers
take a long search-engine report as it comes, and ``to_peptide`` summarises across runs
later. APB2's ion level is the wide (runs x precursors) matrix of the same report, so this
module unpivots its observed cells back into rows. Everything here is Polars; containers
are built in :mod:`apb_msmu.container`.
"""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl
from apb2.parserV2.parse_quant.data.parsed import (
    CategoricalLayerSemantics,
    FinalLayerTable,
    ParsedLevel,
)

from apb_msmu.confidence import Confidence

# APB2 rules name the ProForma columns of the ion level by this convention: the precursor
# ("SEQ[UNIMOD:4]K/2"), its modified sequence, and its unmodified sequence.
PRECURSOR = "ProForma_ion"
PEPTIDOFORM = "ProForma_peptidoform"
PEPTIDE = "ProForma_peptide"
ACCESSIONS_ROLE = "fasta_accessions"

# msmu's readers canonicalise protein groups this way: members split on ";", the accession
# taken from a UniProt "db|ACC|NAME" entry, a contaminant written "Cont_ACC" whichever marker
# its FASTA used. Repeated here so an export needs no msmu install and matches its readers.
GROUP_SEPARATOR = ";"
CONTAMINANT_MARKERS = ("contam_", "Cont_", "CON__")
CONTAMINANT_PREFIX = "Cont_"

_OBS = "obs"
_VALUE = "abundance"


@dataclass(frozen=True, slots=True)
class PsmRows:
    """One row per observed (run, precursor) cell, in msmu's vocabulary.

    Attributes:
        var: msmu's psm ``var`` columns, then ``run_index`` and ``abundance``.
        search_result: Every APB2 feature column and layer value of the same cells, in
            ``var`` order, for msmu's ``varm["search_result"]``.
        runs: Run names in APB2 observation order, which ``run_index`` refers to.
    """

    var: pl.DataFrame
    search_result: pl.DataFrame
    runs: tuple[str, ...]


def _runs(level: ParsedLevel) -> tuple[str, ...]:
    keys = level.obs.key_columns
    if len(keys) != 1:
        raise ValueError(
            f"msmu export needs one observation per run; this result keys samples by {keys}, "
            "as multiplexed designs do"
        )
    return tuple(str(run) for run in level.obs.frame[keys[0]].to_list())


def _cells(layer: FinalLayerTable, name: str) -> pl.DataFrame:
    """Unpivot one layer to (precursor, obs, value), the observation as its integer index."""
    if layer.var_key_columns != (PRECURSOR,):
        raise ValueError(f"layer {layer.layer_name!r} is not keyed by {PRECURSOR}")
    cells = layer.values.unpivot(index=PRECURSOR, variable_name=_OBS, value_name=name)
    cells = cells.with_columns(pl.col(_OBS).str.strip_prefix("obs_").cast(pl.Int64))
    semantics = layer.semantics
    if isinstance(semantics, CategoricalLayerSemantics):
        labels = {code: label for label, code in semantics.categories}
        cells = cells.with_columns(
            pl.col(name).replace_strict(labels, default=None, return_dtype=pl.Utf8)
        )
    return cells


def _is_contaminant(entry: pl.Expr) -> pl.Expr:
    return pl.any_horizontal(*(entry.str.contains(m, literal=True) for m in CONTAMINANT_MARKERS))


def _accession(entry: pl.Expr) -> pl.Expr:
    """One group member as msmu's readers write it: ``ACC`` or ``Cont_ACC``."""
    fields = entry.str.split("|")
    middle = fields.list.get(1, null_on_oob=True)
    accession = pl.when(fields.list.len() == 3).then(middle).otherwise(entry)
    for marker in CONTAMINANT_MARKERS:
        accession = accession.str.replace(marker, "", literal=True)
    return pl.when(_is_contaminant(entry)).then(CONTAMINANT_PREFIX + accession).otherwise(accession)


def _protein_columns(column: str) -> list[pl.Expr]:
    members = pl.col(column).cast(pl.Utf8).fill_null("").str.split(GROUP_SEPARATOR)
    return [
        members.list.eval(_accession(pl.element())).list.join(GROUP_SEPARATOR).alias("proteins"),
        members.list.eval(_is_contaminant(pl.element()))
        .list.any()
        .cast(pl.Int64)
        .alias("contaminant"),
    ]


def _confidence_columns(cells: pl.DataFrame, confidence: Confidence) -> pl.DataFrame:
    for name, layer in (("q_value", confidence.q_value), ("PEP", confidence.pep)):
        if layer is None:
            # msmu's to_peptide aggregates PEP unconditionally; NaN says "not reported".
            if name == "PEP":
                cells = cells.with_columns(pl.lit(float("nan"), pl.Float32).alias(name))
            continue
        values = _cells(layer, name).with_columns(pl.col(name).cast(pl.Float32))
        cells = cells.join(values, on=[PRECURSOR, _OBS], how="left")
    return cells


def _search_result(level: ParsedLevel, cells: pl.DataFrame) -> pl.DataFrame:
    """Every feature column and every layer's value at the exported cells."""
    result = cells.select(PRECURSOR, _OBS).join(level.var.frame, on=PRECURSOR, how="left")
    for name, layer in level.layers.items():
        result = result.join(_cells(layer, name), on=[PRECURSOR, _OBS], how="left")
    return result.drop(_OBS)


def psm_rows(level: ParsedLevel, *, abundance: str, confidence: Confidence) -> PsmRows:
    """Build msmu's psm rows from APB2's ion level.

    Args:
        level: The APB2 ion level.
        abundance: The layer that becomes msmu's ``X``.
        confidence: The q-value and PEP layers to carry.

    Returns:
        The rows of msmu's psm modality and the cells' complete APB2 values.
    """
    runs = _runs(level)
    missing = [c for c in (PRECURSOR, PEPTIDOFORM, PEPTIDE) if c not in level.var.frame.columns]
    if missing:
        raise ValueError(f"ion level lacks the ProForma columns {missing}")
    if abundance not in level.layers:
        raise ValueError(f"no abundance layer {abundance!r}; layers are {sorted(level.layers)}")
    accessions = level.uns.get("column_roles", {})
    if not isinstance(accessions, dict) or ACCESSIONS_ROLE not in accessions:
        raise ValueError(f"ion level declares no {ACCESSIONS_ROLE!r} column role")
    accession_column = str(accessions[ACCESSIONS_ROLE])

    cells = _cells(level.layers[abundance], _VALUE).filter(
        pl.col(_VALUE).is_finite() & (pl.col(_VALUE) != 0)
    )
    cells = cells.join(
        level.var.frame.select(PRECURSOR, PEPTIDOFORM, PEPTIDE, accession_column),
        on=PRECURSOR,
        how="left",
    ).sort(_OBS, PRECURSOR, maintain_order=True)
    cells = _confidence_columns(cells, confidence)
    run_frame = pl.DataFrame(
        {_OBS: range(len(runs)), "filename": runs}, schema={_OBS: pl.Int64, "filename": pl.Utf8}
    )
    cells = cells.join(run_frame, on=_OBS, how="left", maintain_order="left")
    proteins, contaminant = _protein_columns(accession_column)
    # Column order follows msmu's DIA-NN reader.
    var = cells.select(
        (pl.col("filename") + "." + pl.col(PRECURSOR)).alias("feature_id"),
        proteins,
        pl.col(PEPTIDOFORM).alias("peptide"),
        pl.col(PEPTIDE).alias("stripped_peptide"),
        pl.col("filename"),
        pl.col(PRECURSOR).str.split("/").list.last().cast(pl.Int64).alias("charge"),
        pl.col(PEPTIDE).str.len_chars().cast(pl.Int64).alias("peptide_length"),
        pl.lit(0, pl.Int64).alias("decoy"),
        contaminant,
        *(pl.col(c) for c in ("PEP", "q_value") if c in cells.columns),
        pl.col(_OBS).alias("run_index"),
        pl.col(_VALUE).cast(pl.Float32).alias(_VALUE),
    )
    return PsmRows(var=var, search_result=_search_result(level, cells), runs=runs)
