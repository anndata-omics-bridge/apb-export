"""APB references resolved against one level: what an export rule's sources point at.

A reference resolves to an observation column, a variable column, a layer or a scalar, or to
nothing when the result does not carry it. Catalogue references ask apb-catalog by meaning, so
no vendor column name appears in a rule.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from typing import Literal

import numpy as np
import polars as pl
from apb2.api import FinalLayerTable, ParsedLevel, ParsedLevels
from apb_catalog.api import DataCatalog, UnresolvedField

from apb_export.export_rules.schema import (
    AxisSource,
    CatalogueSource,
    FastaSource,
    MeasurementSource,
    ObsSource,
    Reference,
    RuleSource,
    Scalar,
    Source,
)

type Provenance = dict[str, str]

_FASTA_TABLE = "fasta_validation"
_FASTA_PROTEIN_TABLE = "fasta"
_FASTA_CONTAMINANT = "fasta_matches_contaminant"

_MEMBERS = ";"
"""How APB2 joins a protein group's members, as the vendors write them."""


@dataclass(frozen=True, slots=True)
class ObsValues:
    """One value per observation, in the level's observation order."""

    values: pl.Series
    provenance: Provenance


@dataclass(frozen=True, slots=True)
class VarValues:
    """One value per variable, in the level's variable order."""

    values: pl.Series
    provenance: Provenance


@dataclass(frozen=True, slots=True)
class LayerValues:
    """One row per variable and one column per observation, by position."""

    values: pl.DataFrame
    provenance: Provenance


@dataclass(frozen=True, slots=True)
class ScalarValue:
    """One value for the whole level."""

    value: Scalar
    provenance: Provenance


type Resolved = ObsValues | VarValues | LayerValues | ScalarValue


def _has_values(layer: FinalLayerTable) -> bool:
    """Say whether a layer holds any finite, nonzero value."""
    values = layer.quantitative_values()
    total = values.select(pl.sum_horizontal(pl.all().fill_nan(0).fill_null(0)).sum()).item()
    return bool(total)


def _column_has_values(column: pl.Series) -> bool:
    """Say whether a var column holds any finite, nonzero value."""
    return bool(column.cast(pl.Float64).fill_nan(0).fill_null(0).sum())


def rule_field(level: ParsedLevel, name: str) -> object:
    """A top-level field of the apb2 rule stored with a level; ``None`` when it has none."""
    provenance = level.uns.get("provenance")
    stored = provenance.get("rule_json") if isinstance(provenance, dict) else None
    rule = json.loads(stored) if isinstance(stored, str) else {}
    return rule.get(name) if isinstance(rule, dict) else None


def _scalar(value: object) -> Scalar:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    return json.dumps(value)


class LevelSources:
    """Resolve references against one level of one result; catalogues are opened once."""

    def __init__(self, parsed: ParsedLevels, level: str, abundance: str | None) -> None:
        """Bind the result, the level read, and the caller's abundance layer, if any."""
        self._parsed = parsed
        self._name = level
        self._level = parsed.levels[level]
        self._abundance = abundance
        self._catalogues: dict[tuple[str, str | None], DataCatalog] = {}

    def annotations(self) -> list[str]:
        """The level's observation columns other than its keys, such as ``apb2 annotate`` adds."""
        keys = self._level.obs.key_columns
        return [name for name in self._level.obs.frame.columns if name not in keys]

    def annotation(self, name: str) -> pl.Series:
        """One observation annotation column."""
        return self._level.obs.frame.get_column(name)

    @property
    def abundance_chosen(self) -> bool:
        """Whether the caller named the abundance layer, which then becomes the target's X."""
        return self._abundance is not None

    def abundance_kind(self, layer: str) -> str | None:
        """What the level's layer measures, as the ``abundance`` catalogue's kind; ``None`` for a
        layer the catalogue does not describe."""
        fields = self._catalog("abundance", None).fields()
        kinds = fields.filter(
            (pl.col("level") == self._name)
            & (pl.col("location") == "layers")
            & (pl.col("name") == layer)
            & pl.col("present")
        ).get_column("qualifier.kind")
        return kinds.item(0) if kinds.len() else None

    def resolve(self, source: Source) -> Resolved | None:
        """The first candidate the result carries, or ``None`` when it carries none.

        A catalogue candidate the catalogue reviewed but cannot answer with one field, because
        several match or the entries leave a requested qualifier unknown, passes to the next
        candidate of a list; the answer used records it under ``passed``.

        Raises:
            ValueError: A reference cannot be answered, such as a composite axis key or an
                abundance layer that is absent or carries no abundance role.
            UnresolvedField: apb-catalog cannot answer the last candidate with one field, or
                finds the rule unreviewed.
        """
        candidates: list[Reference] = source if isinstance(source, list) else [source]
        passed: list[str] = []
        for reference in candidates:
            try:
                resolved = self._reference(reference)
            except UnresolvedField as error:
                reviewed = bool(error.resolution.candidates)
                if (
                    not isinstance(reference, CatalogueSource)
                    or not reviewed
                    or reference is candidates[-1]
                ):
                    raise
                passed.append(f"{reference.catalog} {reference.kind} {error.resolution.status}")
                continue
            if resolved is not None:
                if not passed:
                    return resolved
                return replace(
                    resolved, provenance={**resolved.provenance, "passed": "; ".join(passed)}
                )
        return None

    def contaminants(self) -> VarValues:
        """apb2's vendor contaminant marking, merged with apb-fasta's FASTA match if present.

        The marking is the ``apb`` set's kind ``contaminant``, which every apb2 rule writes.

        Raises:
            ValueError: The level lacks the marking.
        """
        marked = self._catalog("apb", None).var(self._name, kind="contaminant")
        if marked is None:
            raise ValueError(f"export needs apb2's contaminant marking on level {self._name!r}")
        provenance: Provenance = {"location": "var", "name": marked.name, "catalog": "apb"}
        fasta = self._level.varm.get(_FASTA_TABLE)
        if fasta is not None and _FASTA_CONTAMINANT in fasta.columns:
            marked = marked | fasta.get_column(_FASTA_CONTAMINANT)
            provenance["fasta"] = f"varm.{_FASTA_TABLE}.{_FASTA_CONTAMINANT}"
        return VarValues(marked, provenance)

    def _reference(self, reference: Reference) -> Resolved | None:
        match reference:
            case AxisSource():
                return self._axis(reference)
            case MeasurementSource():
                return self._measurement()
            case RuleSource():
                return self._rule_field(reference)
            case CatalogueSource():
                return self._catalogue(reference)
            case FastaSource():
                return self._fasta(reference)
            case ObsSource():
                return self._obs(reference)

    def _obs(self, reference: ObsSource) -> ObsValues | None:
        if reference.obs not in self.annotations():
            return None
        values = self.annotation(reference.obs)
        return ObsValues(values, {"location": "obs", "name": reference.obs})

    def _fasta(self, reference: FastaSource) -> VarValues | None:
        for table in (_FASTA_TABLE, _FASTA_PROTEIN_TABLE):
            fasta = self._level.varm.get(table)
            if fasta is not None and reference.fasta in fasta.columns:
                provenance: Provenance = {"location": "varm", "name": f"{table}.{reference.fasta}"}
                return VarValues(fasta.get_column(reference.fasta), provenance)
        return None

    def _axis(self, reference: AxisSource) -> ObsValues | VarValues:
        if reference.axis == "obs_keys":
            keys = self._level.obs.key_columns
            if len(keys) != 1:
                raise ValueError(
                    f"export needs one observation per run; this result keys samples by "
                    f"{keys}, as multiplexed designs do"
                )
            values = self._level.obs.frame.get_column(keys[0]).cast(pl.Utf8)
            return ObsValues(values, {"location": "obs", "name": keys[0]})
        keys = self._level.var.key_columns
        if len(keys) != 1:
            raise ValueError(f"export needs one variable key; this level keys by {keys}")
        values = self._level.var.frame.get_column(keys[0])
        return VarValues(values, {"location": "var", "name": keys[0]})

    def _measurement(self) -> LayerValues:
        name = self._abundance or self._level.primary_layer_name
        self._level.abundance_layers((name,))
        values = self._level.layers[name].decoded_values()
        return LayerValues(values, {"location": "layers", "name": name})

    def _rule_field(self, reference: RuleSource) -> ScalarValue | None:
        value = rule_field(self._level, reference.rule)
        if value is None:
            return None
        return ScalarValue(_scalar(value), {"location": "rule", "name": reference.rule})

    def _catalog(self, name: str, enforce: Literal["MUST", "SHOULD"] | None) -> DataCatalog:
        catalog = self._catalogues.get((name, enforce))
        if catalog is None:
            catalog = DataCatalog(self._parsed, name, enforce=enforce)
            self._catalogues[(name, enforce)] = catalog
        return catalog

    def _catalogue(self, reference: CatalogueSource) -> Resolved | None:
        if reference.level is not None:
            return self._through_protein_group(reference, reference.level)
        catalog = self._catalog(reference.catalog, reference.enforce)
        qualifiers = reference.qualifiers
        found: Provenance = {"catalog": reference.catalog, "kind": reference.kind, **qualifiers}
        if reference.location == "obs":
            values = catalog.obs(self._name, kind=reference.kind, **qualifiers)
            if values is None:
                return None
            return ObsValues(values, {"location": "obs", "name": values.name, **found})
        if reference.location == "var":
            column = catalog.var(self._name, kind=reference.kind, **qualifiers)
            if column is None or (reference.accept == "nonzero" and not _column_has_values(column)):
                return None
            return VarValues(column, {"location": "var", "name": column.name, **found})
        layer = catalog.layer(self._name, kind=reference.kind, **qualifiers)
        if layer is None or (reference.accept == "nonzero" and not _has_values(layer)):
            return None
        provenance = {"location": "layers", "name": layer.layer_name, **found}
        return LayerValues(layer.decoded_values(), provenance)

    def _through_protein_group(self, reference: CatalogueSource, level: str) -> LayerValues | None:
        """A catalogued field of ``level``, seen from this level's cells through protein groups.

        Each variable takes the group its ``protein_assignment`` names; each observation the
        column of the same run. A var column has one value per group, the same in every run,
        so its runs need no counterpart in ``level``.
        """
        assignment = self._level.var.roles.get("protein_assignment")
        if level not in self._parsed.levels or assignment is None:
            return None
        coarse = self._parsed.levels[level]
        catalog = self._catalog(reference.catalog, reference.enforce)
        if reference.location != "layers":
            column = catalog.var(level, kind=reference.kind)
            if column is None or (reference.accept == "nonzero" and not _column_has_values(column)):
                return None
            table = column.cast(pl.Float64).to_numpy()[:, np.newaxis]
            runs = np.zeros(self._level.obs.frame.height, dtype=np.int64)
            name = column.name
        else:
            layer = catalog.layer(level, kind=reference.kind)
            if layer is None or (reference.accept == "nonzero" and not _has_values(layer)):
                return None
            table = layer.decoded_values().select(pl.all().cast(pl.Float64)).to_numpy()
            runs = _positions(
                self._level.obs.frame.get_column(
                    _key(self._level.obs.key_columns, self._name, "run")
                ),
                coarse.obs.frame.get_column(_key(coarse.obs.key_columns, level, "run")),
            )
            name = layer.layer_name
        rows = _groups(
            self._level.var.frame.get_column(assignment),
            coarse.var.frame.get_column(_key(coarse.var.key_columns, level, "variable")),
        )
        gathered = np.full((rows.size, runs.size), np.nan)
        known_rows, known_runs = rows >= 0, runs >= 0
        gathered[np.ix_(known_rows, known_runs)] = table[np.ix_(rows[known_rows], runs[known_runs])]
        values = pl.DataFrame(gathered, schema=[f"obs_{i}" for i in range(runs.size)], orient="row")
        provenance: Provenance = {
            "location": reference.location,
            "name": name,
            "level": level,
            "through": assignment,
            "catalog": reference.catalog,
            "kind": reference.kind,
        }
        return LayerValues(values, provenance)


def _key(columns: tuple[str, ...] | list[str], level: str, axis: str) -> str:
    """A level's one key column on an axis; a composite key cannot be matched across levels."""
    if len(columns) != 1:
        raise ValueError(f"the {level} level keys each {axis} by {list(columns)}, not one column")
    return columns[0]


def _positions(
    wanted: pl.Series, available: pl.Series
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    """Each wanted value's position in ``available``; -1 where it is absent."""
    index = {value: position for position, value in enumerate(available.cast(pl.Utf8).to_list())}
    return np.array(
        [index.get(value, -1) for value in wanted.cast(pl.Utf8).to_list()], dtype=np.int64
    )


def _groups(
    assignments: pl.Series, groups: pl.Series
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    """Each assignment's protein group: the one it names whole, else the one it leads; or -1.

    A group leads with its first member, as MaxQuant lists the razor protein first. A member
    leading several groups names none of them.
    """
    keys: list[str | None] = groups.cast(pl.Utf8).to_list()
    whole = {key: position for position, key in enumerate(keys) if key is not None}
    leads: dict[str, int] = {}
    shared: set[str] = set()
    for position, key in enumerate(keys):
        if key is None:
            continue
        lead = key.split(_MEMBERS, 1)[0]
        if lead in leads:
            shared.add(lead)
        leads[lead] = position
    for lead in shared:
        del leads[lead]
    wanted: list[str | None] = assignments.cast(pl.Utf8).to_list()
    return np.array(
        [-1 if value is None else whole.get(value, leads.get(value, -1)) for value in wanted],
        dtype=np.int64,
    )
