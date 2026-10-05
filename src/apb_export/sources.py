"""APB references resolved against one level: what an export rule's sources point at.

A reference resolves to an observation column, a variable column, a layer or a scalar, or to
nothing when the result does not carry it. Catalogue references ask apb-catalog by meaning, so
no vendor column name appears in a rule.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import polars as pl
from apb2.api import FinalLayerTable, ParsedLevels
from apb_catalog.api import Catalog

from apb_export.export_rules.schema import (
    AxisSource,
    CatalogueSource,
    MeasurementSource,
    Reference,
    RoleSource,
    RuleSource,
    Scalar,
    Source,
)

type Provenance = dict[str, str]


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
        self._catalogues: dict[str, Catalog] = {}

    def annotations(self) -> list[str]:
        """The level's observation columns other than its keys, such as ``apb2 annotate`` adds."""
        keys = self._level.obs.key_columns
        return [name for name in self._level.obs.frame.columns if name not in keys]

    def annotation(self, name: str) -> pl.Series:
        """One observation annotation column."""
        return self._level.obs.frame.get_column(name)

    def resolve(self, source: Source) -> Resolved | None:
        """The first candidate the result carries, or ``None`` when it carries none.

        Raises:
            ValueError: A reference cannot be answered, such as a composite axis key or an
                abundance layer that is absent or carries no abundance role.
            UnresolvedField: apb-catalog finds the answer ambiguous or the rule unreviewed.
        """
        candidates: list[Reference] = source if isinstance(source, list) else [source]
        for reference in candidates:
            resolved = self._reference(reference)
            if resolved is not None:
                return resolved
        return None

    def _reference(self, reference: Reference) -> Resolved | None:
        match reference:
            case str():
                return self._column(reference)
            case RoleSource():
                return self._role(reference)
            case AxisSource():
                return self._axis(reference)
            case MeasurementSource():
                return self._measurement()
            case RuleSource():
                return self._rule_field(reference)
            case CatalogueSource():
                return self._catalogue(reference)

    def _column(self, name: str) -> VarValues | None:
        frame = self._level.var.frame
        if name not in frame.columns:
            return None
        return VarValues(frame.get_column(name), {"location": "var", "name": name})

    def _role(self, reference: RoleSource) -> VarValues | None:
        column = self._level.var.roles.get(reference.role)
        if column is None:
            return None
        values = self._level.var.frame.get_column(column)
        return VarValues(values, {"location": "var", "name": column, "role": reference.role})

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
        stored = self._level.uns.get("rule_json")
        rule = json.loads(stored) if isinstance(stored, str) else {}
        value = rule.get(reference.rule) if isinstance(rule, dict) else None
        if value is None:
            return None
        return ScalarValue(_scalar(value), {"location": "rule", "name": reference.rule})

    def _catalogue(self, reference: CatalogueSource) -> LayerValues | VarValues | None:
        catalog = self._catalogues.get(reference.catalogue)
        if catalog is None:
            catalog = Catalog(self._parsed, reference.catalogue)
            self._catalogues[reference.catalogue] = catalog
        found: Provenance = {"catalogue": reference.catalogue, "kind": reference.kind}
        if reference.location == "var":
            column = catalog.var(self._name, concept=reference.concept, kind=reference.kind)
            if column is None:
                return None
            return VarValues(column, {"location": "var", "name": column.name, **found})
        layer = catalog.layer(self._name, concept=reference.concept, kind=reference.kind)
        if layer is None or (reference.accept == "nonzero" and not _has_values(layer)):
            return None
        provenance = {"location": "layers", "name": layer.layer_name, **found}
        return LayerValues(layer.decoded_values(), provenance)
