"""Compile an export rule once, then export APB results with it.

Every entry's sources are resolved and checked against what its group can take before any
row is built, so a result the rule cannot export fails without a partial container. Each
level the rule reads becomes one AnnData; several become one MuData, one modality each.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib.metadata import version
from typing import TYPE_CHECKING, cast

import numpy as np
import polars as pl
from apb2.api import JsonValue, ParsedLevels, UnsJsonCodec
from loguru import logger

from apb_export.computed import compute, dtype
from apb_export.container import Index, Matrix, long_anndata, mudata, wide_anndata
from apb_export.export_rules.loader import EffectiveRule
from apb_export.export_rules.schema import (
    AccessionSyntax,
    AllValues,
    AnnotationMap,
    CatalogueSource,
    ColumnEntry,
    Constant,
    Contaminant,
    EntryBase,
    LayerEntry,
    SourcedEntry,
)
from apb_export.rows import (
    CellRows,
    FeatureRows,
    ObservationRows,
    Rows,
    ScalarRow,
    all_values,
    select_cells,
)
from apb_export.sources import (
    LayerValues,
    LevelSources,
    ObsValues,
    Provenance,
    Resolved,
    ScalarValue,
    VarValues,
    rule_field,
)

if TYPE_CHECKING:
    import anndata as ad
    import mudata as md

NAMESPACE = "apb"
EXPORT_SCHEMA_VERSION = "1"
_OBS_TAKES = (ObsValues, ScalarValue)
_UNS_TAKES = (ScalarValue,)
_VAR_TAKES = {
    "long": (ObsValues, VarValues, LayerValues, ScalarValue),
    "wide": (VarValues, ScalarValue),
}

type Entries = Sequence[ColumnEntry] | Sequence[ColumnEntry | AllValues]


@dataclass(frozen=True, slots=True)
class _Bound:
    """An entry with its sources resolved; ``filled`` writes ``fill_value`` instead."""

    entry: ColumnEntry
    inputs: tuple[Resolved, ...]
    provenance: Provenance
    filled: bool = False


def _describe(provenance: Provenance) -> str:
    return " ".join(provenance[key] for key in ("location", "name") if key in provenance)


def _absent(entry: ColumnEntry, group: str) -> _Bound | None:
    if entry.required:
        raise ValueError(f"{group} entry {entry.name!r}: the result carries no source for it")
    if entry.fill_value is None:
        return None
    return _Bound(entry, (), {"location": "absent", "filled": str(entry.fill_value)}, True)


def _bind(
    entry: ColumnEntry, sources: LevelSources, group: str, takes: tuple[type[object], ...]
) -> _Bound | None:
    """Resolve an entry's sources; ``None`` skips an optional entry with nothing to fill."""
    if isinstance(entry, Constant | AnnotationMap):
        return _Bound(entry, (), {"how": entry.how})
    if isinstance(entry, Contaminant):
        merged = sources.contaminants()
        return _Bound(entry, (merged,), {"how": entry.how, "inputs": _describe(merged.provenance)})
    found = (
        [sources.resolve(entry.source)]
        if isinstance(entry, SourcedEntry)
        else [sources.resolve(reference) for reference in entry.inputs]
    )
    resolved = [item for item in found if item is not None]
    if len(resolved) < len(found):
        return _absent(entry, group)
    for item in resolved:
        if not isinstance(item, takes):
            raise ValueError(
                f"{group} entry {entry.name!r} cannot take {_describe(item.provenance)}"
            )
    if isinstance(entry, SourcedEntry):
        return _Bound(entry, tuple(resolved), resolved[0].provenance)
    inputs = ", ".join(_describe(item.provenance) for item in resolved)
    return _Bound(entry, tuple(resolved), {"how": entry.how, "inputs": inputs})


def _fill(entry: EntryBase) -> object:
    if entry.type == "number" and isinstance(entry.fill_value, str):
        return float(entry.fill_value)
    return entry.fill_value


def _series(bound: _Bound, rows: Rows, syntax: Mapping[str, AccessionSyntax]) -> pl.Series:
    entry = bound.entry
    if bound.filled:
        values = pl.Series([_fill(entry)] * rows.count())
    elif isinstance(entry, Constant):
        values = pl.Series([entry.value] * rows.count())
    elif isinstance(entry, SourcedEntry):
        values = rows.align(bound.inputs[0])
    elif isinstance(entry, AnnotationMap):
        raise ValueError(f"entry {entry.name!r}: annotation_map writes only uns")
    else:
        values = compute(entry, [rows.align(item) for item in bound.inputs], syntax)
    return values.cast(dtype(entry.type, entry.width)).alias(entry.name)


def _nest(values: Sequence[tuple[str, object]]) -> dict[str, object]:
    """Dotted entry names become nested mappings: ``a.b`` is ``{"a": {"b": ...}}``."""
    root: dict[str, object] = {}
    for name, value in values:
        *parents, leaf = name.split(".")
        node = root
        for part in parents:
            child = node.setdefault(part, {})
            if not isinstance(child, dict):
                raise ValueError(f"uns entry {name!r} nests under a value")
            node = child
        node[leaf] = value
    return root


def _repeated(values: VarValues, runs: int) -> LayerValues:
    """A catalogued value per variable as a layer: the same value in every run."""
    column = values.values.cast(pl.Float64)
    frame = pl.DataFrame([column.alias(f"obs_{run}") for run in range(runs)])
    return LayerValues(frame, values.provenance)


def _software(parsed: ParsedLevels) -> set[str]:
    """The software apb2's rules name for the result's levels."""
    names = (rule_field(level, "software_name") for level in parsed.levels.values())
    return {name for name in names if isinstance(name, str)}


def _matrix(values: pl.DataFrame, entry: LayerEntry) -> Matrix:
    """A layer as observations by features, declared missing values blanked."""
    matrix = values.select(pl.all().cast(pl.Float64)).to_numpy().T.copy()
    if entry.missing_values:
        matrix[np.isin(matrix, np.asarray(entry.missing_values, dtype=np.float64))] = np.nan
    return matrix.astype(np.float32 if entry.width == 32 else np.float64)


class _LevelExport:
    """One level of one result, exported by one effective rule."""

    def __init__(self, rule: EffectiveRule, parsed: ParsedLevels, abundance: str | None) -> None:
        self._rule = rule
        self._level = parsed.levels[rule.level]
        self._sources = LevelSources(parsed, rule.level, abundance)
        self._provenance: dict[str, Provenance] = {}

    def long(self) -> ad.AnnData:
        """One row per selected cell, each a feature; X block-diagonal."""
        rule = self._rule.rule
        (entry,) = rule.measurements.layers
        x = self._layer(entry, primary=True)
        if x is None:
            raise ValueError(f"layer {entry.name!r}: the result carries no layer for it")
        bound = self._bind_columns("long")
        uns = self._bind_group("uns", rule.uns, _UNS_TAKES)
        cells = select_cells(self._level, x.values, entry.missing_values)
        obs = self._columns(bound["obs"], ObservationRows(self._level.obs.frame.height))
        var = self._columns(bound["var"], CellRows(cells))
        return self._with_level_part(
            long_anndata(
                obs_index=self._index(obs, rule.axis.obs_keys[0], rule.columns.obs, "obs"),
                obs=self._visible(obs, rule.columns.obs),
                var_index=self._index(var, rule.axis.var_keys[0], rule.columns.var, "var"),
                var=self._visible(var, rule.columns.var),
                x=pl.Series(cells.values).cast(dtype("number", entry.width)).to_numpy(),
                x_obs=cells.obs,
                varm={table.name: all_values(self._level, cells) for table in rule.columns.varm},
                uns=self._uns(uns),
            )
        )

    def wide(self) -> ad.AnnData:
        """APB2's samples-by-features matrices; the primary layer becomes X."""
        rule = self._rule.rule
        primary = rule.measurements.primary_layer
        found = [
            (entry, self._layer(entry, primary=entry.name == primary))
            for entry in rule.measurements.layers
        ]
        bound = self._bind_columns("wide")
        uns = self._bind_group("uns", rule.uns, _UNS_TAKES)
        obs = self._columns(bound["obs"], ObservationRows(self._level.obs.frame.height))
        var = self._columns(bound["var"], FeatureRows(self._level.var.frame.height))
        kept = self._kept(var, rule.columns.var)
        var = [column.filter(kept) for column in var]
        matrices = {
            entry.name: _matrix(values.values.filter(kept), entry)
            for entry, values in found
            if values is not None
        }
        x = matrices.pop(primary)
        return self._with_level_part(
            wide_anndata(
                obs_index=self._index(obs, rule.axis.obs_keys[0], rule.columns.obs, "obs"),
                obs=self._visible(obs, rule.columns.obs),
                var_index=self._index(var, rule.axis.var_keys[0], rule.columns.var, "var"),
                var=self._visible(var, rule.columns.var),
                x=x,
                layers=matrices,
                uns=self._uns(uns),
            )
        )

    def _layer(self, entry: LayerEntry, *, primary: bool) -> LayerValues | None:
        """Resolve a layer; ``None`` for an absent optional one."""
        resolved = self._sources.resolve(entry.source)
        key = f"layers.{entry.name}"
        if resolved is None:
            if primary or entry.required:
                raise ValueError(f"layer {entry.name!r}: the result carries no layer for it")
            self._provenance[key] = {"location": "absent"}
            return None
        if isinstance(resolved, VarValues) and "catalogue" in resolved.provenance:
            resolved = _repeated(resolved, self._level.obs.frame.height)
        if not isinstance(resolved, LayerValues):
            raise ValueError(f"layer {entry.name!r} cannot take {_describe(resolved.provenance)}")
        self._provenance[key] = resolved.provenance
        return resolved

    def _bind_columns(self, shape: str) -> dict[str, list[_Bound | AllValues]]:
        rule = self._rule.rule
        return {
            "obs": self._bind_group("obs", rule.columns.obs, _OBS_TAKES),
            "var": self._bind_group("var", rule.columns.var, _VAR_TAKES[shape]),
        }

    def _bind_group(
        self, group: str, entries: Entries, takes: tuple[type[object], ...]
    ) -> list[_Bound | AllValues]:
        bound: list[_Bound | AllValues] = []
        for entry in entries:
            key = f"{group}.{entry.name}"
            if isinstance(entry, AllValues):
                self._provenance[key] = {"how": entry.how}
                bound.append(entry)
                continue
            item = _bind(entry, self._sources, group, takes)
            self._provenance[key] = {"location": "absent"} if item is None else item.provenance
            if item is not None:
                bound.append(item)
        return bound

    def _columns(self, bound: Sequence[_Bound | AllValues], rows: Rows) -> list[pl.Series]:
        """Evaluate bound entries; an obs ``all_values`` entry adds every annotation column."""
        columns: list[pl.Series] = []
        for item in bound:
            if isinstance(item, AllValues):
                names = self._sources.annotations()
                columns.extend(self._sources.annotation(name) for name in names)
            else:
                columns.append(_series(item, rows, self._rule.rule.accession_syntax))
        names = [column.name for column in columns]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise ValueError(f"level {self._rule.level!r} writes columns {repeated} twice")
        return columns

    def _kept(self, columns: list[pl.Series], entries: Entries) -> pl.Series:
        """Features holding a value in every ``drop_missing`` column; the others are counted."""
        dropping = {e.name for e in entries if not isinstance(e, AllValues) and e.drop_missing}
        kept = pl.Series([True] * self._level.var.frame.height)
        for column in (column for column in columns if column.name in dropping):
            kept &= (column.is_not_null() & (column.cast(pl.String) != "")).fill_null(value=False)
            key = f"var.{column.name}"
            dropped = kept.len() - kept.sum()
            self._provenance[key] = {**self._provenance[key], "dropped": str(dropped)}
            logger.info(
                f"level {self._rule.level}: {dropped} features without {column.name} left out"
            )
        return kept

    @staticmethod
    def _index(columns: list[pl.Series], key: str, entries: Entries, group: str) -> Index:
        (values,) = [column for column in columns if column.name == key]
        if values.is_duplicated().any():
            raise ValueError(f"{group} index {key!r} repeats values")
        (entry,) = [e for e in entries if not isinstance(e, AllValues) and e.name == key]
        named = not isinstance(entry, AllValues) and entry.named_index
        return Index([str(value) for value in values.to_list()], key if named else None)

    @staticmethod
    def _visible(columns: list[pl.Series], entries: Entries) -> pl.DataFrame:
        hidden = {e.name for e in entries if not isinstance(e, AllValues) and e.index_only}
        return pl.DataFrame([column for column in columns if column.name not in hidden])

    def _uns(self, bound: Sequence[_Bound | AllValues]) -> dict[str, object]:
        values: list[tuple[str, object]] = []
        for item in bound:
            if isinstance(item, _Bound):
                values.append((item.entry.name, self._uns_value(item)))
        return _nest(values)

    def _uns_value(self, bound: _Bound) -> object:
        entry = bound.entry
        if isinstance(entry, AnnotationMap):
            return {name: name for name in self._sources.annotations()}
        if isinstance(entry, Constant) and isinstance(entry.value, dict | list):
            # A copy: nested entries are written into it, and the rule is reused.
            return copy.deepcopy(entry.value)
        return _series(bound, ScalarRow(), self._rule.rule.accession_syntax).item()

    def _with_level_part(self, adata: ad.AnnData) -> ad.AnnData:
        """Store the level's APB part, built from what the container actually holds."""
        adata.uns[NAMESPACE] = UnsJsonCodec().encode(self._level_part(adata), {})
        return adata

    def _level_part(self, adata: ad.AnnData) -> dict[str, JsonValue]:
        """The exported level's apb2 records, unchanged, and the export's own record."""
        matrices = (adata.X is not None) + len(adata.layers)
        absent = sum(source.get("location") == "absent" for source in self._provenance.values())
        return {
            "parse": dict(self._level.uns),
            **self._level.metadata,
            "export": {
                "provenance": {
                    "source_level": self._rule.level,
                    "sources": cast(JsonValue, self._provenance),
                },
                "summary": [
                    {
                        "name": "exported_matrices",
                        "label": "Exported matrices",
                        "value": matrices,
                        "unit": "matrices",
                        "status": "ok",
                    },
                    {
                        "name": "absent_entries",
                        "label": "Rule entries the result does not hold",
                        "value": absent,
                        "unit": "entries",
                        "status": "ok",
                    },
                ],
            },
        }


class CompiledExport:
    """One target's effective rules, checked once and applied to any number of results."""

    def __init__(self, rules: Sequence[EffectiveRule]) -> None:
        """Check what Pydantic cannot see across one table's levels.

        Raises:
            ValueError: No level, a long level with more than its primary layer, or a wide
                level with ``varm`` tables.
        """
        if not rules:
            raise ValueError("an export rule needs at least one level")
        output = rules[0].output
        for rule in rules:
            if output.shape == "long" and len(rule.rule.measurements.layers) != 1:
                raise ValueError(f"{rule.document}: long output writes only its primary layer")
            if output.shape == "wide" and rule.rule.columns.varm:
                raise ValueError(f"{rule.document}: wide output writes no varm tables")
        self.rules = tuple(rules)
        self.output = output

    @property
    def levels(self) -> tuple[str, ...]:
        """The APB2 levels the rule reads: those it writes, then those its layers reach into."""
        reached = [
            reference.level
            for rule in self.rules
            for entry in rule.rule.measurements.layers
            for reference in (entry.source if isinstance(entry.source, list) else [entry.source])
            if isinstance(reference, CatalogueSource) and reference.level is not None
        ]
        return tuple(dict.fromkeys([*self.written, *reached]))

    @property
    def written(self) -> tuple[str, ...]:
        """The APB2 levels the rule writes, in rule order; an ``.h5ad`` export writes one."""
        return tuple(rule.level for rule in self.rules)

    @property
    def writes_one_level(self) -> bool:
        """Whether every export writes one level: an ``.h5ad`` rule, or a rule of one level."""
        return self.extension == ".h5ad" or len(self.rules) == 1

    @property
    def extension(self) -> str:
        """The file the target reads: ``.h5ad`` for one AnnData, ``.h5mu`` for a MuData."""
        return self.output.extensions[0]

    def export(
        self, parsed: ParsedLevels, abundance: Mapping[str, str] | None = None
    ) -> ad.AnnData | md.MuData:
        """Write the target's container from one APB2 result.

        Args:
            parsed: An APB2 result holding the rule's required levels.
            abundance: Per level, the layer that becomes X; APB2's primary layer otherwise.

        Raises:
            ValueError: The result lacks a required level or source, a source cannot be taken
                where the rule puts it, or ``abundance`` names a level the rule does not read.
            UnresolvedField: apb-catalog finds an answer ambiguous or the rule unreviewed.
        """
        chosen = dict(abundance or {})
        unknown = sorted(set(chosen) - set(self.written))
        if unknown:
            raise ValueError(f"abundance names levels the rule does not write: {unknown}")
        modalities: dict[str, ad.AnnData] = {}
        for rule in self._exported_by(parsed):
            level = _LevelExport(rule, parsed, chosen.get(rule.level))
            built = level.long() if self.output.shape == "long" else level.wide()
            modalities[rule.rule.modality or rule.level] = built
        root = UnsJsonCodec().encode(self._root_part(parsed), {})
        if self.extension == ".h5ad":
            # As in an apb2 h5ad: the level part under uns[<name>], the root part in uns["apb"].
            ((name, adata),) = modalities.items()
            if name in adata.uns:
                raise ValueError(f"the export rule writes uns[{name!r}] itself")
            adata.uns[name] = {NAMESPACE: adata.uns.pop(NAMESPACE)}
            adata.uns[NAMESPACE] = root
            return adata
        container = mudata(modalities)
        container.uns[NAMESPACE] = root
        return container

    def _root_part(self, parsed: ParsedLevels) -> dict[str, JsonValue]:
        """The result's root apb2 records, unchanged, and what wrote the container from which rule."""
        rule = self.rules[0]
        return {
            "parse": dict(parsed.uns),
            **parsed.metadata,
            "export": {
                "schema_version": EXPORT_SCHEMA_VERSION,
                "provenance": {
                    "package": "apb-export",
                    "package_version": version("apb-export"),
                    "apb2_version": version("apb2"),
                    "target_name": rule.target_name,
                    "target_version_pattern": rule.target_version_pattern,
                    "export_rule": rule.document,
                    "export_rule_file_version": rule.file_version,
                },
            },
        }

    def _exported_by(self, parsed: ParsedLevels) -> list[EffectiveRule]:
        """The levels that export this result: those for its software that it holds.

        An ``.h5ad`` export takes the first of them; an ``.h5mu`` export takes all of them and
        needs every required one.
        """
        software = _software(parsed)
        rules = [
            rule
            for rule in self.rules
            if rule.rule.software is None or software & set(rule.rule.software)
        ]
        held = [rule for rule in rules if rule.level in parsed.levels]
        if self.extension == ".h5ad":
            missing = [] if held else [rule.level for rule in rules]
        else:
            missing = [rule.level for rule in rules if rule not in held and rule.rule.required]
        target = self.rules[0].target_name
        if missing:
            raise ValueError(
                f"{target} export needs APB2's {' or '.join(missing)} level; "
                f"this result has {sorted(parsed.levels)}"
            )
        if not held:
            raise ValueError(f"{target} export reads no level of {sorted(software)} results")
        return held[:1] if self.extension == ".h5ad" else held
