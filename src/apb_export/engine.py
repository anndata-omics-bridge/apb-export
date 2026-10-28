"""Compile an export rule once, then export APB results with it.

Every entry's sources are resolved and checked against what its group can take before any
row is built, so a result the rule cannot export fails without a partial container. Each
level the rule writes becomes one level of a ``ParsedLevels`` for APB2's writer, and a rule's
feature links become one of its directed feature relations.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from importlib.metadata import version
from pathlib import Path
from typing import cast

import polars as pl
from apb2.api import JsonValue, ParsedLevel, ParsedLevels
from loguru import logger

from apb_export.computed import compute, dtype
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
    FeatureLink,
    LayerEntry,
    LayerMinimum,
    MeasurementSource,
    SourcedEntry,
)
from apb_export.outputs import output_for
from apb_export.rows import FeatureRows, ObservationRows, Rows, ScalarRow
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

EXPORT_SCHEMA_VERSION = "1"
_CALLERS_LAYER = MeasurementSource(measurements="primary_layer")
"""The primary layer's source when the caller names the abundance layer."""
_OBS_TAKES = (ObsValues, ScalarValue)
_UNS_TAKES = (ScalarValue,)
_VAR_TAKES = (VarValues, ScalarValue)
_LINKS = pl.Schema(
    {
        "source": pl.String,
        "row": pl.Int64,
        "target": pl.String,
        "column": pl.Int64,
        "value": pl.Float64,
    }
)

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
    if isinstance(entry, LayerMinimum):
        found_layer = sources.resolve(entry.source)
        if found_layer is None:
            return _absent(entry, group)
        reduced = _minimum(found_layer)
        if not isinstance(reduced, takes):
            raise ValueError(
                f"{group} entry {entry.name!r} cannot take {_describe(reduced.provenance)}"
            )
        return _Bound(entry, (reduced,), {**reduced.provenance, "how": entry.how})
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


def _minimum(resolved: Resolved) -> Resolved:
    """A layer reduced to its smallest value per variable; anything else as it is."""
    if not isinstance(resolved, LayerValues):
        return resolved
    values = resolved.values.select(pl.min_horizontal(pl.all().cast(pl.Float64))).to_series()
    return VarValues(values, resolved.provenance)


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


def _layer_frame(values: pl.DataFrame, entry: LayerEntry) -> pl.DataFrame:
    """A layer as APB2 keeps it, features by observations, declared missing values blanked."""
    numbers = values.select(pl.all().cast(pl.Float64))
    if entry.missing_values:
        absent = pl.all().is_in(entry.missing_values)
        numbers = numbers.select(pl.when(absent).then(None).otherwise(pl.all()).name.keep())
    return numbers.cast(pl.Float32 if entry.width == 32 else pl.Float64)


class _LevelExport:
    """One level of one result, exported by one effective rule."""

    def __init__(self, rule: EffectiveRule, parsed: ParsedLevels, abundance: str | None) -> None:
        self._rule = rule
        self._level = parsed.levels[rule.level]
        self._sources = LevelSources(parsed, rule.level, abundance)
        self._provenance: dict[str, Provenance] = {}

    def has_index(self) -> bool:
        """Whether the result carries the source of the level's variable index."""
        key = self._rule.rule.axis.var_keys[0]
        (entry,) = [entry for entry in self._rule.rule.columns.var if entry.name == key]
        return (
            not isinstance(entry, SourcedEntry) or self._sources.resolve(entry.source) is not None
        )

    def parsed_level(self) -> ParsedLevel:
        """APB2's samples-by-features level for APB2's writer; key columns stay columns."""
        rule = self._rule.rule
        primary = rule.measurements.primary_layer
        found = [
            (entry, self._layer(entry, primary=entry.name == primary))
            for entry in rule.measurements.layers
        ]
        bound = self._bind_columns()
        uns = self._bind_group("uns", rule.uns, _UNS_TAKES)
        obs = self._columns(bound["obs"], ObservationRows(self._level.obs.frame.height))
        rows = FeatureRows(self._level.var.frame.height)
        var = self._columns(bound["var"], rows)
        kept = self._kept(var, rule.columns.var, rows.count())
        var = [column.filter(kept) for column in var]
        obs_key, var_key = rule.axis.obs_keys[0], rule.axis.var_keys[0]
        self._key(obs, obs_key, "obs")
        self._key(var, var_key, "var")
        layers = {
            entry.name: _layer_frame(values.values.filter(kept), entry)
            for entry, values in found
            if values is not None
        }
        return ParsedLevel.build(
            obs=pl.DataFrame(obs),
            obs_keys=[obs_key],
            var=pl.DataFrame(var),
            var_keys=[var_key],
            var_roles={},
            primary_layer=primary,
            abundance={primary: layers.pop(primary)},
            auxiliary=layers,
            uns=dict(self._level.uns),
            metadata={**self._level.metadata, "export": self._export_part(len(layers) + 1)},
            export_uns=cast(dict[str, JsonValue], self._uns(uns)),
        )

    def _layer(self, entry: LayerEntry, *, primary: bool) -> LayerValues | None:
        """Resolve a layer; ``None`` for an absent optional one. A caller's abundance layer
        replaces the rule's source of the primary layer."""
        chosen = primary and self._sources.abundance_chosen
        resolved = self._sources.resolve(_CALLERS_LAYER if chosen else entry.source)
        key = f"layers.{entry.name}"
        if resolved is None:
            if primary or entry.required:
                raise ValueError(f"layer {entry.name!r}: the result carries no layer for it")
            self._provenance[key] = {"location": "absent"}
            return None
        if isinstance(resolved, VarValues) and "catalog" in resolved.provenance:
            resolved = _repeated(resolved, self._level.obs.frame.height)
        if not isinstance(resolved, LayerValues):
            raise ValueError(f"layer {entry.name!r} cannot take {_describe(resolved.provenance)}")
        self._provenance[key] = resolved.provenance
        return resolved

    def _bind_columns(self) -> dict[str, list[_Bound | AllValues]]:
        rule = self._rule.rule
        return {
            "obs": self._bind_group("obs", rule.columns.obs, _OBS_TAKES),
            "var": self._bind_group("var", rule.columns.var, _VAR_TAKES),
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

    def _kept(self, columns: list[pl.Series], entries: Entries, count: int) -> pl.Series:
        """Features holding a value in every ``drop_missing`` column; the others are counted.

        A feature is one of the ``count`` variables written.
        """
        dropping = {e.name for e in entries if not isinstance(e, AllValues) and e.drop_missing}
        kept = pl.Series([True] * count)
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
    def _key(columns: list[pl.Series], key: str, group: str) -> pl.Series:
        (values,) = [column for column in columns if column.name == key]
        if values.is_duplicated().any():
            raise ValueError(f"{group} index {key!r} repeats values")
        return values

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

    def _export_part(self, matrices: int) -> dict[str, JsonValue]:
        """The export's own record of one level: where each entry came from, and a summary."""
        absent = sum(source.get("location") == "absent" for source in self._provenance.values())
        return {
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
        }


def _labels(values: pl.Series) -> pl.Series:
    """Values as the labels links match on; null, NaN and empty text match nothing."""
    if values.dtype.is_float():
        values = values.fill_nan(None)
    return values.cast(pl.String).replace("", None)


class FeatureLinks:
    """One feature link of the written levels, as a directed feature relation.

    Every written level with a var column ``column`` maps each feature to the target level's
    feature whose key equals that value; a value matching no key is counted as unlinked.
    """

    def __init__(self, link: FeatureLink, levels: Mapping[str, ParsedLevel]) -> None:
        """Match each level's ``column`` against the target level's keys.

        ``coordinates`` is ``None`` when the target level is not written.
        """
        self._name = link.name
        self.coordinates: pl.DataFrame | None = None
        self._unlinked = 0
        if link.target not in levels:
            return
        target = levels[link.target].var
        keys = _labels(target.frame.get_column(target.key_columns[0])).to_frame("label")
        keys = keys.with_row_index("column")
        blocks: list[pl.DataFrame] = []
        for name, level in levels.items():
            if name == link.target or link.column not in level.var.frame.columns:
                continue
            values = _labels(level.var.frame.get_column(link.column)).to_frame("label")
            values = values.with_row_index("row").drop_nulls("label")
            matched = values.join(keys, on="label", how="left", maintain_order="left")
            self._unlinked += matched.get_column("column").null_count()
            blocks.append(
                matched.drop_nulls("column").select(
                    source=pl.lit(name),
                    row="row",
                    target=pl.lit(link.target),
                    column="column",
                    value=pl.lit(1.0),
                )
            )
        self.coordinates = pl.concat(blocks).cast(_LINKS) if blocks else pl.DataFrame(schema=_LINKS)

    def counts(self) -> JsonValue:
        """The export record's entry: features linked and values matching no target feature."""
        linked = 0 if self.coordinates is None else self.coordinates.height
        return {"name": self._name, "linked": linked, "unlinked": self._unlinked}


class CompiledExport:
    """One target's effective rules, checked once and applied to any number of results."""

    def __init__(self, rules: Sequence[EffectiveRule]) -> None:
        """Keep one target's effective rules.

        Raises:
            ValueError: No level.
        """
        if not rules:
            raise ValueError("an export rule needs at least one level")
        self.rules = tuple(rules)
        self.output = rules[0].output
        self._format = output_for(self.extension)

    def write(self, exported: ParsedLevels, output: Path) -> None:
        """Write the already exported levels through the selected output format."""
        self._format.write(exported, output)

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
        """Whether every export writes one level: a one-level format, or a rule of one level."""
        return self._format.one_level or len(self.rules) == 1

    @property
    def extension(self) -> str:
        """The target's file suffix: ``.h5ad``, ``.h5mu``, ``.csv``, or ``.qpx``."""
        return self.output.extensions[0]

    def export(
        self, parsed: ParsedLevels, abundance: Mapping[str, str] | None = None
    ) -> ParsedLevels:
        """Build the target's export from one APB2 result.

        Args:
            parsed: An APB2 result holding the rule's required levels.
            abundance: Per level, the layer that becomes X; APB2's primary layer otherwise.

        Returns:
            A ``ParsedLevels`` for APB2's writer, one level per level written, and one feature
            relation per feature link.

        Raises:
            ValueError: The result lacks a required level or source, a source cannot be taken
                where the rule puts it, the levels share no observations, or ``abundance`` names
                a level the rule does not read.
            UnresolvedField: apb-catalog finds an answer ambiguous or the rule unreviewed.
        """
        chosen = dict(abundance or {})
        unknown = sorted(set(chosen) - set(self.written))
        if unknown:
            raise ValueError(f"abundance names levels the rule does not write: {unknown}")
        rules = self._exported_by(parsed)
        aligned = self._aligned(parsed, rules)
        exports, skipped = self._held(aligned, rules, chosen)
        dataset, dataset_sources = self._dataset_uns(aligned, exports[0][0].level)
        levels = {rule.rule.modality or rule.level: level.parsed_level() for rule, level in exports}
        record = self._export_record(dataset_sources, skipped)
        relations = {link.name: FeatureLinks(link, levels) for link in self.rules[0].feature_links}
        if relations:
            record["feature_links"] = [links.counts() for links in relations.values()]
        result = ParsedLevels(
            levels=levels,
            uns=dict(aligned.uns),
            metadata={
                **aligned.metadata,
                "export": {"schema_version": EXPORT_SCHEMA_VERSION, "provenance": record},
            },
            export_uns=cast(dict[str, JsonValue], dataset),
        )
        for name, links in relations.items():
            if links.coordinates is not None:
                result = result.with_feature_relation(
                    name, links.coordinates, {"producer": "apb-export"}
                )
        return result

    def _aligned(self, parsed: ParsedLevels, rules: Sequence[EffectiveRule]) -> ParsedLevels:
        """Several levels written, with the levels their layers reach into, on one observation
        axis, as APB2's writer stores them; one level written stays as it is.

        APB2 relabels levels keyed by experiment onto the raw files when the two map one to one,
        as in MaxQuant; the group holding the first level written is kept.
        """
        if len(rules) == 1:
            return parsed
        names = dict.fromkeys([*(rule.level for rule in rules), *self.levels])
        levels = {name: parsed.levels[name] for name in names if name in parsed.levels}
        return self._format.align(replace(parsed, levels=levels), rules[0].level)

    def _held(
        self, parsed: ParsedLevels, rules: Sequence[EffectiveRule], chosen: Mapping[str, str]
    ) -> tuple[list[tuple[EffectiveRule, _LevelExport]], list[str]]:
        """Each level's export; an optional level left out when the result carries no index for
        it, or when its observations cannot be aligned with the first level's."""
        exports: list[tuple[EffectiveRule, _LevelExport]] = []
        skipped: list[str] = []
        for rule in rules:
            if rule.level not in parsed.levels:
                if rule.rule.required:
                    raise ValueError(
                        f"level {rule.level!r} does not share the observations of "
                        f"{sorted(parsed.levels)}"
                    )
                logger.info(f"level {rule.level}: its observations are not aligned; left out")
                skipped.append(rule.level)
                continue
            level = _LevelExport(rule, parsed, chosen.get(rule.level))
            if not rule.rule.required and not level.has_index():
                logger.info(f"level {rule.level}: the result carries no index for it; left out")
                skipped.append(rule.level)
                continue
            exports.append((rule, level))
        if not exports:
            raise ValueError(f"no level of this result carries its index: {skipped}")
        return exports, skipped

    def _dataset_uns(
        self, parsed: ParsedLevels, level: str
    ) -> tuple[dict[str, object], dict[str, Provenance]]:
        """The table's own uns entries, rule fields read from the first level written."""
        sources = LevelSources(parsed, level, None)
        values: list[tuple[str, object]] = []
        provenance: dict[str, Provenance] = {}
        for entry in self.rules[0].dataset_uns:
            bound = _bind(entry, sources, "uns", _UNS_TAKES)
            provenance[f"uns.{entry.name}"] = (
                {"location": "absent"} if bound is None else bound.provenance
            )
            if bound is None:
                continue
            if isinstance(entry, Constant) and isinstance(entry.value, dict | list):
                values.append((entry.name, copy.deepcopy(entry.value)))
            else:
                values.append((entry.name, _series(bound, ScalarRow(), {}).item()))
        return _nest(values), provenance

    def _export_record(
        self, sources: dict[str, Provenance], skipped: list[str]
    ) -> dict[str, JsonValue]:
        """What wrote the container from which rule, its own uns sources and left-out levels."""
        rule = self.rules[0]
        record: dict[str, JsonValue] = {
            "package": "apb-export",
            "package_version": version("apb-export"),
            "apb2_version": version("apb2"),
            "target_name": rule.target_name,
            "target_version_pattern": rule.target_version_pattern,
            "export_rule": rule.document,
            "export_rule_file_version": rule.file_version,
        }
        if sources:
            record["sources"] = cast(JsonValue, sources)
        if skipped:
            record["skipped_levels"] = cast(JsonValue, skipped)
        return record

    def _exported_by(self, parsed: ParsedLevels) -> list[EffectiveRule]:
        """The levels that export this result: those for its software that it holds.

        A one-level format (``.h5ad``, ``.csv``) takes the first of them; the others take all
        of them and need every required one.
        """
        software = _software(parsed)
        rules = [
            rule
            for rule in self.rules
            if rule.rule.software is None or software & set(rule.rule.software)
        ]
        held = [rule for rule in rules if rule.level in parsed.levels]
        if self._format.one_level:
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
        return held[:1] if self._format.one_level else held
