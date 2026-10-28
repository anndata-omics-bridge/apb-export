"""Packaged and plugin export rules and their effective levels."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import cache
from importlib import resources
from importlib.metadata import EntryPoints, entry_points
from importlib.resources.abc import Traversable

from pydantic import JsonValue

from apb_export.export_rules.schema import (
    ColumnEntry,
    ExportRuleDocument,
    FeatureLink,
    LevelRule,
    Output,
)

_DOCUMENTS = resources.files("apb_export.export_rules").joinpath("documents")
_RULES_FILE = "rules.json"

_ONE_LEVEL = (".h5ad", ".csv")
"""Outputs that hold one level, the first of the table's levels that a result holds."""

PLUGIN_GROUP = "apb_export.targets"
"""The entry-point group through which installed packages add export targets."""


@dataclass(frozen=True, slots=True)
class EffectiveRule:
    """One table's level with its base merged in: what one export writes.

    Attributes:
        document: The rule file, such as ``msmu/v0_4/rules.json``.
        target_name: The tool the export is read by.
        target_version_pattern: Target releases the rule was checked against.
        file_version: The rule file's own version.
        output: What the table writes.
        level: The APB2 level it reads.
        rule: The level's validated declarations.
        dataset_uns: The table's own uns entries, for the whole container.
        feature_links: The table's relations between modalities.
    """

    document: str
    target_name: str
    target_version_pattern: str
    file_version: str
    output: Output
    level: str
    rule: LevelRule
    dataset_uns: tuple[ColumnEntry, ...] = ()
    feature_links: tuple[FeatureLink, ...] = ()


def _merge(base: dict[str, JsonValue], level: dict[str, JsonValue]) -> dict[str, JsonValue]:
    """The level overrides base keys; named ``accession_syntax`` blocks combine."""
    merged = {**base, **level}
    syntax: dict[str, JsonValue] = {}
    for block in (base.get("accession_syntax"), level.get("accession_syntax")):
        if isinstance(block, dict):
            syntax.update(block)
    if syntax:
        merged["accession_syntax"] = syntax
    return merged


def effective_rules(document: ExportRuleDocument, name: str) -> tuple[EffectiveRule, ...]:
    """Validate every table level of a rule document.

    An ``.h5ad`` or ``.csv`` table's levels are alternatives in order: each before the last
    names the software it is for, so a software-specific level precedes the one every software
    takes.

    Raises:
        ValueError: The document has more than one table, such a level before the last
            names no software, a level's declarations are invalid, or an ``.h5mu`` level lacks
            a modality.
    """
    if len(document.tables) != 1:
        raise ValueError(f"{name}: an export rule writes one table, not {len(document.tables)}")
    rules: list[EffectiveRule] = []
    for table in document.tables:
        for level, block in table.levels.items():
            rule = LevelRule.model_validate(_merge(table.base, block))
            if ".h5mu" in table.output.extensions and rule.modality is None:
                raise ValueError(f"{name}: level {level!r} writes .h5mu but names no modality")
            rules.append(
                EffectiveRule(
                    document=name,
                    target_name=document.target_name,
                    target_version_pattern=document.target_version_pattern,
                    file_version=document.file_version,
                    output=table.output,
                    level=level,
                    rule=rule,
                    dataset_uns=tuple(table.uns),
                    feature_links=tuple(table.feature_links),
                )
            )
        unnamed = [rule.level for rule in rules[:-1] if rule.rule.software is None]
        one_level = [ext for ext in table.output.extensions if ext in _ONE_LEVEL]
        if one_level and unnamed:
            raise ValueError(
                f"{name}: {one_level[0]} levels {unnamed} precede another but name no software"
            )
    return tuple(rules)


def _versions(target: str, folder: Traversable) -> dict[str, Traversable]:
    """Each ``<version>/rules.json`` under one target's folder, keyed as its document name."""
    files: dict[str, Traversable] = {}
    for version in sorted(folder.iterdir(), key=lambda entry: entry.name):
        candidate = version.joinpath(_RULES_FILE)
        if version.is_dir() and candidate.is_file():
            files[f"{target}/{version.name}/{_RULES_FILE}"] = candidate
    return files


def plugin_folders(points: EntryPoints) -> dict[str, Traversable]:
    """The rule folder each entry point names, by target.

    Raises:
        ValueError: Two entry points name one target.
    """
    repeated = sorted(
        name for name, count in Counter(point.name for point in points).items() if count > 1
    )
    if repeated:
        raise ValueError(f"several installed plugins name the export targets {repeated}")
    return {point.name: resources.files(point.value) for point in points}


@cache
def _plugin_folders() -> dict[str, Traversable]:
    """Every installed plugin's rule folder, by target."""
    return plugin_folders(entry_points(group=PLUGIN_GROUP))


def _rule_files() -> dict[str, Traversable]:
    """Packaged rule files, then those of installed plugins.

    Raises:
        ValueError: A plugin target has a packaged target's name.
    """
    files: dict[str, Traversable] = {}
    for target in sorted(_DOCUMENTS.iterdir(), key=lambda entry: entry.name):
        if target.is_dir() and not target.name.startswith("_"):
            files |= _versions(target.name, target)
    plugins = _plugin_folders()
    clashing = sorted({name.split("/", 1)[0] for name in files} & set(plugins))
    if clashing:
        raise ValueError(f"plugin export targets {clashing} have the names of packaged targets")
    for target, folder in sorted(plugins.items()):
        files |= _versions(target, folder)
    return files


def packaged_documents() -> tuple[str, ...]:
    """Every packaged and plugin rule file, as ``<target>/<version>/rules.json``."""
    return tuple(_rule_files())


def packaged_rules(target: str) -> tuple[EffectiveRule, ...]:
    """The effective rules of one target's packaged or plugin rule file.

    Raises:
        ValueError: No rule, or more than one rule file, is packaged or installed for the target.
    """
    names = [name for name in _rule_files() if name.split("/", 1)[0] == target]
    if len(names) != 1:
        raise ValueError(f"expected one export rule for {target!r}, found {names}")
    text = _rule_files()[names[0]].read_text(encoding="utf-8")
    return effective_rules(ExportRuleDocument.model_validate_json(text), names[0])
