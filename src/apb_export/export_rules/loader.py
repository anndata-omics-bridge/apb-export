"""Packaged export rules and their effective levels."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from importlib.resources.abc import Traversable

from pydantic import JsonValue

from apb_export.export_rules.schema import ExportRuleDocument, LevelRule, Output

_DOCUMENTS = resources.files("apb_export.export_rules").joinpath("documents")
_RULES_FILE = "rules.json"


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
    """

    document: str
    target_name: str
    target_version_pattern: str
    file_version: str
    output: Output
    level: str
    rule: LevelRule


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

    An ``.h5ad`` table's levels are alternatives in order: each before the last names the
    software it is for, so a software-specific level precedes the one every software takes.

    Raises:
        ValueError: The document has more than one table, an ``.h5ad`` level before the last
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
                )
            )
        unnamed = [rule.level for rule in rules[:-1] if rule.rule.software is None]
        if ".h5ad" in table.output.extensions and unnamed:
            raise ValueError(f"{name}: .h5ad levels {unnamed} precede another but name no software")
    return tuple(rules)


def _rule_files() -> dict[str, Traversable]:
    files: dict[str, Traversable] = {}
    for target in sorted(_DOCUMENTS.iterdir(), key=lambda entry: entry.name):
        if not target.is_dir() or target.name.startswith("_"):
            continue
        for version in sorted(target.iterdir(), key=lambda entry: entry.name):
            candidate = version.joinpath(_RULES_FILE)
            if candidate.is_file():
                files[f"{target.name}/{version.name}/{_RULES_FILE}"] = candidate
    return files


def packaged_documents() -> tuple[str, ...]:
    """Every packaged rule file, as ``<target>/<version>/rules.json``."""
    return tuple(_rule_files())


def packaged_rules(target: str) -> tuple[EffectiveRule, ...]:
    """The effective rules of one target's packaged rule file.

    Raises:
        ValueError: No rule, or more than one rule file, is packaged for the target.
    """
    names = [name for name in _rule_files() if name.split("/", 1)[0] == target]
    if len(names) != 1:
        raise ValueError(f"expected one packaged export rule for {target!r}, found {names}")
    text = _rule_files()[names[0]].read_text(encoding="utf-8")
    return effective_rules(ExportRuleDocument.model_validate_json(text), names[0])
