"""Which files each software writes, as its own documentation names them.

``hints.json`` lists, per software, the exports the converter reads: the native result file
names, the parameter file, and the apb2 rule (and so the version range) each export belongs
to. The page shows them as upload hints; examples find their own export by version, file
kind and corpus label, so they can show and serve the names the tool wrote.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from apb2.api import get_rules

HINTS = Path(__file__).parent / "hints.json"


def _pattern(name: str) -> bool:
    """A name the user or run chooses, such as "<report name>.tsv" or "*.txt"."""
    return any(mark in name for mark in "*<>")


@dataclass(frozen=True, slots=True)
class Export:
    """One kind of output a software writes and the converter reads."""

    # The apb2 rules that read it, such as "alphadia/v1_10/rules.json"; their version
    # patterns say which versions write it.
    rules: tuple[str, ...]
    # The versions it covers, for people: "1.10".
    versions: str
    # The result file or files, as the tool names them.
    result: tuple[str, ...]
    # The parameter or log file, as the tool names it; None when there is none.
    params: str | None
    # How it arrives: file suffixes such as ".tsv", or "folder" for a directory.
    kinds: tuple[str, ...]
    # Corpus names that choose it whatever the version, such as "FragPipe (DIA-NN quant)".
    labels: tuple[str, ...] = ()
    note: str | None = None

    def concrete_name(self) -> str | None:
        """The file name an example of this export goes by; None for a pattern or several files."""
        if len(self.result) != 1 or _pattern(self.result[0]):
            return None
        return Path(self.result[0]).name

    def result_names(self, count: int) -> list[str] | None:
        """Names for ``count`` result files in the tool's order; None when they are not fixed."""
        if count == 1:
            name = self.concrete_name()
            return None if name is None else [name]
        if count != len(self.result) or any(_pattern(name) for name in self.result):
            return None
        return [Path(name).name for name in self.result]

    def concrete_params(self) -> str | None:
        """The parameter file name an example goes by; None for a pattern."""
        if self.params is None or _pattern(self.params):
            return None
        return Path(self.params).name

    def describe(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class SoftwareHint:
    exports: tuple[Export, ...]
    params_required: bool
    sources: tuple[str, ...]


def load_hints(path: Path = HINTS) -> tuple[dict[str, object], dict[str, SoftwareHint]]:
    """The raw document for the page, and the parsed hints by software."""
    document = json.loads(path.read_text(encoding="utf-8"))
    hints = {
        name: SoftwareHint(
            exports=tuple(
                Export(
                    rules=tuple(export["rules"]),
                    versions=export["versions"],
                    result=tuple(export["result"]),
                    params=export.get("params"),
                    kinds=tuple(export["kinds"]),
                    labels=tuple(export.get("labels", ())),
                    note=export.get("note"),
                )
                for export in entry["exports"]
            ),
            params_required=entry.get("params_required", False),
            sources=tuple(entry.get("sources", ())),
        )
        for name, entry in document["software"].items()
    }
    return document, hints


def rule_versions() -> dict[str, str]:
    """Each packaged apb2 rule and the version pattern it accepts."""
    return {
        variant.rule: variant.software_version_pattern
        for category in ("DDA", "DIA")
        for rules in get_rules(category)
        for variant in rules.variants
    }


def export_for(
    hint: SoftwareHint,
    patterns: dict[str, str],
    label: str,
    version: str | None,
    data: Path,
) -> Export | None:
    """The one export an example is: by corpus label, else by version, then by file kind.

    Without a version, as for an example without a parameter file, the file kind alone must
    single one out, as it does for pb_custom and MetaMorpheus.
    """
    kind = "folder" if data.is_dir() else data.suffix.lower()
    labelled = [export for export in hint.exports if label in export.labels]
    candidates = labelled or [
        export
        for export in hint.exports
        if not export.labels
        and (
            version is None
            or any(re.match(patterns.get(rule, r"(?!)"), version.strip()) for rule in export.rules)
        )
    ]
    matching = [export for export in candidates if kind in export.kinds]
    return matching[0] if len(matching) == 1 else None
