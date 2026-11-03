"""The layout of export rule files: an object or list on one line where it fits, else one item
per line.

``python -m apb_export.export_rules.layout`` rewrites the packaged rules after an edit.
"""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import JsonValue

WIDTH = 100
"""The widest line an object or list stays on."""
_INDENT = 2
DOCUMENTS = Path(__file__).parent / "documents"


def format_rule(document: JsonValue) -> str:
    """A rule document as its file holds it."""
    return _render(document, 0, 0) + "\n"


def _render(value: JsonValue, indent: int, start: int) -> str:
    """One value whose first character sits in column ``start``, nested ``indent`` deep."""
    compact = json.dumps(value, ensure_ascii=False, separators=(", ", ": "))
    if not isinstance(value, dict | list) or not value or start + len(compact) < WIDTH:
        return compact
    pad = " " * (indent + _INDENT)
    if isinstance(value, dict):
        lines = []
        for key, item in value.items():
            prefix = f"{pad}{json.dumps(key, ensure_ascii=False)}: "
            lines.append(prefix + _render(item, indent + _INDENT, len(prefix)))
        return "{\n" + ",\n".join(lines) + "\n" + " " * indent + "}"
    items = [pad + _render(item, indent + _INDENT, len(pad)) for item in value]
    return "[\n" + ",\n".join(items) + "\n" + " " * indent + "]"


def packaged_rule_files() -> tuple[Path, ...]:
    """Every packaged ``<target>/<version>/rules.json``."""
    return tuple(sorted(DOCUMENTS.glob("*/*/rules.json")))


def write_layout(paths: tuple[Path, ...]) -> tuple[Path, ...]:
    """Rewrite each rule file in its layout; the files that changed."""
    changed: list[Path] = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        formatted = format_rule(json.loads(text))
        if formatted != text:
            path.write_text(formatted, encoding="utf-8")
            changed.append(path)
    return tuple(changed)


if __name__ == "__main__":
    for written in write_layout(packaged_rule_files()):
        print(written)
