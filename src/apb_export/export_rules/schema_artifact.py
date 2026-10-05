"""The published JSON Schemas of export rules, written from the models.

``python -m apb_export.export_rules.schema_artifact`` rewrites them after a model change.
"""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import BaseModel

from apb_export.export_rules.schema import ExportRuleDocument, LevelRule

SCHEMA_DIRECTORY = Path(__file__).parent / "documents" / "_schema"
_MODELS: dict[str, type[BaseModel]] = {
    "document.schema.json": ExportRuleDocument,
    "rule.schema.json": LevelRule,
}


def json_schemas() -> dict[str, dict[str, object]]:
    """The document and effective-rule schemas the models declare, by file name."""
    return {name: model.model_json_schema() for name, model in _MODELS.items()}


def write_artifacts() -> tuple[Path, ...]:
    """Write both schemas beside the packaged rules."""
    SCHEMA_DIRECTORY.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for name, schema in json_schemas().items():
        path = SCHEMA_DIRECTORY / name
        path.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")
        paths.append(path)
    return tuple(paths)


if __name__ == "__main__":
    for written in write_artifacts():
        print(written)
