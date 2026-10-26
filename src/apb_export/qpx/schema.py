"""The QPX 1.1 Arrow contract, frozen from the qpx 1.1.6 schemas."""

from __future__ import annotations

import json
from importlib.resources import files
from typing import Any, cast

import pyarrow as pa

_SPEC = cast(dict[str, Any], json.loads(files(__package__).joinpath("schema.json").read_text()))


def _type(name: str) -> pa.DataType:
    if name.startswith("list<"):
        return pa.list_(_type(name[5:-1]))
    if name in _SPEC["types"]:
        return pa.struct(_fields(_SPEC["types"][name]["fields"]))
    return pa.type_for_alias(name)


def _fields(fields: dict[str, Any]) -> list[pa.Field]:
    return [
        pa.field(name, _type(field["type"]), nullable=field.get("nullable", True))
        for name, field in fields.items()
        if not field.get("optional", False)
    ]


def arrow_schema(view: str) -> pa.Schema:
    """Required columns and their exact nested types; unknown values may be null."""
    return pa.schema(_fields(_SPEC[view]["fields"]))
