"""What the page offers: vendor software, outputs, and the upload limit; written once."""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

from apb2.api import get_rules

from apb_export.api import Exporter
from apb_export.web.exports import load_hints
from apb_export.web.store import Output, write_json

# apb2's --format values and what each writes.
_APB2_FORMATS = (
    ("hdf5", "APB2 result: h5ad or h5mu"),
    ("parquet", "APB2 result: Parquet folder (zip)"),
    ("duckdb", "APB2 result: DuckDB"),
)


def command(name: str) -> Path:
    """A console script installed beside this interpreter.

    Raises:
        FileNotFoundError: The script is not installed in this environment.
    """
    path = Path(sys.executable).with_name(name)
    if not path.is_file():
        raise FileNotFoundError(f"{name} is not installed beside {sys.executable}")
    return path


def software() -> list[str]:
    """Every software with a packaged quant-result rule, DDA and DIA together."""
    names = {rules.software_name for category in ("DDA", "DIA") for rules in get_rules(category)}
    return sorted(names, key=str.casefold)


def _accepts_annotation(target: str) -> bool:
    """Whether the target's subcommand takes ``--annotation``; its own help says so."""
    shown = subprocess.run(
        [command("apb-export"), target, "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    return "--annotation" in shown.stdout


def outputs() -> list[Output]:
    """apb-export's targets first, then APB2's own result formats."""
    targets = [
        Output(
            id=f"apb-export:{target}",
            label=f"{target} ({Exporter(target).extension})",
            command="apb-export",
            name=target,
            extension=Exporter(target).extension,
            annotation=_accepts_annotation(target),
        )
        for target in Exporter.targets()
    ]
    formats = [
        Output(
            id=f"apb2:{name}",
            label=label,
            command="apb2",
            name=name,
            extension="",
            annotation=True,
        )
        for name, label in _APB2_FORMATS
    ]
    return targets + formats


ABOUT = Path(__file__).parent / "outputs.json"


def about() -> dict[str, dict[str, object]]:
    """What each output gives, which tool opens it and how, with links, by output id."""
    return json.loads(ABOUT.read_text(encoding="utf-8"))


def write_options(path: Path, max_upload_bytes: int) -> tuple[list[str], list[Output]]:
    """Write ``options.json`` and return the software and outputs it lists.

    Raises:
        ValueError: An offered output has no description in ``outputs.json``.
    """
    names, offered = software(), outputs()
    described = about()
    missing = [output.id for output in offered if output.id not in described]
    if missing:
        raise ValueError(f"{ABOUT} describes no output {missing}")
    write_json(
        path,
        {
            "software": names,
            "hints": load_hints()[0],
            "outputs": [{**asdict(output), "about": described[output.id]} for output in offered],
            "max_upload_bytes": max_upload_bytes,
        },
    )
    return names, offered
