"""Collect the web app's examples into one folder a deployed ``apb-export-web`` mounts.

The folder holds ``examples.csv``, each example's files under ``data/`` at its corpus path,
and the module SDRFs under ``sdrf/``. Its rows are a corpus's, such as APB Studio's routine
corpus, plus extra rows. An extra row is a corpus row under the same root, such as a
ProteoBench submission outside the routine corpus, unless its ``source`` column names a file
to copy from the workspace, such as a sample apb2 commits; a ``.gz`` source is unpacked. Its
``sdrf`` column, when set, names an SDRF placed beside the input as ``sdrf.tsv``, for a tool
that named its samples otherwise than the module SDRF's raw files, such as WOMBAT.
Serve it with
``--examples <out>/examples.csv --examples-root <out>/data --examples-sdrf <out>/sdrf``.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import sys
from pathlib import Path

from cyclopts import App

app = App(name="web_examples", help=__doc__)

COLUMNS = ("input_file", "vendor_parameter_file", "module", "software_name")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _copy(source: Path, target: Path) -> None:
    """Copy one file, unpacking a ``.gz`` source whose target is not one."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".gz" and target.suffix != ".gz":
        with gzip.open(source, "rb") as packed, target.open("wb") as plain:
            shutil.copyfileobj(packed, plain)
    else:
        shutil.copy2(source, target)


@app.default
def collect(
    corpus: Path,
    root: Path,
    sdrf: Path,
    out: Path,
    *,
    extra: Path | None = None,
    workspace: Path = Path(".."),
) -> None:
    """Copy the corpus rows' folders, the extra rows' sources and the SDRFs into ``out``."""
    rows = _rows(corpus)
    added = [] if extra is None else _rows(extra)
    folders: set[Path] = set()
    for row in [*rows, *(row for row in added if not row["source"])]:
        for column in ("input_file", "vendor_parameter_file"):
            if row[column]:
                path = root / row[column]
                folders.add(path if path.is_dir() else path.parent)
    for folder in sorted(folders):
        shutil.copytree(folder, out / "data" / folder.relative_to(root), dirs_exist_ok=True)
    for row in added:
        if row["source"]:
            _copy(workspace / row["source"], out / "data" / row["input_file"])
        if row["sdrf"]:
            _copy(workspace / row["sdrf"], (out / "data" / row["input_file"]).parent / "sdrf.tsv")
    with (out / "examples.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows([*rows, *added])
    (out / "sdrf").mkdir(exist_ok=True)
    for path in sorted(sdrf.glob("*.sdrf.tsv")):
        shutil.copy2(path, out / "sdrf" / path.name)
    print(f"{len(rows)} corpus and {len(added)} extra examples -> {out}")


if __name__ == "__main__":
    sys.exit(app(sys.argv[1:]))
