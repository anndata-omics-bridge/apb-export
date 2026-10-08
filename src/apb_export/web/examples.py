"""Example inputs served beside the page: the rows of one corpus CSV, as APB Studio writes it.

A row names an input file or folder, its parameter file, its ProteoBench module and its
software, relative to a data root. An ``sdrf.tsv`` inside an input folder, or beside an input
file, is the example's sample annotation; otherwise ``<module>.sdrf.tsv`` in an SDRF folder,
such as apb-proteobench's module data, is.

ProteoBench stores uploads as ``input_file.*`` and ``param_0.*``. An example's version, read
from its parameter file, picks the export it is, so the page shows and serves the names the
tool itself wrote.
"""

from __future__ import annotations

import csv
import hashlib
import re
import shutil
from dataclasses import dataclass
from itertools import islice
from pathlib import Path

import polars as pl
from apb2.api import parse_search_parameters
from loguru import logger

from apb_export.web.exports import Export, SoftwareHint, export_for
from apb_export.web.store import write_json

_COLUMNS = ("input_file", "vendor_parameter_file", "module", "software_name")
# What a preview shows: the first lines, each cut to a width a page can still scroll.
HEAD_LINES = 20
_HEAD_WIDTH = 4000


@dataclass(frozen=True, slots=True)
class Example:
    id: str
    # The corpus's name, such as "FragPipe (DIA-NN quant)".
    label: str
    # The software whose parameter file it carries, by apb2's name: "FragPipe".
    software: str
    # The software that wrote the result, read alone without parameters: "DIA-NN".
    result_software: str
    module: str
    data: Path
    params: Path | None
    annotation: Path | None
    # The version its parameter file names; None without one.
    version: str | None = None
    # The export it is, when its version and file kind identify one.
    export: Export | None = None
    # Further result files stored beside ``data`` as ``<stem>_*``, such as AlphaDIA 1.12's
    # precursor.matrix.tsv; APB Studio passes them with it as one folder.
    secondary: tuple[Path, ...] = ()
    # The FASTA its module prescribes, from the corpus's optional ``fasta`` column.
    fasta: Path | None = None

    def inputs(self) -> dict[str, list[tuple[str, Path]]]:
        """Each input field's files: the names the tool wrote, else the stored names."""
        data = [self.data, *self.secondary]
        native = (
            self.export.result_names(len(data)) if self.export and self.data.is_file() else None
        )
        named = {"data": list(zip(native or [path.name for path in data], data, strict=True))}
        if self.params is not None:
            name = self.export.concrete_params() if self.export else None
            named["params"] = [(name or self.params.name, self.params)]
        if self.annotation is not None:
            named["annotation"] = [(self.annotation.name, self.annotation)]
        if self.fasta is not None:
            named["fasta"] = [(self.fasta.name, self.fasta)]
        return named

    def folder_files(self) -> list[Path]:
        """The files inside a folder result, such as MaxQuant's txt folder; none for a file."""
        if not self.data.is_dir():
            return []
        return sorted(path for path in self.data.iterdir() if path.is_file())

    def files(self) -> dict[str, Path]:
        """The example's downloadable files by the name they are served under.

        A folder result's files go by their own names; one that is also the parameter file or
        annotation, such as MaxQuant's mqpar.xml, is the same file.
        """
        named = {path.name: path for path in self.folder_files()}
        named.update(
            (name, path)
            for files in self.inputs().values()
            for name, path in files
            if path.is_file()
        )
        return named

    def describe(self) -> dict[str, object]:
        inputs = self.inputs()
        params = inputs.get("params")
        annotation = inputs.get("annotation")
        return {
            "id": self.id,
            "label": self.label,
            "software": self.software,
            "version": self.version,
            "module": self.module,
            "data": [name for name, _ in inputs["data"]],
            "stored_data": [path.name for _, path in inputs["data"]],
            "folder": self.data.is_dir(),
            "folder_files": [path.name for path in self.folder_files()],
            "params": None if params is None else params[0][0],
            "stored_params": None if self.params is None else self.params.name,
            "annotation": None if annotation is None else annotation[0][0],
            "fasta": None if self.fasta is None else self.fasta.name,
            "export": None if self.export is None else self.export.describe(),
            "bytes": sum(_size(path) for path in (self.data, *self.secondary)),
        }


def _size(path: Path) -> int:
    if path.is_dir():
        return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
    return path.stat().st_size


def _secondary(data: Path) -> tuple[Path, ...]:
    """Files stored beside a result as ``<stem>_*``, as APB Studio finds them."""
    if data.is_dir():
        return ()
    return tuple(path for path in sorted(data.parent.glob(f"{data.stem}_*")) if path.is_file())


def _annotation(data: Path, module: str, sdrf: Path | None) -> Path | None:
    candidates = [(data if data.is_dir() else data.parent) / "sdrf.tsv"]
    if sdrf is not None:
        candidates.append(sdrf / f"{module}.sdrf.tsv")
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def _key(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _names(label: str, software: list[str]) -> tuple[str, str]:
    """The parameter and result software of a corpus name, as apb2 names them.

    A compound name such as "FragPipe (DIA-NN quant)" carries FragPipe's parameters and
    DIA-NN's result, as APB Studio reads it.
    """
    by_key = {_key(name): name for name in software}
    parameters = label.split("(", 1)[0].strip()
    compound = re.search(r"\((.+?)\s+quant\)", label)
    result = compound.group(1) if compound else parameters
    return by_key.get(_key(parameters), parameters), by_key.get(_key(result), result)


def _version(params: Path | None, software: str) -> str | None:
    """The version a parameter file names, as apb2 reads it; None when it cannot."""
    if params is None:
        return None
    try:
        version = parse_search_parameters(params, software).software_version
    except (OSError, ValueError) as error:
        # Display only: the example still runs; the page shows the stored name instead.
        logger.warning(f"no version for example parameters {params}: {error}")
        return None
    return None if version is None else str(version).strip()


def load_examples(
    corpus: Path,
    root: Path,
    sdrf: Path | None = None,
    software: list[str] | None = None,
    hints: dict[str, SoftwareHint] | None = None,
    patterns: dict[str, str] | None = None,
) -> list[Example]:
    """Every row of the corpus CSV whose input exists under ``root``.

    With hints, a row no export of its software describes stays in the corpus but is no
    example: the page offers only what it says the converter reads. AlphaDIA 1.10's rows are
    such, a table ProteoBench joined from AlphaDIA's two outputs.

    Paths are absolute: jobs link to them from their own folders.

    Raises:
        ValueError: The CSV lacks a corpus column, or names an input that does not exist.
    """
    with corpus.open(encoding="utf-8", newline="") as table:
        rows = list(csv.DictReader(table))
    root = root.absolute()
    sdrf = None if sdrf is None else sdrf.absolute()
    examples: list[Example] = []
    for row in rows:
        absent = [column for column in _COLUMNS if column not in row]
        if absent:
            raise ValueError(f"{corpus} lacks columns {absent}")
        data = root / row["input_file"]
        if not data.exists():
            raise ValueError(f"{corpus} names {data}, which does not exist")
        params = root / row["vendor_parameter_file"] if row["vendor_parameter_file"] else None
        parameter_software, result_software = _names(row["software_name"], software or [])
        hint = (hints or {}).get(parameter_software)
        version = _version(params, parameter_software) if hint is not None else None
        export = (
            None
            if hint is None
            else export_for(hint, patterns or {}, row["software_name"], version, data)
        )
        if hint is not None and export is None:
            logger.warning(f"{row['input_file']}: no {parameter_software} hint describes it")
            continue
        examples.append(
            Example(
                id=hashlib.sha256(row["input_file"].encode()).hexdigest()[:12],
                label=row["software_name"],
                software=parameter_software,
                result_software=result_software,
                module=row["module"],
                data=data,
                params=params,
                annotation=_annotation(data, row["module"], sdrf),
                version=version,
                export=export,
                secondary=_secondary(data),
                fasta=root / row["fasta"] if row.get("fasta") else None,
            )
        )
    return examples


def write_examples(path: Path, examples: list[Example]) -> None:
    write_json(path, {"examples": [example.describe() for example in examples]})


def head(path: Path) -> dict[str, object]:
    """The first lines of a file as the page previews it; a Parquet file's first rows as TSV."""
    if path.suffix.lower() == ".parquet":
        table = pl.read_parquet(path, n_rows=HEAD_LINES)
        lines = table.write_csv(separator="\t").splitlines()
        complete = table.height < HEAD_LINES
        kind = "parquet"
    else:
        with path.open(encoding="utf-8", errors="replace", newline="") as text:
            lines = [line.rstrip("\r\n") for line in islice(text, HEAD_LINES + 1)]
        complete = len(lines) <= HEAD_LINES
        lines, kind = lines[:HEAD_LINES], "text"
    return {
        "name": path.name,
        "bytes": path.stat().st_size,
        "format": kind,
        "lines": [line[:_HEAD_WIDTH] for line in lines],
        "complete": complete,
    }


def write_previews(folder: Path, examples: list[Example]) -> None:
    """Each example file's head as ``<folder>/<example id>/<served name>.json``."""
    if folder.exists():
        shutil.rmtree(folder)
    for example in examples:
        (folder / example.id).mkdir(parents=True)
        for name, path in example.files().items():
            write_json(folder / example.id / f"{name}.json", {**head(path), "name": name})
