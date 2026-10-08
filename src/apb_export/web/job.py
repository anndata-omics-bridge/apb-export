"""``python -m apb_export.web.job <job-dir>``: convert, annotate, QC and package one upload.

Runs in its own process, one per job, so a crash or memory blow-up fails the job, not the
server. Conversion is the ``apb-export`` or ``apb2`` command itself; this module only orders
the commands, reads their results for QC, and writes ``status.json``, ``qc.json`` and the
downloads.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
import traceback
from collections.abc import Callable
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import polars as pl
from apb2.api import ParsedLevel, read_parsed_levels, sidecar_path

from apb_export.web.options import command
from apb_export.web.qc import Matrix, grouping_column, summarize
from apb_export.web.store import DownloadFile, Job, Request, Step, write_json


class StepFailed(Exception):
    """A command exited non-zero; the message is what it reported."""


def error_text(output: str) -> str:
    """The ERROR messages a command logged, else the tail of its output."""
    errors = [line.split(" - ", 1)[-1] for line in output.splitlines() if "| ERROR" in line]
    if errors:
        return "\n".join(errors)
    lines = [line for line in output.splitlines() if line.strip()]
    return "\n".join(lines[-5:]) or "the command failed without output"


class Runner:
    """One job's steps, each recorded in its status as it starts and ends."""

    __slots__ = ("job", "request", "status")

    def __init__(self, job: Job) -> None:
        self.job = job
        self.request = Request.read(job.request_path)
        self.status = job.status()

    def _path(self, relative: str | None) -> Path | None:
        return None if relative is None else self.job.root / relative

    def _step[T](self, name: str, action: Callable[[], T]) -> T:
        step = Step(name, "running")
        self.status.state = "running"
        self.status.steps.append(step)
        self.job.write_status(self.status)
        started = time.monotonic()
        try:
            result = action()
        except BaseException:
            step.state = "failed"
            step.seconds = round(time.monotonic() - started, 2)
            raise
        step.state = "done"
        step.seconds = round(time.monotonic() - started, 2)
        self.job.write_status(self.status)
        return result

    def _command(self, argv: list[str | Path]) -> None:
        with self.job.log_path.open("a", encoding="utf-8") as log:
            log.write(f"$ {' '.join(str(part) for part in argv)}\n")
        ran = subprocess.run(
            [str(part) for part in argv], capture_output=True, text=True, check=False
        )
        with self.job.log_path.open("a", encoding="utf-8") as log:
            log.write(ran.stdout + ran.stderr + "\n")
        if ran.returncode != 0:
            raise StepFailed(error_text(ran.stdout + ran.stderr))

    def _vendor_options(self) -> list[str | Path]:
        params, software = self._path(self.request.params), self.request.software
        return [
            *(["--params", params] if params is not None else []),
            *(["--software", software] if software is not None else []),
        ]

    def convert(self, results: Path) -> None:
        """Write the chosen output into ``results``."""
        output, data = self.request.output, self.job.root / self.request.data
        # A folder input such as maxquant_v2.8.1.0 has no suffix to drop.
        base = data.name if data.is_dir() else data.stem
        annotation = self._path(self.request.annotation)
        fasta = self._path(self.request.fasta)
        results.mkdir()
        if output.command == "apb-export":
            target = results / f"{base}_{output.name}{output.extension}"
            annotated = (
                ["--annotation", annotation] if annotation is not None and output.annotation else []
            )
            checked = ["--fasta", fasta] if fasta is not None and output.fasta else []
            self._step(
                "convert",
                lambda: self._command(
                    [
                        command("apb-export"),
                        output.name,
                        data,
                        target,
                        *self._vendor_options(),
                        *annotated,
                        *checked,
                    ]
                ),
            )
            return
        annotated_dir = results if fasta is None else self.job.root / "annotated"
        converted = annotated_dir if annotation is None else self.job.root / "converted"
        converted.mkdir(exist_ok=True)
        self._step(
            "convert",
            lambda: self._command(
                [
                    command("apb2"),
                    "convert",
                    data,
                    "--format",
                    output.name,
                    "--output",
                    converted / base,
                    *self._vendor_options(),
                ]
            ),
        )
        if annotation is not None:
            annotated_dir.mkdir(exist_ok=True)
            self._step("annotate", lambda: self._annotate(converted, annotation, annotated_dir))
        if fasta is not None:
            self._step("verify-peptides", lambda: self._verify(annotated_dir, fasta, results))

    def _annotate(self, converted: Path, annotation: Path, results: Path) -> None:
        for source in result_paths(converted):
            self._command([command("apb2"), "annotate", source, annotation, results / source.name])

    def _verify(self, annotated: Path, fasta: Path, results: Path) -> None:
        """Check each result's peptides against the FASTA, as ``apb-fasta verify-peptides``."""
        for source in result_paths(annotated):
            self._command(
                [
                    command("apb-fasta"),
                    "verify-peptides",
                    source,
                    fasta,
                    "--output",
                    results / source.name,
                ]
            )

    def run(self) -> None:
        results = self.job.root / "results"
        self.convert(results)
        self._step("qc", lambda: write_json(self.job.qc_path, self.qc(results)))
        self._step("package", lambda: self.package(results))
        self.status.state = "done"
        self.job.write_status(self.status)

    def qc(self, results: Path) -> dict[str, object]:
        columns = annotation_columns(self._path(self.request.annotation))
        matrices = [
            {"file": path.name, **summarize(matrix)}
            for path in result_paths(results)
            for matrix in (
                anndata_matrices(path, columns)
                if self.request.output.command == "apb-export"
                else parsed_matrices(path, columns)
            )
        ]
        return {"matrices": matrices}

    def package(self, results: Path) -> None:
        """Move every file to the downloads, zipping folder results such as Parquet."""
        downloads = self.job.downloads
        downloads.mkdir()
        for entry in sorted(results.iterdir()):
            if entry.is_dir():
                shutil.make_archive(
                    str(downloads / entry.name), "zip", root_dir=results, base_dir=entry.name
                )
            else:
                entry.replace(downloads / entry.name)
        self.status.files = [
            DownloadFile(path.name, path.stat().st_size) for path in sorted(downloads.iterdir())
        ]


def result_paths(folder: Path) -> list[Path]:
    """The results in a folder, without their APB sidecars."""
    entries = sorted(folder.iterdir())
    sidecars = {sidecar_path(entry) for entry in entries}
    return [entry for entry in entries if entry not in sidecars]


def annotation_columns(annotation: Path | None) -> list[str]:
    """The header of an uploaded CSV or TSV sample table, in file order."""
    if annotation is None:
        return []
    with annotation.open(encoding="utf-8-sig") as table:
        header = table.readline().rstrip("\r\n")
    separator = "\t" if "\t" in header else ","
    return [column.strip().strip('"') for column in header.split(separator)]


def _matrix(
    level: str,
    layer: str,
    samples: list[str],
    values: np.ndarray,
    sample_table: dict[str, list[str]],
    columns: list[str],
) -> Matrix:
    chosen = grouping_column(columns, sample_table)
    return Matrix(
        level=level,
        layer=layer,
        samples=samples,
        values=np.asarray(values, dtype=np.float64),
        groups=None if chosen is None else chosen[1],
        grouping=None if chosen is None else chosen[0],
    )


def _anndata_matrix(level: str, adata: ad.AnnData, columns: list[str]) -> Matrix:
    dense = adata.to_df().to_numpy()
    table = {str(column): adata.obs[column].astype(str).tolist() for column in adata.obs.columns}
    return _matrix(level, "X", [str(name) for name in adata.obs_names], dense.T, table, columns)


def anndata_matrices(path: Path, columns: list[str]) -> list[Matrix]:
    """X of an exported ``.h5ad``, or of each modality of an ``.h5mu``."""
    if path.suffix == ".h5mu":
        mdata = md.read_h5mu(path)
        return [
            _anndata_matrix(name, modality, columns)
            for name, modality in mdata.mod.items()
            if isinstance(modality, ad.AnnData)
        ]
    return [_anndata_matrix("X", ad.read_h5ad(path), columns)]


def _parsed_matrix(name: str, level: ParsedLevel, columns: list[str]) -> Matrix:
    obs = level.obs.frame.with_columns(pl.all().cast(pl.String).fill_null(""))
    samples = obs.select(pl.concat_str(level.obs.key_columns, separator=" | ")).to_series()
    values = level.layers[level.primary_layer_name].quantitative_values().to_numpy()
    table = {column: obs[column].to_list() for column in obs.columns}
    return _matrix(name, level.primary_layer_name, samples.to_list(), values, table, columns)


def parsed_matrices(path: Path, columns: list[str]) -> list[Matrix]:
    """The primary layer of every level of an APB2 result."""
    return [
        _parsed_matrix(name, level, columns)
        for name, level in read_parsed_levels(path).levels.items()
    ]


def main() -> None:
    """Run the job in ``sys.argv[1]`` and record its outcome; exit non-zero on failure."""
    job = Job(Path(sys.argv[1]))
    runner = Runner(job)
    try:
        runner.run()
    except StepFailed as error:
        runner.status.state, runner.status.error = "failed", str(error)
        job.write_status(runner.status)
        sys.exit(1)
    except Exception as error:
        # The job's boundary: record the failure for the page, keep the trace in the log.
        with job.log_path.open("a", encoding="utf-8") as log:
            log.write(traceback.format_exc())
        runner.status.state, runner.status.error = "failed", f"{type(error).__name__}: {error}"
        job.write_status(runner.status)
        sys.exit(1)


if __name__ == "__main__":
    main()
