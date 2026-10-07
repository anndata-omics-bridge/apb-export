"""``apb-export-web``: serve the upload-and-convert page on one host and port."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated

import uvicorn
from cyclopts import App, Parameter

from apb_export.web.examples import load_examples, write_examples
from apb_export.web.exports import load_hints, rule_versions
from apb_export.web.options import write_options
from apb_export.web.server import Worker, create_app
from apb_export.web.store import JobStore

app = App(
    name="apb-export-web",
    help="Serve a page that converts uploaded vendor results and shows their QC.",
)


@app.default
def serve(
    *,
    store: Annotated[Path, Parameter(help="Directory for options and job folders")] = Path(
        "apb-export-web"
    ),
    host: Annotated[str, Parameter(help="Interface to listen on")] = "127.0.0.1",
    port: Annotated[int, Parameter(help="Port to listen on")] = 8770,
    max_upload_mb: Annotated[int, Parameter(help="Largest upload, all files together")] = 4096,
    ttl_hours: Annotated[float, Parameter(help="Delete finished jobs after this long")] = 24.0,
    timeout_minutes: Annotated[float, Parameter(help="Fail a job running longer")] = 60.0,
    examples: Annotated[
        Path | None, Parameter(help="Corpus CSV whose rows the page offers as examples")
    ] = None,
    examples_root: Annotated[
        Path | None, Parameter(help="Data root the example corpus paths are relative to")
    ] = None,
    examples_sdrf: Annotated[
        Path | None, Parameter(help="Folder of <module>.sdrf.tsv files annotating the examples")
    ] = None,
) -> None:
    """Write the options, fail jobs a previous server left unfinished, and serve."""
    jobs = JobStore(store, ttl_seconds=ttl_hours * 3600)
    jobs.interrupt_unfinished()
    jobs.expire()
    max_upload_bytes = max_upload_mb * 1024 * 1024
    if (examples is None) != (examples_root is None):
        raise ValueError("--examples and --examples-root go together")
    software, outputs = write_options(jobs.options_path, max_upload_bytes)
    offered = (
        []
        if examples is None or examples_root is None
        else load_examples(
            examples, examples_root, examples_sdrf, software, load_hints()[1], rule_versions()
        )
    )
    write_examples(jobs.examples_path, offered)
    worker = Worker(jobs, timeout_seconds=timeout_minutes * 60)
    app = create_app(jobs, worker, outputs, software, max_upload_bytes, offered)
    uvicorn.run(app, host=host, port=port)


def main() -> None:
    """Run the command line and exit with its status."""
    sys.exit(app(sys.argv[1:]))
