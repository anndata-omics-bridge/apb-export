"""The REST server: stores uploads, queues jobs, serves the files jobs write.

No route computes an answer. ``POST /api/jobs`` writes the upload and a queued status; one
worker thread runs each job as its own process; every ``GET`` serves a file from disk.
"""

from __future__ import annotations

import os
import queue
import signal
import subprocess
import sys
import threading
from dataclasses import asdict
from pathlib import Path
from typing import Annotated, BinaryIO

from fastapi import FastAPI, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from apb_export.web.examples import Example
from apb_export.web.store import (
    TERMINAL,
    Job,
    JobStore,
    Output,
    Request,
    Status,
    safe_name,
    write_json,
)

STATIC = Path(__file__).parent / "static"
_JOB_FILES = frozenset({"status.json", "qc.json", "job.log"})
_CHUNK = 1 << 20


def fail(job: Job, message: str) -> None:
    status = job.status()
    status.state, status.error = "failed", message
    job.write_status(status)


class Worker:
    """One thread running queued jobs one at a time, each in its own process group."""

    __slots__ = ("_queue", "_thread", "running", "store", "timeout_seconds")

    def __init__(self, store: JobStore, timeout_seconds: float) -> None:
        self.store = store
        self.timeout_seconds = timeout_seconds
        self.running: str | None = None
        self._queue: queue.Queue[Job] = queue.Queue()
        self._thread = threading.Thread(target=self._loop, name="apb-export-web-jobs", daemon=True)
        self._thread.start()

    def submit(self, job: Job) -> None:
        self._queue.put(job)

    def join(self) -> None:
        """Wait until every submitted job has finished."""
        self._queue.join()

    def _loop(self) -> None:
        while True:
            job = self._queue.get()
            self.running = job.id
            try:
                self.run(job)
            finally:
                self.running = None
                self._queue.task_done()

    def run(self, job: Job) -> None:
        with job.log_path.open("a", encoding="utf-8") as log:
            process = subprocess.Popen(
                [sys.executable, "-m", "apb_export.web.job", str(job.root)],
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                code = process.wait(timeout=self.timeout_seconds)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                fail(job, f"timed out after {self.timeout_seconds:.0f} s")
                return
        if job.status().state not in TERMINAL:
            fail(job, f"the job process exited with code {code} before it finished")


def _save(upload: UploadFile, target: Path, budget: int) -> int:
    """Copy an upload to disk and return the bytes written; stop past ``budget``."""
    written = 0
    source: BinaryIO = upload.file
    with target.open("wb") as sink:
        while chunk := source.read(_CHUNK):
            written += len(chunk)
            if written > budget:
                raise HTTPException(413, "the upload exceeds this server's limit")
            sink.write(chunk)
    return written


class Api:
    """The route handlers, bound to one store, worker and offer."""

    __slots__ = ("examples", "max_upload_bytes", "offered", "software", "store", "worker")

    def __init__(
        self,
        store: JobStore,
        worker: Worker,
        outputs: list[Output],
        software: list[str],
        max_upload_bytes: int,
        examples: list[Example],
    ) -> None:
        self.store = store
        self.worker = worker
        self.offered = {output.id: output for output in outputs}
        self.software = software
        self.max_upload_bytes = max_upload_bytes
        self.examples = {example.id: example for example in examples}

    def _job(self, job_id: str) -> Job:
        job = self.store.get(job_id)
        if job is None:
            raise HTTPException(404, "no such job")
        return job

    def options(self) -> FileResponse:
        return FileResponse(self.store.options_path, headers={"Cache-Control": "no-store"})

    def _chosen(self, output: str, software: str, has_params: bool) -> Output:
        chosen = self.offered.get(output)
        if chosen is None:
            raise HTTPException(400, f"unknown output {output!r}")
        if software and software not in self.software:
            raise HTTPException(400, f"unknown software {software!r}")
        if not has_params and not software:
            raise HTTPException(400, "choose the software or upload its parameter file")
        return chosen

    def _save_all(self, job: Job, uploads: dict[str, list[UploadFile]]) -> dict[str, str]:
        """Each field's uploads under their own names, within the size limit.

        Several result files go into one folder named after the first, which the
        converter reads as one vendor result.
        """
        budget = self.max_upload_bytes
        saved: dict[str, str] = {}
        for field, files in uploads.items():
            if not files:
                continue
            names = [safe_name(upload.filename, field) for upload in files]
            if len(set(names)) != len(names):
                raise HTTPException(400, f"two {field} files share a name: {sorted(names)}")
            folder = job.inputs / field
            if len(files) > 1:
                folder = folder / Path(names[0]).stem
            folder.mkdir(parents=True)
            for upload, name in zip(files, names, strict=True):
                budget -= _save(upload, folder / name, budget)
            target = folder if len(files) > 1 else folder / names[0]
            saved[field] = str(target.relative_to(job.root))
        return saved

    def _example(self, example_id: str) -> Example:
        example = self.examples.get(example_id)
        if example is None:
            raise HTTPException(404, "no such example")
        return example

    def examples_file(self) -> FileResponse:
        return FileResponse(self.store.examples_path, headers={"Cache-Control": "no-store"})

    def example_file(self, example_id: str, name: str) -> FileResponse:
        path = self._example(example_id).files().get(name)
        if path is None:
            raise HTTPException(404, "no such file")
        return FileResponse(path, filename=name)

    def example_head(self, example_id: str, name: str) -> FileResponse:
        """The head of an example file, written at startup."""
        if name not in self._example(example_id).files():
            raise HTTPException(404, "no such file")
        path = self.store.previews / example_id / f"{name}.json"
        if not path.is_file():
            raise HTTPException(404, "no preview of this file")
        return FileResponse(path, headers={"Cache-Control": "no-cache"})

    def _link_example(self, job: Job, example: Example, fields: set[str]) -> dict[str, str]:
        """The chosen example files, linked into the job under the names the tool wrote.

        Several result files go into one folder named after the first, which the
        converter reads as one vendor result.
        """
        saved: dict[str, str] = {}
        for field, named in example.inputs().items():
            if field not in fields:
                continue
            folder = job.inputs / field
            if len(named) > 1:
                folder = folder / Path(named[0][0]).stem
            folder.mkdir(parents=True)
            for name, path in named:
                (folder / name).symlink_to(path, target_is_directory=path.is_dir())
            target = folder if len(named) > 1 else folder / named[0][0]
            saved[field] = str(target.relative_to(job.root))
        return saved

    def _inputs(
        self,
        job: Job,
        example: str,
        fields: set[str],
        uploads: dict[str, list[UploadFile]],
    ) -> dict[str, str]:
        if example:
            return self._link_example(job, self._example(example), fields)
        return self._save_all(job, uploads)

    def _example_fields(self, example: Example, example_files: str) -> set[str]:
        """Which of the example's files to use: the result, then parameters, then annotation."""
        fields = {field for field in example_files.split(",") if field} or {
            "data",
            "params",
            "annotation",
        }
        held = {"data"} | {
            field
            for field, path in (("params", example.params), ("annotation", example.annotation))
            if path is not None
        }
        if "data" not in fields or not fields <= {"data", "params", "annotation"}:
            raise HTTPException(400, "example files are data, then params, then annotation")
        if not fields <= held:
            raise HTTPException(400, f"this example has no {sorted(fields - held)}")
        return fields

    def submit(
        self,
        output: Annotated[str, Form()],
        data: list[UploadFile] | None = None,
        example: Annotated[str, Form()] = "",
        example_files: Annotated[str, Form()] = "",
        software: Annotated[str, Form()] = "",
        params: UploadFile | None = None,
        annotation: UploadFile | None = None,
    ) -> dict[str, str]:
        uploads = {
            "data": data or [],
            "params": [] if params is None else [params],
            "annotation": [] if annotation is None else [annotation],
        }
        fields: set[str] = set()
        if example:
            if any(uploads.values()):
                raise HTTPException(400, "an example brings its own files; upload none")
            chosen_example = self._example(example)
            fields = self._example_fields(chosen_example, example_files)
            # Read alone, a result names its own producer; with parameters, their software.
            named = (
                chosen_example.software if "params" in fields else chosen_example.result_software
            )
            software = named if named in self.software else ""
            has_params = "params" in fields
        elif not data:
            raise HTTPException(400, "upload a vendor result or choose an example")
        else:
            has_params = params is not None
        chosen = self._chosen(output, software, has_params)
        self.store.expire()
        job = self.store.create()
        try:
            saved = self._inputs(job, example, fields, uploads)
        except HTTPException:
            self.store.remove(job)
            raise
        request = Request(
            data=saved["data"],
            output=chosen,
            params=saved.get("params"),
            annotation=saved.get("annotation"),
            software=software or None,
        )
        write_json(job.request_path, asdict(request))
        job.write_status(Status(id=job.id, state="queued"))
        self.worker.submit(job)
        return {"id": job.id}

    def download(self, job_id: str, name: str) -> FileResponse:
        path = self._job(job_id).download(name)
        if path is None or not path.is_file():
            raise HTTPException(404, "no such file")
        return FileResponse(path, filename=name)

    def job_file(self, job_id: str, name: str) -> FileResponse:
        path = self._job(job_id).root / name
        if name not in _JOB_FILES or not path.is_file():
            raise HTTPException(404, "no such file")
        return FileResponse(path, headers={"Cache-Control": "no-store"})

    def remove(self, job_id: str) -> dict[str, str]:
        job = self._job(job_id)
        if self.worker.running == job.id or job.status().state == "queued":
            raise HTTPException(409, "the job has not finished")
        self.store.remove(job)
        return {"id": job_id}


def create_app(
    store: JobStore,
    worker: Worker,
    outputs: list[Output],
    software: list[str],
    max_upload_bytes: int,
    examples: list[Example] | None = None,
    static: Path = STATIC,
) -> FastAPI:
    """The API routes, then the built page at ``/``."""
    api = Api(store, worker, outputs, software, max_upload_bytes, examples or [])
    app = FastAPI(title="apb-export-web", docs_url=None, redoc_url=None)
    app.add_api_route("/api/options.json", api.options, methods=["GET"])
    app.add_api_route("/api/examples.json", api.examples_file, methods=["GET"])
    app.add_api_route("/api/examples/{example_id}/{name}", api.example_file, methods=["GET"])
    app.add_api_route("/api/examples/{example_id}/head/{name}", api.example_head, methods=["GET"])
    app.add_api_route("/api/jobs", api.submit, methods=["POST"])
    app.add_api_route("/api/jobs/{job_id}/files/{name}", api.download, methods=["GET"])
    app.add_api_route("/api/jobs/{job_id}/{name}", api.job_file, methods=["GET"])
    app.add_api_route("/api/jobs/{job_id}", api.remove, methods=["DELETE"])

    def page() -> FileResponse:
        # Revalidated on every load: a rebuild replaces the hashed scripts it names.
        return FileResponse(static / "index.html", headers={"Cache-Control": "no-cache"})

    app.add_api_route("/", page, methods=["GET"])
    app.mount("/assets", StaticFiles(directory=static / "assets"), name="assets")
    return app
