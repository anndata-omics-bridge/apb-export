"""Job directories: the request, the status, the inputs, the outputs, all as files."""

from __future__ import annotations

import json
import os
import re
import shutil
import time
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

type State = Literal["queued", "running", "done", "failed"]

_JOB_ID = re.compile(r"^[0-9a-f]{32}$")
TERMINAL: frozenset[State] = frozenset({"done", "failed"})


def write_json(path: Path, value: object) -> None:
    """Write JSON atomically, so a reader never sees half a file."""
    partial = path.with_name(f".{path.name}.partial")
    partial.write_text(json.dumps(value, indent=1), encoding="utf-8")
    partial.replace(path)


@dataclass(frozen=True, slots=True)
class Output:
    """One output the page offers: an apb-export target or an APB2 result format."""

    id: str
    label: str
    command: Literal["apb-export", "apb2"]
    # The target name for apb-export, the --format value for apb2.
    name: str
    # The target's file suffix; empty for apb2, which picks its own.
    extension: str
    annotation: bool
    # Whether a FASTA changes the result: msmu's contaminant flags, APB2's peptide check.
    fasta: bool = False


@dataclass(frozen=True, slots=True)
class Request:
    """What the user uploaded and chose; paths are relative to the job directory."""

    data: str
    output: Output
    params: str | None = None
    annotation: str | None = None
    software: str | None = None
    fasta: str | None = None

    @classmethod
    def read(cls, path: Path) -> Request:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            data=raw["data"],
            output=Output(**raw["output"]),
            params=raw["params"],
            annotation=raw["annotation"],
            software=raw["software"],
            fasta=raw.get("fasta"),
        )


@dataclass(slots=True)
class Step:
    name: str
    state: State
    seconds: float | None = None


@dataclass(slots=True)
class DownloadFile:
    name: str
    bytes: int


@dataclass(slots=True)
class Status:
    """The one file the page polls."""

    id: str
    state: State
    steps: list[Step] = field(default_factory=lambda: list[Step]())
    error: str | None = None
    files: list[DownloadFile] = field(default_factory=lambda: list[DownloadFile]())

    @classmethod
    def read(cls, path: Path) -> Status:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            id=raw["id"],
            state=raw["state"],
            steps=[Step(**step) for step in raw["steps"]],
            error=raw["error"],
            files=[DownloadFile(**item) for item in raw["files"]],
        )


class Job:
    """One job directory."""

    __slots__ = ("root",)

    def __init__(self, root: Path) -> None:
        self.root = root

    @property
    def id(self) -> str:
        return self.root.name

    @property
    def inputs(self) -> Path:
        return self.root / "input"

    @property
    def downloads(self) -> Path:
        return self.root / "downloads"

    @property
    def request_path(self) -> Path:
        return self.root / "request.json"

    @property
    def status_path(self) -> Path:
        return self.root / "status.json"

    @property
    def qc_path(self) -> Path:
        return self.root / "qc.json"

    @property
    def log_path(self) -> Path:
        return self.root / "job.log"

    def status(self) -> Status:
        return Status.read(self.status_path)

    def write_status(self, status: Status) -> None:
        write_json(self.status_path, asdict(status))

    def download(self, name: str) -> Path | None:
        """A finished file by its listed name; ``None`` for anything else."""
        listed = {item.name for item in self.status().files}
        return self.downloads / name if name in listed else None


class JobStore:
    """Every job directory under one root, with a time-to-live for finished jobs."""

    __slots__ = ("jobs", "root", "ttl_seconds")

    def __init__(self, root: Path, ttl_seconds: float) -> None:
        self.root = root
        self.jobs = root / "jobs"
        self.ttl_seconds = ttl_seconds
        self.jobs.mkdir(parents=True, exist_ok=True)

    @property
    def options_path(self) -> Path:
        return self.root / "options.json"

    @property
    def examples_path(self) -> Path:
        return self.root / "examples.json"

    @property
    def previews(self) -> Path:
        """One folder per example, holding ``<name>.json``, the head of each of its files."""
        return self.root / "previews"

    def create(self) -> Job:
        job = Job(self.jobs / uuid.uuid4().hex)
        job.inputs.mkdir(parents=True)
        return job

    def get(self, job_id: str) -> Job | None:
        """The job, when the id is well formed and its status exists."""
        if not _JOB_ID.match(job_id):
            return None
        job = Job(self.jobs / job_id)
        return job if job.status_path.is_file() else None

    def remove(self, job: Job) -> None:
        shutil.rmtree(job.root)

    def expire(self) -> None:
        """Delete finished jobs older than the time-to-live, and debris without a status."""
        cutoff = time.time() - self.ttl_seconds
        for root in self.jobs.iterdir():
            job = Job(root)
            if not job.status_path.is_file():
                if root.stat().st_mtime < cutoff:
                    shutil.rmtree(root)
            elif job.status().state in TERMINAL and job.status_path.stat().st_mtime < cutoff:
                shutil.rmtree(root)

    def interrupt_unfinished(self) -> None:
        """Fail every job a previous server left queued or running; nothing will resume it."""
        for root in self.jobs.iterdir():
            job = Job(root)
            if job.status_path.is_file():
                status = job.status()
                if status.state not in TERMINAL:
                    status.state = "failed"
                    status.error = "interrupted: the server restarted before this job finished"
                    job.write_status(status)


def safe_name(name: str | None, fallback: str) -> str:
    """An uploaded file's own name, without any directory part; vendor detection reads it."""
    base = os.path.basename((name or "").replace("\\", "/")).strip()
    return base if base not in {"", ".", ".."} else fallback
