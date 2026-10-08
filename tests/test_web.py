"""apb-export-web: upload, job, QC and downloads, on synthetic DIA-NN input."""

from __future__ import annotations

import json
import re
import zipfile
from collections.abc import Iterator
from dataclasses import asdict
from io import BytesIO
from pathlib import Path

import httpx
import mudata as md
import numpy as np
import pandas as pd
import polars as pl
import pytest
from fastapi.testclient import TestClient

from apb_export.web.examples import (
    HEAD_LINES,
    Example,
    head,
    load_examples,
    write_examples,
    write_previews,
)
from apb_export.web.exports import Export, SoftwareHint, export_for, load_hints, rule_versions
from apb_export.web.job import error_text
from apb_export.web.options import about, write_options
from apb_export.web.qc import ALL_SAMPLES, Matrix, grouping_column, summarize
from apb_export.web.server import Worker, create_app
from apb_export.web.store import JobStore, Output, Status, safe_name, write_json
from conftest import RUNS, DiannInput, write_diann

type Options = tuple[list[str], list[Output]]


@pytest.fixture(scope="module")
def options(tmp_path_factory: pytest.TempPathFactory) -> Options:
    """What the server offers; computed once, since it asks each target's --help."""
    return write_options(tmp_path_factory.mktemp("options") / "options.json", 1 << 30)


@pytest.fixture
def client(tmp_path: Path, options: Options) -> Iterator[tuple[TestClient, Worker, JobStore]]:
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    write_json(store.options_path, {"software": software, "outputs": [asdict(o) for o in outputs]})
    worker = Worker(store, timeout_seconds=600)
    with TestClient(create_app(store, worker, outputs, software, 1 << 30)) as test_client:
        yield test_client, worker, store


def _annotation(folder: Path) -> Path:
    """A prolfquapp-style sample table: the run key, then one condition per run."""
    table = folder / "samples.tsv"
    rows = "".join(f"{run}\t{run[4]}\n" for run in RUNS)
    table.write_text("raw_file\tcondition\n" + rows, encoding="utf-8")
    return table


def _fasta(folder: Path) -> Path:
    """Two proteins: PEPTIDEK's, and a ProteoBench-style contaminant holding CONTPEPK."""
    fasta = folder / "proteins.fasta"
    fasta.write_text(
        ">sp|P1|ONE_HUMAN One OS=Homo sapiens OX=9606 GN=ONE\nMPEPTIDEKAACLLKR\n"
        ">sp|Cont_P4|FOUR_BOVIN Four OS=Bos taurus OX=9913 GN=FOUR\nMCONTPEPKR\n",
        encoding="utf-8",
    )
    return fasta


def _post(
    test_client: TestClient,
    diann: DiannInput,
    output: str,
    annotation: Path | None = None,
    fasta: Path | None = None,
) -> httpx.Response:
    files = {
        "data": (diann.report.name, diann.report.read_bytes()),
        "params": (diann.log.name, diann.log.read_bytes()),
    }
    for field, path in (("annotation", annotation), ("fasta", fasta)):
        if path is not None:
            files[field] = (path.name, path.read_bytes())
    return test_client.post("/api/jobs", data={"output": output}, files=files)


def _submit(
    test_client: TestClient,
    diann: DiannInput,
    output: str,
    annotation: Path | None = None,
    fasta: Path | None = None,
) -> str:
    response = _post(test_client, diann, output, annotation, fasta)
    assert response.status_code == 200, response.text
    return response.json()["id"]


def test_options_list_targets_formats_and_annotation_support(options: Options) -> None:
    software, outputs = options
    by_id = {output.id: output for output in outputs}

    assert "DIA-NN" in software
    assert by_id["apb-export:prolfqua"].extension == ".h5ad"
    assert by_id["apb-export:prolfqua"].annotation
    assert not by_id["apb-export:msmu"].annotation, "msmu's command takes no --annotation"
    assert {"apb2:hdf5", "apb2:parquet", "apb2:duckdb"} <= by_id.keys()
    assert {output for output, offered in by_id.items() if offered.fasta} == {
        "apb-export:msmu",
        "apb2:hdf5",
        "apb2:parquet",
        "apb2:duckdb",
    }, "only msmu reads the FASTA check; APB2 results keep it"


def test_an_msmu_job_flags_the_uploaded_fastas_contaminants(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput, tmp_path: Path
) -> None:
    test_client, worker, _ = client
    job = _submit(test_client, diann, "apb-export:msmu", fasta=_fasta(tmp_path))
    worker.join()

    assert test_client.get(f"/api/jobs/{job}/status.json").json()["state"] == "done"
    content = test_client.get(f"/api/jobs/{job}/files/report_msmu.h5mu").content
    target = tmp_path / "out.h5mu"
    target.write_bytes(content)
    var = md.read_h5mu(target)["psm"].var
    assert isinstance(var, pd.DataFrame)
    assert var.loc["run_A1.CONTPEPK/2", "contaminant"] == 1


def test_an_apb2_job_checks_the_peptides_against_the_uploaded_fasta(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput, tmp_path: Path
) -> None:
    test_client, worker, _ = client
    job = _submit(test_client, diann, "apb2:hdf5", _annotation(tmp_path), _fasta(tmp_path))
    worker.join()

    status = test_client.get(f"/api/jobs/{job}/status.json").json()
    assert status["state"] == "done", status
    assert [step["name"] for step in status["steps"]] == [
        "convert",
        "annotate",
        "verify-peptides",
        "qc",
        "package",
    ]


def test_an_output_that_reads_no_fasta_refuses_one(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput, tmp_path: Path
) -> None:
    test_client, _, store = client
    response = _post(test_client, diann, "apb-export:prolfqua", fasta=_fasta(tmp_path))

    assert response.status_code == 400
    assert "takes no FASTA" in response.text
    assert list(store.jobs.iterdir()) == []


def test_prolfqua_job_groups_cv_by_the_annotation(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput, tmp_path: Path
) -> None:
    test_client, worker, _ = client
    job = _submit(test_client, diann, "apb-export:prolfqua", _annotation(tmp_path))
    worker.join()

    status = test_client.get(f"/api/jobs/{job}/status.json").json()
    assert status["state"] == "done", status
    assert [step["name"] for step in status["steps"]] == ["convert", "qc", "package"]
    names = [item["name"] for item in status["files"]]
    assert names == ["report_prolfqua.h5ad", "report_prolfqua.h5ad.apb.json"]
    qc = test_client.get(f"/api/jobs/{job}/qc.json").json()
    (matrix,) = qc["matrices"]
    assert matrix["cv"]["grouping"] == "condition"
    assert [series["group"] for series in matrix["cv"]["series"]] == ["A", "B"]
    assert [sample["group"] for sample in matrix["samples"]] == ["A", "A", "B", "B"]
    download = test_client.get(f"/api/jobs/{job}/files/report_prolfqua.h5ad")
    assert download.status_code == 200
    assert download.content[:8] == b"\x89HDF\r\n\x1a\n"


def test_msmu_ignores_the_annotation_and_uses_all_samples(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput, tmp_path: Path
) -> None:
    test_client, worker, _ = client
    job = _submit(test_client, diann, "apb-export:msmu", _annotation(tmp_path))
    worker.join()

    assert test_client.get(f"/api/jobs/{job}/status.json").json()["state"] == "done"
    qc = test_client.get(f"/api/jobs/{job}/qc.json").json()
    assert {matrix["cv"]["grouping"] for matrix in qc["matrices"]} == {ALL_SAMPLES}


def test_apb2_parquet_is_annotated_and_zipped(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput, tmp_path: Path
) -> None:
    test_client, worker, _ = client
    job = _submit(test_client, diann, "apb2:parquet", _annotation(tmp_path))
    worker.join()

    status = test_client.get(f"/api/jobs/{job}/status.json").json()
    assert status["state"] == "done", status
    assert [step["name"] for step in status["steps"]] == ["convert", "annotate", "qc", "package"]
    archive = next(item["name"] for item in status["files"] if item["name"].endswith(".zip"))
    content = test_client.get(f"/api/jobs/{job}/files/{archive}").content
    assert any(
        name.endswith("manifest.json") for name in zipfile.ZipFile(BytesIO(content)).namelist()
    )
    qc = test_client.get(f"/api/jobs/{job}/qc.json").json()
    assert {matrix["cv"]["grouping"] for matrix in qc["matrices"]} == {"condition"}


def test_a_failing_conversion_reports_the_commands_error(
    client: tuple[TestClient, Worker, JobStore], tmp_path: Path
) -> None:
    test_client, worker, _ = client
    data = tmp_path / "nonsense.tsv"
    data.write_text("a\tb\n1\t2\n", encoding="utf-8")
    response = test_client.post(
        "/api/jobs",
        data={"output": "apb2:duckdb", "software": "DIA-NN"},
        files={"data": (data.name, data.read_bytes())},
    )
    job = response.json()["id"]
    worker.join()

    status = test_client.get(f"/api/jobs/{job}/status.json").json()
    assert status["state"] == "failed"
    assert status["error"]
    assert status["files"] == []
    assert test_client.get(f"/api/jobs/{job}/job.log").status_code == 200


def test_uploads_are_validated_before_a_job_exists(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput
) -> None:
    test_client, _, store = client
    data = {"data": (diann.report.name, diann.report.read_bytes())}

    assert test_client.post("/api/jobs", data={"output": "nope"}, files=data).status_code == 400
    no_params = test_client.post("/api/jobs", data={"output": "apb2:hdf5"}, files=data)
    assert no_params.status_code == 400
    unknown = {"output": "apb2:hdf5", "software": "Excel"}
    assert test_client.post("/api/jobs", data=unknown, files=data).status_code == 400
    assert list(store.jobs.iterdir()) == []


def test_an_oversized_upload_leaves_no_job(
    tmp_path: Path, options: Options, diann: DiannInput
) -> None:
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    worker = Worker(store, timeout_seconds=60)
    with TestClient(create_app(store, worker, outputs, software, 10)) as test_client:
        response = test_client.post(
            "/api/jobs",
            data={"output": "apb2:hdf5", "software": "DIA-NN"},
            files={"data": (diann.report.name, diann.report.read_bytes())},
        )
    assert response.status_code == 413
    assert list(store.jobs.iterdir()) == []


def test_files_are_served_only_by_listed_name(
    client: tuple[TestClient, Worker, JobStore], diann: DiannInput
) -> None:
    test_client, worker, _ = client
    job = _submit(test_client, diann, "apb-export:proteopy")
    worker.join()

    assert test_client.get(f"/api/jobs/{job}/request.json").status_code == 404
    assert test_client.get(f"/api/jobs/{job}/files/..%2Fstatus.json").status_code == 404
    assert test_client.get("/api/jobs/not-a-job/status.json").status_code == 404
    assert test_client.delete(f"/api/jobs/{job}").status_code == 200
    assert test_client.get(f"/api/jobs/{job}/status.json").status_code == 404


def test_the_page_and_options_are_served(client: tuple[TestClient, Worker, JobStore]) -> None:
    test_client, _, _ = client
    assert "DIA-NN" in test_client.get("/api/options.json").json()["software"]
    page = test_client.get("/")
    assert page.status_code == 200
    assert "<apb-upload-form>" in page.text
    assert page.headers["cache-control"] == "no-cache", "a rebuild must never strand an old page"
    script = re.search(r'src="\./(assets/[^"]+\.js)"', page.text)
    assert script is not None
    assert test_client.get(f"/{script.group(1)}").status_code == 200


def test_unfinished_jobs_fail_when_a_new_server_starts(tmp_path: Path) -> None:
    store = JobStore(tmp_path, ttl_seconds=0)
    job = store.create()
    job.write_status(Status(id=job.id, state="running"))
    store.interrupt_unfinished()
    assert job.status().state == "failed"
    store.expire()
    assert not job.root.exists(), "finished jobs past their time-to-live are deleted"


def test_qc_counts_densities_and_cv_by_group() -> None:
    values = np.array(
        [
            [100.0, 110.0, 1000.0, np.nan],
            [200.0, 0.0, 2000.0, 2200.0],
            [300.0, 330.0, np.nan, np.nan],
        ]
    )
    matrix = Matrix(
        "ion", "Intensity", ["a1", "a2", "b1", "b2"], values, ["A", "A", "B", "B"], "condition"
    )

    summary = summarize(matrix)

    assert summary["features"] == 3
    samples = summary["samples"]
    assert isinstance(samples, list)
    assert [sample["detected"] for sample in samples] == [3, 2, 2, 1]
    cv = summary["cv"]
    assert isinstance(cv, dict)
    a, b = cv["series"]
    assert (a["group"], a["features"], b["features"]) == ("A", 2, 1), "zero and NaN are missing"
    assert a["median"] == pytest.approx(100 * np.std([100, 110], ddof=1) / 105, abs=1e-3)
    assert sum(a["counts"]) == 2


def test_cv_falls_back_to_all_samples() -> None:
    matrix = Matrix("ion", "Intensity", ["s1", "s2"], np.array([[1.0, 3.0]]))
    cv = summarize(matrix)["cv"]
    assert isinstance(cv, dict)
    assert cv["grouping"] == ALL_SAMPLES
    assert cv["series"][0]["features"] == 1


def test_grouping_needs_a_column_that_splits_samples() -> None:
    table = {"raw_file": ["r1", "r2", "r3"], "condition": ["A", "A", "B"]}
    assert grouping_column(["raw_file", "condition"], table) == ("condition", ["A", "A", "B"])
    assert grouping_column(["raw_file", "absent"], table) is None
    constant = {**table, "disease": ["not available"] * 3}
    assert grouping_column(["disease", "condition"], constant) == ("condition", ["A", "A", "B"])
    written = {**table, "factor_value_spike": ["x", "y", "y"]}
    chosen = grouping_column(["condition", "factor value[spike]"], written)
    assert chosen == ("factor_value_spike", ["x", "y", "y"]), "matched as apb2 writes it"


def test_upload_names_lose_their_directories() -> None:
    assert safe_name("../../etc/passwd", "data") == "passwd"
    assert safe_name("C:\\runs\\report.tsv", "data") == "report.tsv"
    assert safe_name("..", "data") == "data"
    assert safe_name(None, "params") == "params"


def test_error_text_prefers_logged_errors() -> None:
    logged = "x | INFO | m - fine\ny | ERROR | m:f:1 - no rule matches\n"
    assert error_text(logged) == "no rule matches"
    assert error_text("") == "the command failed without output"


def _corpus(root: Path) -> Path:
    """A one-row corpus CSV: a synthetic DIA-NN 1.8 run stored as ProteoBench stores uploads."""
    folder = root / "submissions" / "one"
    folder.mkdir(parents=True)
    diann = write_diann(folder)
    diann.report.rename(folder / "input_file.tsv")
    diann.log.rename(folder / "param_0..txt")
    _annotation(folder).rename(folder / "sdrf.tsv")
    corpus = root / "corpus.csv"
    corpus.write_text(
        "input_file,vendor_parameter_file,module,software_name\n"
        "submissions/one/input_file.tsv,submissions/one/param_0..txt,dia_test,DIA-NN\n",
        encoding="utf-8",
    )
    return corpus


def _load(corpus: Path, root: Path, sdrf: Path | None = None) -> list[Example]:
    """The corpus's examples, named as the page names them."""
    return load_examples(
        corpus, root, sdrf, ["DIA-NN", "FragPipe"], load_hints()[1], rule_versions()
    )


@pytest.fixture
def examples(tmp_path: Path) -> list[Example]:
    root = tmp_path / "data"
    return _load(_corpus(root), root)


def test_an_example_brings_the_fasta_its_corpus_row_names(tmp_path: Path) -> None:
    root = tmp_path / "data"
    (example,) = _load(_with_fasta(_corpus(root), root), root)

    assert example.fasta == root / "module.fasta"
    assert example.describe()["fasta"] == "module.fasta"
    assert example.inputs()["fasta"] == [("module.fasta", root / "module.fasta")]


def test_examples_read_a_corpus_and_find_its_sdrf(examples: list[Example]) -> None:
    (example,) = examples
    described = example.describe()
    keys = ("software", "version", "module", "data", "stored_data", "params", "annotation")
    assert {key: described[key] for key in keys} == {
        "software": "DIA-NN",
        "version": "1.8.1",
        "module": "dia_test",
        "data": ["report.tsv"],
        "stored_data": ["input_file.tsv"],
        "params": "report.log.txt",
        "annotation": "sdrf.tsv",
    }, "the names DIA-NN wrote, not ProteoBench's storage names"
    assert set(example.files()) == {"report.tsv", "report.log.txt", "sdrf.tsv"}


def test_examples_under_a_relative_root_have_absolute_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    corpus = _corpus(tmp_path / "data")
    monkeypatch.chdir(tmp_path)

    (example,) = _load(corpus, Path("data"))

    assert example.data.is_absolute(), "a job's link to a relative path would dangle"
    assert example.annotation is not None and example.annotation.is_absolute()


def test_a_corpus_naming_a_missing_input_is_refused(tmp_path: Path) -> None:
    corpus = tmp_path / "corpus.csv"
    corpus.write_text(
        "input_file,vendor_parameter_file,module,software_name\nabsent.tsv,,m,Sage\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="does not exist"):
        load_examples(corpus, tmp_path)


def test_an_example_runs_without_upload_and_brings_its_annotation(
    tmp_path: Path, options: Options, examples: list[Example]
) -> None:
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    write_examples(store.examples_path, examples)
    worker = Worker(store, timeout_seconds=600)
    (example,) = examples
    with TestClient(create_app(store, worker, outputs, software, 1 << 30, examples)) as test_client:
        listed = test_client.get("/api/examples.json").json()["examples"]
        assert [item["id"] for item in listed] == [example.id]
        sdrf = test_client.get(f"/api/examples/{example.id}/sdrf.tsv")
        assert sdrf.text.startswith("raw_file\tcondition")
        assert test_client.get(f"/api/examples/{example.id}/corpus.csv").status_code == 404

        response = test_client.post(
            "/api/jobs", data={"output": "apb-export:prolfqua", "example": example.id}
        )
        job = response.json()["id"]
        worker.join()
        status = test_client.get(f"/api/jobs/{job}/status.json").json()
        assert status["state"] == "done", status
        qc = test_client.get(f"/api/jobs/{job}/qc.json").json()
        assert {matrix["cv"]["grouping"] for matrix in qc["matrices"]} == {"condition"}
        assert (store.jobs / job / "input" / "data" / "report.tsv").is_symlink()

        mixed = test_client.post(
            "/api/jobs",
            data={"output": "apb2:hdf5", "example": example.id},
            files={"data": ("x.tsv", b"a\n")},
        )
        assert mixed.status_code == 400
        unknown = test_client.post("/api/jobs", data={"output": "apb2:hdf5", "example": "nope"})
        assert unknown.status_code == 404
        nothing = test_client.post("/api/jobs", data={"output": "apb2:hdf5"})
        assert nothing.status_code == 400


def test_a_module_sdrf_annotates_examples_and_groups_by_its_factor(
    tmp_path: Path, options: Options
) -> None:
    root = tmp_path / "data"
    corpus = _corpus(root)
    (root / "submissions" / "one" / "sdrf.tsv").unlink()
    modules = tmp_path / "modules"
    modules.mkdir()
    rows = "".join(f"{run}\t{run}.raw\t{run[4]}\n" for run in RUNS)
    (modules / "dia_test.sdrf.tsv").write_text(
        "source name\tcomment[data file]\tfactor value[condition]\n" + rows, encoding="utf-8"
    )
    (example,) = _load(corpus, root, modules)
    assert example.annotation == modules / "dia_test.sdrf.tsv"

    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    worker = Worker(store, timeout_seconds=600)
    with TestClient(
        create_app(store, worker, outputs, software, 1 << 30, [example])
    ) as test_client:
        job = test_client.post(
            "/api/jobs", data={"output": "apb-export:prolfqua", "example": example.id}
        ).json()["id"]
        worker.join()
        status = test_client.get(f"/api/jobs/{job}/status.json").json()
        assert status["state"] == "done", status
        qc = test_client.get(f"/api/jobs/{job}/qc.json").json()
    (matrix,) = qc["matrices"]
    assert matrix["cv"]["grouping"] == "factor_value_condition"
    assert [series["group"] for series in matrix["cv"]["series"]] == ["A", "B"]


def test_hints_name_every_offered_software_and_packaged_rule(options: Options) -> None:
    software, _ = options
    _, hints = load_hints()
    patterns = rule_versions()
    named = {rule for hint in hints.values() for export in hint.exports for rule in export.rules}
    # An export may name a rule apb2 is renaming; one that names no packaged rule is stale.
    stale = [
        export.result
        for hint in hints.values()
        for export in hint.exports
        if not set(export.rules) & patterns.keys()
    ]

    assert set(software) <= set(hints)
    assert stale == [], f"exports read by no packaged apb2 rule: {stale}"
    # AlphaDIA 1.10's rule reads only a table ProteoBench joined; AlphaDIA never writes it.
    unoffered = {"alphadia/v1_10/rules.json"}
    missing = patterns.keys() - named - unoffered
    assert missing == set(), f"apb2 rules without a hint: {missing}"
    assert not unoffered & named, "the page claims a rule it leaves out"
    assert hints["Sage"].params_required, "apb2 refuses Sage without its parameter file"
    assert all(hint.sources for hint in hints.values()), (
        "every hint cites where its names come from"
    )


def _with_fasta(corpus: Path, root: Path) -> Path:
    """The corpus with its one row naming a module FASTA under ``root``."""
    _fasta(root).rename(root / "module.fasta")
    rows = corpus.read_text(encoding="utf-8").splitlines()
    corpus.write_text(f"{rows[0]},fasta\n{rows[1]},module.fasta\n", encoding="utf-8")
    return corpus


def _run_example(tmp_path: Path, options: Options, files: str) -> tuple[JobStore, str]:
    root = tmp_path / "data"
    (example,) = _load(_with_fasta(_corpus(root), root), root)
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    worker = Worker(store, timeout_seconds=600)
    with TestClient(
        create_app(store, worker, outputs, software, 1 << 30, [example])
    ) as test_client:
        form = {"output": "apb2:hdf5", "example": example.id, "example_files": files}
        job = test_client.post("/api/jobs", data=form).json()["id"]
        worker.join()
    return store, job


def test_a_result_only_example_sends_the_result_and_the_software(
    tmp_path: Path, options: Options
) -> None:
    store, job = _run_example(tmp_path, options, "data")

    request = json.loads((store.jobs / job / "request.json").read_text(encoding="utf-8"))
    assert (request["params"], request["annotation"], request["software"]) == (None, None, "DIA-NN")
    command = (store.jobs / job / "job.log").read_text(encoding="utf-8").splitlines()[0]
    assert "--software DIA-NN" in command
    assert "--params" not in command


@pytest.mark.parametrize(
    ("files", "steps", "grouping"),
    [
        ("data,params", ["convert", "qc", "package"], ALL_SAMPLES),
        ("data,params,annotation", ["convert", "annotate", "qc", "package"], "condition"),
        (
            "data,params,annotation,fasta",
            ["convert", "annotate", "verify-peptides", "qc", "package"],
            "condition",
        ),
    ],
)
def test_an_example_adds_params_then_annotation(
    tmp_path: Path, options: Options, files: str, steps: list[str], grouping: str
) -> None:
    store, job = _run_example(tmp_path, options, files)

    finished = store.get(job)
    assert finished is not None
    status = finished.status()
    assert status.state == "done", status
    assert [step.name for step in status.steps] == steps
    qc = json.loads(finished.qc_path.read_text(encoding="utf-8"))
    assert {matrix["cv"]["grouping"] for matrix in qc["matrices"]} == {grouping}


def test_example_file_choices_must_start_with_the_result(tmp_path: Path, options: Options) -> None:
    root = tmp_path / "data"
    corpus = _corpus(root)
    (root / "submissions" / "one" / "sdrf.tsv").unlink()
    (example,) = _load(corpus, root)
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    worker = Worker(store, timeout_seconds=60)
    with TestClient(
        create_app(store, worker, outputs, software, 1 << 30, [example])
    ) as test_client:
        for files in ("params", "data,fasta", "data,params,annotation"):
            form = {"output": "apb2:hdf5", "example": example.id, "example_files": files}
            assert test_client.post("/api/jobs", data=form).status_code == 400, files
    assert list(store.jobs.iterdir()) == []


def test_a_compound_corpus_name_splits_into_parameter_and_result_software(
    tmp_path: Path,
) -> None:
    root = tmp_path / "data"
    corpus = _corpus(root)
    corpus.write_text(
        corpus.read_text(encoding="utf-8").replace(",DIA-NN\n", ",FragPipe (DIA-NN quant)\n"),
        encoding="utf-8",
    )
    (example,) = load_examples(corpus, root, software=["DIA-NN", "FragPipe"])
    assert (example.label, example.software, example.result_software) == (
        "FragPipe (DIA-NN quant)",
        "FragPipe",
        "DIA-NN",
    )


def test_a_two_file_export_names_its_files_in_the_tools_order() -> None:
    export = Export(
        rules=("alphadia/v1_12/rules.json",),
        versions="1.12",
        result=("precursors.tsv", "precursor.matrix.tsv"),
        params="log.txt",
        kinds=(".tsv",),
    )

    assert export.result_names(2) == ["precursors.tsv", "precursor.matrix.tsv"]
    assert export.result_names(1) is None, "one file of a two-file export has no fixed name"


def test_without_a_version_the_file_kind_alone_must_pick_the_export(tmp_path: Path) -> None:
    hints, patterns = load_hints()[1], rule_versions()
    table = tmp_path / "table.tsv"
    table.write_text("Sequence\n", encoding="utf-8")

    custom = export_for(hints["pb_custom"], patterns, "pb_custom", None, table)
    peaks = export_for(hints["MetaMorpheus"], patterns, "MetaMorpheus", None, table)
    ambiguous = export_for(hints["DIA-NN"], patterns, "DIA-NN", None, tmp_path / "report.parquet")

    assert custom is not None and custom.rules == ("pb_custom/rules.json",)
    assert peaks is not None and peaks.result == ("Task1-SearchTask/AllQuantifiedPeaks.tsv",)
    assert ambiguous is None, "DIA-NN 1.9 and 2.x both write report.parquet"


def test_a_folder_example_lists_and_serves_every_file_in_it(tmp_path: Path) -> None:
    folder = tmp_path / "data" / "txt"
    folder.mkdir(parents=True)
    for name in ("evidence.txt", "peptides.txt", "mqpar.xml"):
        (folder / name).write_text(f"{name}\n", encoding="utf-8")
    corpus = tmp_path / "corpus.csv"
    corpus.write_text(
        "input_file,vendor_parameter_file,module,software_name\ntxt,txt/mqpar.xml,m,MaxQuant\n",
        encoding="utf-8",
    )

    (example,) = load_examples(corpus, tmp_path / "data")
    write_previews(tmp_path / "previews", [example])

    assert example.describe()["folder_files"] == ["evidence.txt", "mqpar.xml", "peptides.txt"]
    assert set(example.files()) == {"evidence.txt", "mqpar.xml", "peptides.txt"}
    assert (tmp_path / "previews" / example.id / "peptides.txt.json").is_file()


def test_example_previews_show_the_first_lines_of_each_file(
    tmp_path: Path, options: Options, examples: list[Example]
) -> None:
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    write_previews(store.previews, examples)
    (example,) = examples
    with TestClient(
        create_app(store, Worker(store, 60), outputs, software, 1 << 30, examples)
    ) as test_client:
        preview = test_client.get(f"/api/examples/{example.id}/head/report.tsv").json()
        sdrf = test_client.get(f"/api/examples/{example.id}/head/sdrf.tsv").json()
        absent = test_client.get(f"/api/examples/{example.id}/head/corpus.csv")

    assert preview["name"] == "report.tsv", "the name the page serves, not the stored one"
    assert preview["format"] == "text"
    assert preview["lines"][0].startswith("Run\tModified.Sequence\t")
    assert len(preview["lines"]) == HEAD_LINES and not preview["complete"]
    assert sdrf["complete"] and sdrf["lines"][0].startswith("raw_file\tcondition")
    assert absent.status_code == 404


def test_a_parquet_head_is_its_first_rows_as_tab_separated_text(tmp_path: Path) -> None:
    path = tmp_path / "report.parquet"
    pl.DataFrame({"Run": ["a", "b", "c"], "Intensity": [1.0, 2.0, None]}).write_parquet(path)

    preview = head(path)

    assert preview["format"] == "parquet"
    assert preview["lines"] == ["Run\tIntensity", "a\t1.0", "b\t2.0", "c\t"]
    assert preview["complete"]


def test_a_row_no_hint_describes_is_no_example(tmp_path: Path) -> None:
    root = tmp_path / "data"
    corpus = _corpus(root)
    diann = load_hints()[1]["DIA-NN"]
    unclaimed = SoftwareHint(exports=(), params_required=False, sources=diann.sources)

    offered = load_examples(corpus, root, None, ["DIA-NN"], {"DIA-NN": unclaimed}, rule_versions())

    assert offered == [], "the row stays in the corpus, but the page claims no export for it"
    assert len(_load(corpus, root)) == 1


def test_an_example_brings_its_secondary_file_as_one_folder(
    tmp_path: Path, options: Options
) -> None:
    root = tmp_path / "data"
    corpus = _corpus(root)
    folder = root / "submissions" / "one"
    (folder / "input_file_secondary.tsv").write_text("x\n", encoding="utf-8")
    (example,) = _load(corpus, root)
    software, outputs = options
    store = JobStore(tmp_path / "store", ttl_seconds=3600)
    write_examples(store.examples_path, [example])
    worker = Worker(store, timeout_seconds=600)
    with TestClient(
        create_app(store, worker, outputs, software, 1 << 30, [example])
    ) as test_client:
        listed = test_client.get("/api/examples.json").json()
        form = {"output": "apb2:hdf5", "example": example.id, "example_files": "data"}
        job = test_client.post("/api/jobs", data=form).json()["id"]
        worker.join()

    assert example.secondary == (folder / "input_file_secondary.tsv",)
    request = json.loads((store.jobs / job / "request.json").read_text(encoding="utf-8"))
    linked = store.jobs / job / request["data"]
    assert linked.is_dir(), "two result files reach the converter as one folder"
    assert sorted(path.name for path in linked.iterdir()) == [
        "input_file.tsv",
        "input_file_secondary.tsv",
    ], "DIA-NN writes one file, so the stored names stay"
    assert listed["examples"][0]["stored_data"] == ["input_file.tsv", "input_file_secondary.tsv"]


def test_several_uploaded_result_files_become_one_folder(
    client: tuple[TestClient, Worker, JobStore], tmp_path: Path
) -> None:
    test_client, worker, store = client
    files = [
        ("data", ("precursors.tsv", b"a\tb\n1\t2\n")),
        ("data", ("precursor.matrix.tsv", b"a\tb\n1\t2\n")),
    ]
    form = {"output": "apb2:hdf5", "software": "AlphaDIA"}
    job = test_client.post("/api/jobs", data=form, files=files).json()["id"]
    worker.join()

    request = json.loads((store.jobs / job / "request.json").read_text(encoding="utf-8"))
    assert request["data"] == "input/data/precursors"
    names = sorted(path.name for path in (store.jobs / job / request["data"]).iterdir())
    assert names == ["precursor.matrix.tsv", "precursors.tsv"]

    twice = [("data", ("lfq.tsv", b"x\n")), ("data", ("lfq.tsv", b"y\n"))]
    assert test_client.post("/api/jobs", data=form, files=twice).status_code == 400
    assert len(list(store.jobs.iterdir())) == 1, "a refused upload leaves no job behind"


def test_every_offered_output_is_described_with_https_links(options: Options) -> None:
    _, outputs = options
    described = about()

    assert {output.id for output in outputs} == set(described), "no stale or missing entries"
    for entry in described.values():
        assert entry["what"] and entry["open"] and entry["links"]
        assert entry["group"] in {"tool", "apb2"}
        links = entry["links"]
        assert isinstance(links, list)
        assert all(str(link["url"]).startswith("https://") for link in links)
