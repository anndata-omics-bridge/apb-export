"""The one command: vendor output in, msmu-ready .h5mu out."""

from __future__ import annotations

from pathlib import Path

import mudata as md
import pytest
from apb_catalog.resolver import UnresolvedField

from apb_msmu.api import MsmuExporter
from apb_msmu.cli import app, main
from conftest import DiannInput


def _run(*argv: str) -> int:
    return app(list(argv), result_action="return_value")


def test_writes_msmu_file(diann: DiannInput, tmp_path: Path) -> None:
    target = tmp_path / "nested" / "out.h5mu"

    assert _run(str(diann.report), str(target), "--params", str(diann.log)) == 0
    psm = md.read_h5mu(target)["psm"]
    assert psm.n_vars == diann.cells
    assert psm.uns["identification_file"] == str(diann.report), "the command names its source"


def test_refuses_overwrite_wrong_suffix_and_bad_input(diann: DiannInput, tmp_path: Path) -> None:
    existing = tmp_path / "out.h5mu"
    existing.write_text("keep", encoding="utf-8")

    assert _run(str(diann.report), str(existing), "--params", str(diann.log)) == 1
    assert existing.read_text(encoding="utf-8") == "keep"
    assert _run(str(diann.report), str(tmp_path / "out.h5ad"), "--params", str(diann.log)) == 1
    assert _run(str(tmp_path / "missing.tsv"), str(tmp_path / "x.h5mu")) == 1


def test_main_exits_with_the_status(
    diann: DiannInput, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("sys.argv", ["apb-msmu", str(diann.report), str(tmp_path / "x.h5ad")])
    with pytest.raises(SystemExit) as stopped:
        main()
    assert stopped.value.code == 1


def test_an_unreviewed_rule_is_reported_not_guessed(
    diann: DiannInput, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _unreviewed(_self: MsmuExporter, _parsed: object, /) -> md.MuData:
        error = UnresolvedField.__new__(UnresolvedField)
        LookupError.__init__(error, "confidence on ion (kind=q_value) is unknown")
        raise error

    monkeypatch.setattr(MsmuExporter, "export", _unreviewed)
    target = tmp_path / "out.h5mu"

    assert _run(str(diann.report), str(target), "--params", str(diann.log)) == 1
    assert not target.exists()
