"""The one command: vendor output in, msmu-ready .h5mu out."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import mudata as md
import pandas as pd
import pytest
from apb2.api import sidecar_path
from apb_catalog.api import UnresolvedField

from apb_export.api import Exporter
from apb_export.cli import app, main
from conftest import RUNS, DiannInput, write_maxquant


def _run(*argv: str) -> int:
    return app(list(argv), result_action="return_value")


def test_writes_msmu_file(diann: DiannInput, tmp_path: Path) -> None:
    target = tmp_path / "nested" / "out.h5mu"

    assert _run("msmu", str(diann.report), str(target), "--params", str(diann.log)) == 0
    psm = md.read_h5mu(target)["psm"]
    assert psm.n_vars == diann.cells
    assert psm.uns["identification_file"] == str(diann.report), "the command names its source"


def test_msmu_flags_the_fasta_contaminants_when_given_a_fasta(
    diann: DiannInput, tmp_path: Path
) -> None:
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(
        ">sp|P1|ONE_HUMAN One OS=Homo sapiens OX=9606 GN=ONE\nMPEPTIDEKAACLLKR\n"
        ">sp|Cont_P4|FOUR_BOVIN Four OS=Bos taurus OX=9913 GN=FOUR\nMCONTPEPKR\n"
    )
    plain, checked = tmp_path / "plain.h5mu", tmp_path / "checked.h5mu"

    assert _run("msmu", str(diann.report), str(plain), "--params", str(diann.log)) == 0
    assert (
        _run(
            "msmu",
            str(diann.report),
            str(checked),
            "--params",
            str(diann.log),
            "--fasta",
            str(fasta),
        )
        == 0
    )

    without = md.read_h5mu(plain)["psm"].var
    flagged = md.read_h5mu(checked)["psm"].var
    assert isinstance(without, pd.DataFrame) and isinstance(flagged, pd.DataFrame)
    assert without["contaminant"].sum() == 0, "DIA-NN marks no contaminant itself"
    assert flagged.loc["run_A1.CONTPEPK/2", "contaminant"] == 1
    assert flagged.loc["run_A1.PEPTIDEK/2", "contaminant"] == 0


def test_alphapepttools_carries_the_fasta_check_beside_its_columns(
    diann: DiannInput, tmp_path: Path
) -> None:
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(
        ">sp|P1|ONE_HUMAN One OS=Homo sapiens OX=9606 GN=ONE\nMPEPTIDEKAACLLKR\n"
        ">sp|Cont_P4|FOUR_BOVIN Four OS=Bos taurus OX=9913 GN=FOUR\nMCONTPEPKR\n"
    )
    plain, checked = tmp_path / "plain.h5mu", tmp_path / "checked.h5mu"
    common = ("--params", str(diann.log))

    assert _run("alphapepttools", str(diann.report), str(plain), *common) == 0
    assert (
        _run("alphapepttools", str(diann.report), str(checked), *common, "--fasta", str(fasta)) == 0
    )

    without = md.read_h5mu(plain)["precursors"].var
    flagged = md.read_h5mu(checked)["precursors"].var
    assert isinstance(without, pd.DataFrame) and isinstance(flagged, pd.DataFrame)
    assert "contaminant" in without.columns and "in_fasta" not in without.columns
    row = flagged.loc["CONTPEPK/2"]
    assert (bool(row["contaminant"]), bool(row["in_fasta"]), row["fasta_proteins"]) == (
        True,
        True,
        "sp|Cont_P4|FOUR_BOVIN",
    )


def test_refuses_overwrite_wrong_suffix_and_bad_input(diann: DiannInput, tmp_path: Path) -> None:
    existing = tmp_path / "out.h5mu"
    existing.write_text("keep", encoding="utf-8")

    assert _run("msmu", str(diann.report), str(existing), "--params", str(diann.log)) == 1
    assert existing.read_text(encoding="utf-8") == "keep"
    assert (
        _run("msmu", str(diann.report), str(tmp_path / "out.h5ad"), "--params", str(diann.log)) == 1
    )
    assert _run("msmu", str(tmp_path / "missing.tsv"), str(tmp_path / "x.h5mu")) == 1


def test_main_exits_with_the_status(
    diann: DiannInput, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "sys.argv", ["apb-export", "msmu", str(diann.report), str(tmp_path / "x.h5ad")]
    )
    with pytest.raises(SystemExit) as stopped:
        main()
    assert stopped.value.code == 1


def test_an_unreviewed_rule_is_reported_not_guessed(
    diann: DiannInput, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _unreviewed(_self: Exporter, _parsed: object, /) -> md.MuData:
        error = UnresolvedField.__new__(UnresolvedField)
        LookupError.__init__(error, "confidence on ion (kind=q_value) is unknown")
        raise error

    monkeypatch.setattr(Exporter, "export", _unreviewed)
    target = tmp_path / "out.h5mu"

    assert _run("msmu", str(diann.report), str(target), "--params", str(diann.log)) == 1
    assert not target.exists()


def _annotation(folder: Path) -> Path:
    """A prolfquapp-style sample table: the run key, then one condition per run."""
    table = folder / "samples.tsv"
    rows = "".join(f"{run}\t{run[4]}\n" for run in RUNS)
    table.write_text("raw_file\tcondition\n" + rows, encoding="utf-8")
    return table


def test_prolfqua_takes_its_factors_from_the_annotation(diann: DiannInput, tmp_path: Path) -> None:
    target = tmp_path / "out.h5ad"
    annotation = str(_annotation(tmp_path))

    assert (
        _run(
            "prolfqua",
            str(diann.report),
            str(target),
            "--params",
            str(diann.log),
            "--annotation",
            annotation,
        )
        == 0
    )
    adata = ad.read_h5ad(target)
    factors = adata.uns["prolfquapp"]["analysis_configuration"]["factors"]
    assert dict(factors) == {"condition": "condition"}
    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]


def test_prolfqua_annotates_proteins_from_the_fasta(diann: DiannInput, tmp_path: Path) -> None:
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(
        ">sp|P1|ONE_HUMAN One OS=Homo sapiens OX=9606 GN=ONE PE=1 SV=1\nMPEPTIDEKAACLLKR\n"
    )
    plain, checked = tmp_path / "plain.h5ad", tmp_path / "checked.h5ad"
    common = ("--params", str(diann.log))

    assert _run("prolfqua", str(diann.report), str(plain), *common) == 0
    assert _run("prolfqua", str(diann.report), str(checked), *common, "--fasta", str(fasta)) == 0

    without = ad.read_h5ad(plain)
    annotated = ad.read_h5ad(checked)
    annotation = annotated.uns["prolfquapp"]["protein_annotation"]
    assert (annotation["description"], annotation["cleaned_ids"]) == ("description", "IDcolumn")
    assert isinstance(without.var, pd.DataFrame) and isinstance(annotated.var, pd.DataFrame)
    assert (without.var["description"] == without.var["protein_Id"]).all()
    assert "gene_name" not in without.var.columns
    one = annotated.var[annotated.var["protein_Id"] == "P1"].iloc[0]
    assert (one["IDcolumn"], one["fasta.id"], one["gene_name"]) == ("P1", "sp|P1|ONE_HUMAN", "ONE")
    assert one["description"] == "One OS=Homo sapiens OX=9606 GN=ONE PE=1 SV=1"
    assert (one["protein_length"], one["nr_tryptic_peptides"]) == (16, 1)
    absent = annotated.var[annotated.var["protein_Id"] == "P7"].iloc[0]
    assert pd.isna(absent["description"]), "a protein the FASTA lacks gets no description"


def test_proteopy_and_alphapepttools_write_their_files_and_sidecars(
    diann: DiannInput, tmp_path: Path
) -> None:
    proteins, linked = tmp_path / "proteins.h5ad", tmp_path / "linked.h5mu"

    assert _run("proteopy", str(diann.report), str(proteins), "--params", str(diann.log)) == 0
    assert _run("alphapepttools", str(diann.report), str(linked), "--params", str(diann.log)) == 0
    assert "protein_id" in ad.read_h5ad(proteins).var.columns
    assert set(md.read_h5mu(linked).mod) == {"precursors", "proteins"}
    for path, levels in ((proteins, ["proteins"]), (linked, ["precursors", "proteins"])):
        sidecar = json.loads(sidecar_path(path).read_text(encoding="utf-8"))
        assert [level["name"] for level in sidecar["levels"]] == levels, "Studio's viewer reads it"


def test_a_target_refuses_the_other_file_type(diann: DiannInput, tmp_path: Path) -> None:
    wrong = str(tmp_path / "out.h5mu")

    assert _run("prolfqua", str(diann.report), wrong, "--params", str(diann.log)) == 1


def test_software_alone_exports_without_a_parameter_file(tmp_path: Path) -> None:
    (tmp_path / "txt").mkdir()
    inputs = write_maxquant(tmp_path / "txt")
    inputs.params.unlink()
    target = tmp_path / "out.h5ad"

    assert _run("prolfqua", str(inputs.folder), str(target), "--software", "MaxQuant") == 0
    assert ad.read_h5ad(target).n_obs == 4
