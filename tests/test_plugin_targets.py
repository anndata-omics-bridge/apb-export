"""Plugin targets: discovery, the generated command, and the rule features only they use so far.

``tests/toy_target`` stands in for an installed plugin: a synthetic multi-level ``.h5mu`` rule
with SDRF ``obs`` columns, catalogue qualifiers, its own ``uns`` and feature links.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from importlib.metadata import EntryPoint, EntryPoints
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Any

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from apb2.api import ParseRuleCompiler
from apb_fasta.api import FastaAnnotator

from apb_export.api import Exporter
from apb_export.cli import add_plugin_commands, app
from apb_export.container import link_features
from apb_export.engine import CompiledExport
from apb_export.export_rules import loader
from apb_export.export_rules.loader import PLUGIN_GROUP, effective_rules, plugin_folders
from apb_export.export_rules.schema import ExportRuleDocument
from apb_export.web import options
from conftest import RUNS, DiannInput, apb_part, export_sources

TOY = Path(__file__).parent / "toy_target"
# Each synthetic precursor's FASTA proteins: PEPTIDEK also occurs in P8, outside its group P1;
# SHAREDK occurs in P5 and P6 while its group names P5 and a CON__P6 the FASTA lacks;
# LASTPEPK occurs nowhere.
FASTA = (
    ">sp|P1|ONE_HUMAN One\nMPEPTIDEKR\n"
    ">sp|P2|P2_HUMAN Two\nMAACLLKR\n"
    ">sp|P3|THREE_YEAST Three\nMYEASTPEPKR\n"
    ">sp|Cont_P4|FOUR_BOVIN Four\nMCONTPEPKR\n"
    ">sp|P5|FIVE_HUMAN Five\nMSHAREDKR\n"
    ">sp|P6|SIX_HUMAN Six\nMSHAREDKW\n"
    ">sp|P8|EIGHT_HUMAN Eight\nMPEPTIDEKW\n"
)
GROUP_SPECIFIC = {
    "PEPTIDEK/2": False,
    "AAC[UNIMOD:4]LLK/2": True,
    "YEASTPEPK/2": True,
    "CONTPEPK/2": True,
    "SHAREDK/2": False,
    "LASTPEPK/2": None,
}
SDRF_COLUMNS = (
    "source name",
    "characteristics[organism]",
    "characteristics[biological replicate]",
    "assay name",
    "technology type",
    "comment[proteomics data acquisition method]",
    "comment[label]",
    "comment[data file]",
    "factor value[spiked compound]",
)


def _plugins(folders: dict[str, Traversable]) -> dict[str, Traversable]:
    return folders


@pytest.fixture
def toy(monkeypatch: pytest.MonkeyPatch) -> Iterator[str]:
    """The toy target, as if a plugin package had installed it, with its generated command."""
    monkeypatch.setattr(loader, "_plugin_folders", lambda: _plugins({"toy": TOY}))
    add_plugin_commands()
    yield "toy"
    del app["toy"]


def _run(*argv: str) -> object:
    return app(list(argv), result_action="return_value")


def _fasta(folder: Path) -> Path:
    path = folder / "proteins.fasta"
    path.write_text(FASTA, encoding="utf-8")
    return path


def _sdrf(folder: Path) -> Path:
    rows = [
        (
            f"sample_{run[4]}",
            "homo sapiens",
            run[5],
            f"assay_{run}",
            "proteomic profiling by mass spectrometry",
            "Data-independent acquisition",
            "label free sample",
            f"{run}.raw",
            run[4],
        )
        for run in RUNS
    ]
    path = folder / "samples.sdrf.tsv"
    lines = ["\t".join(SDRF_COLUMNS), *("\t".join(row) for row in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _export(diann: DiannInput, folder: Path, *extra: str) -> md.MuData:
    target = folder / "toy.h5mu"
    assert _run("toy", str(diann.report), str(target), "--params", str(diann.log), *extra) == 0
    return md.read_h5mu(target)


def _modality(container: md.MuData, name: str) -> ad.AnnData:
    modality = container.mod[name]
    assert isinstance(modality, ad.AnnData)
    return modality


def _var(modality: ad.AnnData) -> pd.DataFrame:
    assert isinstance(modality.var, pd.DataFrame)
    return modality.var


def _obs(modality: ad.AnnData) -> pd.DataFrame:
    assert isinstance(modality.obs, pd.DataFrame)
    return modality.obs


def _document() -> dict[str, Any]:
    text = (TOY / "v1" / "rules.json").read_text(encoding="utf-8")
    document: dict[str, Any] = json.loads(text)
    return document


def test_an_installed_target_is_listed_and_applied(toy: str) -> None:
    """Exporter treats a plugin target like a packaged one."""
    assert toy in Exporter.targets()
    exporter = Exporter(toy)
    assert exporter.levels == ("protein", "peptidoform", "ion")
    assert exporter.extension == ".h5mu"


def test_the_generated_command_takes_the_general_options(
    toy: str, capsys: pytest.CaptureFixture[str]
) -> None:
    """FASTA and SDRF, but no --abundance: each level's primary layer becomes X."""
    _run(toy, "--help")
    shown = capsys.readouterr().out
    assert "--fasta" in shown
    assert "--sdrf" in shown
    assert "--params" in shown
    assert "--abundance" not in shown


def test_a_plugin_named_like_a_packaged_target_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(loader, "_plugin_folders", lambda: _plugins({"msmu": TOY}))
    with pytest.raises(ValueError, match=r"plugin export targets \['msmu'\]"):
        Exporter.targets()


def test_two_plugins_naming_one_target_are_refused() -> None:
    points = EntryPoints(
        EntryPoint(name="toy", value=f"package_{index}", group=PLUGIN_GROUP) for index in range(2)
    )
    with pytest.raises(ValueError, match=r"several installed plugins name the export targets"):
        plugin_folders(points)


def test_the_web_page_offers_only_the_targets_it_describes(toy: str) -> None:
    """An installed plugin target has no description there, so the page leaves it out."""
    offered = {output.id for output in options.outputs()}
    assert f"apb-export:{toy}" not in offered
    assert "apb-export:msmu" in offered


def test_the_target_writes_the_modalities_the_result_holds(
    toy: str, diann: DiannInput, tmp_path: Path
) -> None:
    """Protein and ion modalities; the peptidoform level the result lacks is left out."""
    container = _export(diann, tmp_path, "--sdrf", str(_sdrf(tmp_path)))

    assert sorted(container.mod) == ["ion", "protein"]
    dataset = {key: container.uns[key] for key in container.uns if key != "apb"}
    assert dataset["toy_version"] == "1"
    assert list(dataset["toy_modalities"]) == ["protein", "ion"]
    assert dataset["search_engine"] == "DIA-NN"
    assert _modality(container, "ion").uns["level"] == "ion"
    record = apb_part(container.uns["apb"])["export"]["provenance"]
    assert record["sources"]["uns.search_engine"] == {"location": "rule", "name": "software_name"}


def test_the_sdrf_names_the_assays_and_fills_the_samples(
    toy: str, diann: DiannInput, tmp_path: Path
) -> None:
    """assay_id is the SDRF's assay name; the run stays as raw_file."""
    ion = _modality(_export(diann, tmp_path, "--sdrf", str(_sdrf(tmp_path))), "ion")
    obs = _obs(ion)

    assert list(obs.index) == [f"assay_{run}" for run in RUNS]
    assert list(obs["raw_file"]) == list(RUNS)
    assert list(obs["condition"]) == [run[4] for run in RUNS]
    assert set(obs["organism"]) == {"homo sapiens"}
    assert set(obs["label"]) == {"label free sample"}


def test_without_an_sdrf_the_run_is_the_assay(toy: str, diann: DiannInput, tmp_path: Path) -> None:
    """The obs index falls back to the run; SDRF-only columns are left out."""
    obs = _obs(_modality(_export(diann, tmp_path), "ion"))

    assert list(obs.index) == list(RUNS)
    assert "sample_id" not in obs.columns


def test_qualifiers_choose_the_catalogue_answer(
    toy: str, diann: DiannInput, tmp_path: Path
) -> None:
    """DIA-NN's q-value rates the identification; it reports none for the MS1 peak."""
    ion = _modality(_export(diann, tmp_path), "ion")

    assert "q_value" in ion.layers
    assert "ms1_q_value" not in ion.layers
    source = export_sources(ion)["layers.q_value"]
    assert (source["catalogue"], source["stage"]) == ("identification_confidence", "identification")


def test_apb_fasta_supplies_group_specificity(toy: str, diann: DiannInput, tmp_path: Path) -> None:
    """Specific when the group holds every protein the peptide matches; missing without a match."""
    ion = _modality(_export(diann, tmp_path, "--fasta", str(_fasta(tmp_path))), "ion")
    flags = _var(ion)["group_specific"].astype(object)
    found = {str(key): (None if value is pd.NA else bool(value)) for key, value in flags.items()}

    assert found == GROUP_SPECIFIC


def test_ions_link_to_the_protein_group_they_name(
    toy: str, diann: DiannInput, tmp_path: Path
) -> None:
    """feature_mapping is symmetric over every feature and joins each ion to its group."""
    container = _export(diann, tmp_path)
    protein, ion = _modality(container, "protein"), _modality(container, "ion")
    links = container.varp["feature_mapping"]

    assert links.shape == (protein.n_vars + ion.n_vars,) * 2
    assert (links != links.T).nnz == 0
    rows, columns = links[protein.n_vars :, : protein.n_vars].nonzero()
    pairs = {
        (ion.var_names[row], protein.var_names[column])
        for row, column in zip(rows, columns, strict=True)
    }
    assert pairs == {(name, group) for name, group in _var(ion)["protein_group"].items()}


def _compiled(document: dict[str, Any]) -> CompiledExport:
    return CompiledExport(effective_rules(ExportRuleDocument.model_validate(document), "test"))


def test_an_optional_level_without_its_index_is_left_out(diann: DiannInput, tmp_path: Path) -> None:
    """A protein level whose index the result cannot supply is skipped and recorded."""
    document = _document()
    document["tables"][0]["levels"]["protein"]["columns"]["var"][0]["source"] = "ProForma_fragment"
    compiled = _compiled(document)
    parsed = ParseRuleCompiler(diann.report, diann.log, compiled.levels).compile().parse()
    checked = FastaAnnotator.read((_fasta(tmp_path),)).verify_peptides(parsed).parsed

    container = compiled.export(checked)

    assert isinstance(container, md.MuData)
    assert list(container.mod) == ["ion"]
    assert apb_part(container.uns["apb"])["export"]["provenance"]["skipped_levels"] == ["protein"]
    assert "feature_mapping" not in container.varp


def test_a_required_layer_the_result_lacks_stops_the_export(diann: DiannInput) -> None:
    """DIA-NN reports no experiment-wide precursor q-value, which this variant requires."""
    document = _document()
    layers = document["tables"][0]["base"]["measurements"]["layers"]
    layers[1]["source"]["kind"] = "global_q_value"
    compiled = _compiled(document)
    parsed = ParseRuleCompiler(diann.report, diann.log, compiled.levels).compile().parse()

    with pytest.raises(ValueError, match="layer 'q_value': the result carries no layer"):
        compiled.export(parsed)


def test_links_skip_empty_groups_and_count_unknown_ones() -> None:
    """An empty group links nowhere; a group the target lacks is counted, not guessed."""
    proteins = ad.AnnData(np.zeros((1, 2)), var=pd.DataFrame(index=["P1", "P2"]))
    ions = ad.AnnData(
        np.zeros((1, 4)),
        var=pd.DataFrame({"protein_group": ["P1", "", None, "P9"]}, index=["a", "b", "c", "d"]),
    )
    container = md.MuData({"protein": proteins, "ion": ions})

    assert link_features(container, "feature_mapping", "protein_group", "protein") == (1, 1)
    assert container.varp["feature_mapping"].nnz == 2
