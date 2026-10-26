"""QPX quantities and identities, verified independently of the report renderer."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
import pytest
from apb2.api import AnnotationCompiler, ParsedLevels, ParseRuleCompiler, packaged_rule_declarations

from apb_export.api import Exporter
from conftest import RUNS, DiannInput, write_maxquant


@pytest.fixture
def parsed(diann: DiannInput) -> ParsedLevels:
    return (
        ParseRuleCompiler.detect(diann.report, diann.log, requested_levels=Exporter("qpx").levels)
        .compile()
        .parse()
    )


def test_quantities_membership_and_canonical_nested_types(
    parsed: ParsedLevels, diann: DiannInput, tmp_path: Path
) -> None:
    exporter = Exporter("qpx")
    output = tmp_path / "quantities.qpx"
    exporter.write(exporter.export(parsed), output)
    features = pl.read_parquet(output / "apb.feature.parquet")
    proteins = pl.read_parquet(output / "apb.pg.parquet")
    assert features.height == diann.cells
    assert proteins.height == diann.cells
    assert features.get_column("feature_id").n_unique() == features.height
    assert proteins.get_column("pg_id").n_unique() == proteins.height
    assert set(features.get_column("run_file_name")) == set(RUNS)
    assert features.filter(pl.col("peptidoform").str.contains("UNIMOD", literal=True)).height > 0
    assert proteins.filter(pl.col("pg_accessions").list.len() == 2).get_column("pg_accessions")[
        0
    ].to_list() == ["P5", "CON__P6"]
    value = features.filter(
        (pl.col("sequence") == "PEPTIDEK") & (pl.col("run_file_name") == "run_A1")
    )
    assert value.get_column("intensities")[0].to_list() == [{"label": "LFQ", "intensity": 900.0}]
    assert features.schema["scan"] == pl.List(pl.Int32)
    assert features.schema["observed_mz"] == pl.Float32
    assert features.get_column("observed_mz").null_count() == features.height
    assert features.get_column("rt").min() == 600.0, "catalogued DIA-NN minutes become seconds"
    assert features.get_column("modifications").null_count() == features.height
    metadata = pq.read_schema(output / "apb.feature.parquet").metadata
    assert metadata is not None and metadata[b"qpx_version"] == b"1.1"
    manifest = json.loads((output / "export.json").read_text())
    assert manifest["tables"]["feature"] == diann.cells
    assert manifest["export"]["provenance"]["export_rule"] == "qpx/v1_1/rules.json"


@pytest.mark.parametrize("level,view", [("peptidoform_ion", "feature"), ("protein", "pg")])
def test_either_level_is_sufficient_without_annotations(
    parsed: ParsedLevels, tmp_path: Path, level: str, view: str
) -> None:
    exporter = Exporter("qpx")
    output = tmp_path / f"{view}.qpx"
    exporter.write(exporter.export(replace(parsed, levels={level: parsed.levels[level]})), output)
    assert {path.name for path in output.glob("*.parquet")} == {
        f"apb.{view}.parquet",
        "apb.run.parquet",
        "apb.sample.parquet",
    }
    samples = pl.read_parquet(output / "apb.sample.parquet")
    assert samples.get_column("organism").to_list() == ["not available"] * len(RUNS)
    assert samples.get_column("sample_accession").to_list() == list(RUNS)


def test_zero_is_a_measurement_and_missing_cells_stay_absent(
    parsed: ParsedLevels, tmp_path: Path
) -> None:
    exporter = Exporter("qpx")
    exported = exporter.export(parsed)
    layer = exported.levels["feature"].layers["intensity"]
    count = layer.cells.height
    layer.cells = (
        layer.cells.with_row_index("row")
        .with_columns(value=pl.when(pl.col("row") == 0).then(0.0).otherwise(pl.col("value")))
        .drop("row")
    )
    output = tmp_path / "zero.qpx"
    exporter.write(exported, output)
    features = pl.read_parquet(output / "apb.feature.parquet")
    assert features.height == count
    assert features.get_column("intensities")[0].to_list() == [{"label": "LFQ", "intensity": 0.0}]


def test_full_group_membership_defines_identity(parsed: ParsedLevels, tmp_path: Path) -> None:
    exporter = Exporter("qpx")
    exported = exporter.export(parsed)
    level = exported.levels["pg"]
    level.var.frame = (
        level.var.frame.with_row_index("row")
        .with_columns(
            proteins=pl.when(pl.col("row") == 0)
            .then(pl.lit("P1;P2"))
            .when(pl.col("row") == 1)
            .then(pl.lit("P1;P3"))
            .otherwise(pl.col("proteins"))
        )
        .drop("row")
    )
    output = tmp_path / "groups.qpx"
    exporter.write(exported, output)
    groups = pl.read_parquet(output / "apb.pg.parquet").filter(pl.col("anchor_protein") == "P1")
    assert groups.get_column("pg_accessions").list.join(";").n_unique() == 2
    assert groups.get_column("pg_id").n_unique() == groups.height


@pytest.mark.parametrize("value", [float("inf"), 1e100])
def test_unrepresentable_intensities_fail_without_partial_output(
    parsed: ParsedLevels, tmp_path: Path, value: float
) -> None:
    exporter = Exporter("qpx")
    exported = exporter.export(parsed)
    layer = exported.levels["feature"].layers["intensity"]
    layer.cells = layer.cells.with_columns(value=pl.lit(value))
    output = tmp_path / "invalid.qpx"
    with pytest.raises(ValueError, match="float32 intensity"):
        exporter.write(exported, output)
    assert not output.exists()


def test_refuses_overwrite_and_wrong_suffix(parsed: ParsedLevels, tmp_path: Path) -> None:
    exporter = Exporter("qpx")
    exported = exporter.export(parsed)
    output = tmp_path / "existing.qpx"
    output.mkdir()
    (output / "keep").write_text("preserve")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        exporter.write(exported, output)
    assert (output / "keep").read_text() == "preserve"
    with pytest.raises(ValueError, match=r"\.qpx directory"):
        exporter.write(exported, tmp_path / "wrong.parquet")


def test_maxquant_uses_apb_run_aliases_without_aggregating(tmp_path: Path) -> None:
    inputs = write_maxquant(tmp_path)
    exporter = Exporter("qpx")
    parsed = (
        ParseRuleCompiler.detect(inputs.folder, inputs.params, exporter.levels).compile().parse()
    )
    output = tmp_path / "maxquant.qpx"
    exporter.write(exporter.export(parsed), output)
    features = pl.read_parquet(output / "apb.feature.parquet")
    proteins = pl.read_parquet(output / "apb.pg.parquet")
    assert set(features.get_column("run_file_name")) == {"raw_S1", "raw_S2", "raw_S3", "raw_S4"}
    assert set(proteins.get_column("grouped_runs").list.first()) == set(
        features.get_column("run_file_name")
    )
    assert proteins.height == 16


def test_rt_uses_catalogue_units_without_vendor_names_in_the_serializer(
    parsed: ParsedLevels, tmp_path: Path
) -> None:
    level = parsed.levels["peptidoform_ion"]
    # Same synthetic numeric RT values, held under a reviewed seconds-based rule.
    level.uns["provenance"] = {
        "rule_json": packaged_rule_declarations()[("alphadia/v2/rules.json", "peptidoform_ion")][0]
    }
    level.layers["RT_Observed"] = replace(level.layers.pop("RT"), layer_name="RT_Observed")
    exporter = Exporter("qpx")
    exported = exporter.export(replace(parsed, levels={"peptidoform_ion": level}))
    source = exported.levels["feature"].metadata["export"]
    assert isinstance(source, dict)
    provenance = source["provenance"]
    assert isinstance(provenance, dict)
    sources = provenance["sources"]
    assert isinstance(sources, dict)
    assert sources["layers.rt_second"] == {
        "location": "layers",
        "name": "RT_Observed",
        "catalog": "chromatography",
        "kind": "retention_time",
        "unit": "second",
    }
    output = tmp_path / "seconds.qpx"
    exporter.write(exported, output)
    assert pl.read_parquet(output / "apb.feature.parquet").get_column("rt").min() == 10.0


def test_absent_rt_is_not_invented(parsed: ParsedLevels, tmp_path: Path) -> None:
    del parsed.levels["peptidoform_ion"].layers["RT"]
    exporter = Exporter("qpx")
    output = tmp_path / "no-rt.qpx"
    exporter.write(exporter.export(parsed), output)
    features = pl.read_parquet(output / "apb.feature.parquet")
    assert features.get_column("rt").null_count() == features.height


def test_optional_sample_groups_are_preserved(parsed: ParsedLevels, tmp_path: Path) -> None:
    annotation = tmp_path / "samples.tsv"
    annotation.write_text("raw_file\tcondition\n" + "".join(f"{run}\t{run[4]}\n" for run in RUNS))
    annotated = AnnotationCompiler().compile(annotation).parse(parsed).annotate().parsed
    exporter = Exporter("qpx")
    output = tmp_path / "annotated.qpx"
    exporter.write(exporter.export(annotated), output)
    samples = pl.read_parquet(output / "apb.sample.parquet")
    assert dict(samples.select("sample_accession", "condition").iter_rows()) == {
        "run_A1": "A",
        "run_A2": "A",
        "run_B1": "B",
        "run_B2": "B",
    }


def test_precursor_quantities_do_not_require_a_protein_assignment(
    parsed: ParsedLevels, diann: DiannInput, tmp_path: Path
) -> None:
    level = parsed.levels["peptidoform_ion"]
    del level.var.roles["protein_assignment"]
    exporter = Exporter("qpx")
    output = tmp_path / "unassigned.qpx"
    exporter.write(exporter.export(replace(parsed, levels={"peptidoform_ion": level})), output)
    features = pl.read_parquet(output / "apb.feature.parquet")
    assert features.height == diann.cells
    assert features.get_column("anchor_protein").null_count() == diann.cells
    assert features.get_column("pg_accessions").null_count() == diann.cells


def test_composite_observations_are_refused(parsed: ParsedLevels) -> None:
    level = parsed.levels["peptidoform_ion"]
    level.obs = replace(
        level.obs,
        frame=level.obs.frame.with_columns(channel=pl.lit("126")),
        key_columns=("Run", "channel"),
    )
    with pytest.raises(ValueError, match="one observation per run"):
        Exporter("qpx").export(replace(parsed, levels={"peptidoform_ion": level}))


def test_separate_observation_groups_keep_both_quantified_levels(
    parsed: ParsedLevels, diann: DiannInput, tmp_path: Path
) -> None:
    protein = parsed.levels["protein"]
    protein.obs.frame = protein.obs.frame.with_columns(Run=pl.col("Run") + "_protein")
    exporter = Exporter("qpx")
    output = tmp_path / "separate.qpx"
    exporter.write(exporter.export(parsed), output)
    features = pl.read_parquet(output / "apb.feature.parquet")
    proteins = pl.read_parquet(output / "apb.pg.parquet")
    assert features.height == proteins.height == diann.cells
    assert set(features.get_column("run_file_name")) == set(RUNS)
    assert set(proteins.get_column("grouped_runs").list.first()) == {
        f"{run}_protein" for run in RUNS
    }
