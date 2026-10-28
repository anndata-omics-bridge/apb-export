"""QPX serialization of already exported APB quantities; no vendor parsing or aggregation."""

from __future__ import annotations

import json
import shutil
import tempfile
from dataclasses import replace
from importlib.metadata import version
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
from apb2.api import ParsedLevel, ParsedLevels

from apb_export.qpx.schema import arrow_schema

_MEMBERS = pl.col("proteins").str.split(";")


def _quantities(level: ParsedLevel) -> pl.DataFrame:
    """Join stored primary quantities to their feature and observation rows by position."""
    cells = level.layers[level.primary_layer_name].cells
    observed = cells.filter(pl.col("value").is_not_null() & ~pl.col("value").is_nan())
    if observed.select(pl.col("value").cast(pl.Float32).is_infinite().any()).item():
        raise ValueError("QPX float32 intensity cannot represent an infinite or overflowing value")
    return (
        observed.join(level.var.frame.with_row_index("var_pos"), on="var_pos", how="left")
        .join(level.obs.frame.with_row_index("obs_pos"), on="obs_pos", how="left")
        .sort("var_pos", "obs_pos")
    )


def _features(level: ParsedLevel) -> pl.DataFrame:
    frame = _quantities(level)
    if "proteins" not in frame.columns:
        frame = frame.with_columns(proteins=pl.lit(None, dtype=pl.String))
    members = _MEMBERS.list.eval(
        pl.struct(
            pl.element().alias("accession"),
            pl.lit(None, dtype=pl.Int32).alias("start"),
            pl.lit(None, dtype=pl.Int32).alias("end"),
            pl.lit(None, dtype=pl.String).alias("pre"),
            pl.lit(None, dtype=pl.String).alias("post"),
        )
    )
    return (
        _retention_time(frame, level)
        .with_row_index("feature_id")
        .with_columns(
            pl.col("feature_id").cast(pl.Int64),
            pg_accessions=members,
            anchor_protein=_MEMBERS.list.first(),
            intensities=pl.concat_list(
                pl.struct(pl.lit("LFQ").alias("label"), pl.col("value").alias("intensity"))
            ),
            scan=pl.lit([], dtype=pl.List(pl.Int32)),
        )
    )


def _retention_time(frame: pl.DataFrame, level: ParsedLevel) -> pl.DataFrame:
    """Convert only the original units resolved by the chromatography catalogue."""
    values: list[pl.Expr] = []
    for name, seconds in (("rt_second", 1), ("rt_minute", 60)):
        layer = level.layers.get(name)
        if layer is None:
            continue
        frame = frame.join(
            layer.cells.rename({"value": name}),
            on=["var_pos", "obs_pos"],
            how="left",
            validate="1:1",
        )
        values.append(pl.col(name) * seconds)
    rt = pl.coalesce(values) if values else pl.lit(None, dtype=pl.Float64)
    return frame.with_columns(rt.alias("rt"))


def _proteins(level: ParsedLevel) -> pl.DataFrame:
    return (
        _quantities(level)
        .with_row_index("pg_id")
        .with_columns(
            pl.col("pg_id").cast(pl.Int64),
            pg_accessions=_MEMBERS,
            anchor_protein=_MEMBERS.list.first(),
            grouped_runs=pl.concat_list("run_file_name"),
            label=pl.lit("LFQ"),
            intensity=pl.col("value"),
            peptides=pl.lit(
                [], dtype=pl.List(pl.Struct({"protein_name": pl.String, "peptide_count": pl.Int32}))
            ),
        )
    )


def _design(parsed: ParsedLevels) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Use run identities as technical sample identities; carry supplied annotations."""
    frames = [level.obs.frame for level in parsed.levels.values()]
    obs = pl.concat(frames, how="diagonal_relaxed")
    # A level may carry fewer annotation columns than another; retain known values.
    annotations = [name for name in obs.columns if name != "run_file_name"]
    for name in annotations:
        conflicts = obs.group_by("run_file_name").agg(pl.col(name).drop_nulls().n_unique())
        if conflicts.get_column(name).gt(1).any():
            raise ValueError(f"conflicting observation annotation {name!r} for the same QPX run")
    obs = obs.group_by("run_file_name", maintain_order=True).agg(
        *(pl.col(name).drop_nulls().first() for name in annotations)
    )
    samples = obs.rename({"run_file_name": "sample_accession"})
    # QPX 1.1 requires non-null organism/tissue strings. This is an explicit missing
    # marker, never an inferred organism or experimental grouping.
    samples = samples.with_columns(
        *(
            pl.col(name).fill_null("not available")
            if name in samples.columns
            else pl.lit("not available").alias(name)
            for name in ("organism", "organism_part")
        )
    )
    runs = obs.select("run_file_name").with_columns(
        run_accession=pl.col("run_file_name"),
        samples=pl.concat_list(
            pl.struct(
                pl.col("run_file_name").alias("sample_accession"),
                pl.lit("LFQ").alias("label"),
                pl.lit(None, dtype=pl.Int32).alias("biological_replicate"),
                pl.lit(None, dtype=pl.Int32).alias("technical_replicate"),
            )
        ),
    )
    return runs, samples


def _write(view: str, frame: pl.DataFrame, output: Path, record: dict[str, str]) -> None:
    schema = arrow_schema(view)
    columns = frame.select(
        *(
            pl.col(name) if name in frame.columns else pl.lit(None).alias(name)
            for name in schema.names
        )
    )
    table = columns.to_arrow().cast(schema)
    if view == "sample":
        for name in frame.columns:
            if name not in schema.names:
                table = table.append_column(name, frame.select(name).to_arrow().column(0))
    metadata = {
        **record,
        "file_type": f"{view}_file",
        "compression_format": "zstd",
    }
    table = table.replace_schema_metadata({k.encode(): v.encode() for k, v in metadata.items()})
    pq.write_table(table, output / f"apb.{view}.parquet", compression="zstd")


class QpxOutput:
    """Separate QPX views can represent every observation group without discarding levels."""

    @property
    def one_level(self) -> bool:
        """QPX writes a view for each level it reads."""
        return False

    def align(self, parsed: ParsedLevels, first_level: str) -> ParsedLevels:
        """Use APB2's safe run aliases, retaining levels from every observation group."""
        del first_level
        return replace(
            parsed,
            levels={
                name: level
                for group in parsed.observation_groups()
                for name, level in group.levels.items()
            },
        )

    def write(self, exported: ParsedLevels, output: Path) -> None:
        """Write a complete QPX directory, refusing existing destinations."""
        if output.suffix != ".qpx":
            raise ValueError("QPX export needs a new .qpx directory")
        if output.exists():
            raise FileExistsError(f"refusing to overwrite {output}")
        builders = {"feature": _features, "pg": _proteins}
        tables = {name: builders[name](level) for name, level in exported.levels.items()}
        if not any(frame.height for frame in tables.values()):
            raise ValueError("QPX export has no reported precursor or protein quantities")
        tables["run"], tables["sample"] = _design(exported)
        record = {
            "qpx_version": "1.1",
            "creator": "apb-export",
            "writer_version": version("apb-export"),
            "software_provider": "apb2",
        }
        output.parent.mkdir(parents=True, exist_ok=True)
        # Build beside the destination; a failed serializer never leaves a partial export.
        with tempfile.TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as temporary:
            staging = Path(temporary)
            for name, frame in tables.items():
                _write(name, frame, staging, record)
            manifest = {
                "format": "qpx",
                **record,
                "tables": {name: frame.height for name, frame in tables.items()},
                "export": exported.metadata["export"],
            }
            (staging / "export.json").write_text(json.dumps(manifest, indent=2) + "\n")
            # mkdir is exclusive even if another process created the destination meanwhile.
            output.mkdir()
            try:
                for path in staging.iterdir():
                    path.rename(output / path.name)
            except OSError:
                shutil.rmtree(output)
                raise
