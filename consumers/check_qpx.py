"""Validate synthetic QPX with qpx 1.1.6 and generate real pMultiQC 0.0.48 reports.

Run in a consumer environment containing qpx, pmultiqc and multiqc, not in the exporter.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
from qpx.core.data import FeatureSchema, PgSchema, RunSchema, SampleSchema

SCHEMAS = {"feature": FeatureSchema, "pg": PgSchema, "run": RunSchema, "sample": SampleSchema}


def main(folder: Path) -> None:
    """Require schema-valid tables and populated reports for all three level combinations."""
    for name in ("qpx", "qpx-peptidoform_ion", "qpx-protein"):
        source = (folder / f"{name}.qpx").resolve()
        for path in source.glob("*.parquet"):
            view = path.name.split(".")[-2]
            result = SCHEMAS[view].validate_full(pq.read_table(path), strict=True)
            assert result.is_valid, (path, result.errors)
        output = (folder / f"{name}-report").resolve()
        command = [
            sys.executable,
            "-m",
            "multiqc",
            "--qpx-plugin",
            str(source),
            "--outdir",
            str(output),
        ]
        completed = subprocess.run(command, cwd=folder, capture_output=True, text=True)
        log = completed.stdout + completed.stderr
        (folder / f"{name}-multiqc.log").write_text(log)
        assert completed.returncode == 0, log
        assert "ERROR" not in log and "failed" not in log.lower(), log
        assert (output / "multiqc_report.html").stat().st_size > 100_000
        data = json.loads((output / "multiqc_data/multiqc_data.json").read_text())
        plots = data["report_plot_data"]
        assert "identification_summary_table" in plots, plots.keys()
        if name != "qpx-protein":
            assert "peptide_intensity_distribution_box" in plots, plots.keys()
        if name != "qpx-peptidoform_ion":
            assert "protein_intensity_distribution_box" in plots, plots.keys()
        print(f"{name}: schema-valid, populated report; {len(data['report_plot_data'])} plots")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
