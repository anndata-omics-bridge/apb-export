"""Write the synthetic DIA-NN example and each target's export of it, for the consumer image.

Runs in apb-export's own environment, ``uv run python consumers/export_examples.py OUT``; the
consumer image then checks that each tool reads its file. No real vendor data is involved.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

from apb2.api import AnnotationCompiler, ParseRuleCompiler

from apb_export.api import Exporter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from conftest import RUNS, write_diann  # the tests' synthetic DIA-NN writer

# msmu's own reader quantifies DIA-NN by Precursor.Quantity; compare like with like.
ABUNDANCE = {"msmu": "Precursor_Quantity"}


def main(folder: Path) -> None:
    """Write ``report.tsv``, ``samples.tsv`` and one file per target into ``folder``."""
    folder.mkdir(parents=True, exist_ok=True)
    diann = write_diann(folder)
    annotation = folder / "samples.tsv"
    rows = "".join(f"{run}\t{run[4]}\n" for run in RUNS)
    annotation.write_text("raw_file\tcondition\n" + rows, encoding="utf-8")
    for target in Exporter.targets():
        exporter = Exporter(target, abundance=ABUNDANCE.get(target))
        parsed = (
            ParseRuleCompiler.detect(diann.report, diann.log, requested_levels=exporter.levels)
            .compile()
            .parse()
        )
        if target != "msmu":
            parsed = AnnotationCompiler().compile(annotation).parse(parsed).annotate().parsed
        path = folder / f"{target}{exporter.extension}"
        exporter.write(exporter.export(parsed), path)
        print(path)
        if target == "qpx":
            for level in exporter.levels:
                single = replace(parsed, levels={level: parsed.levels[level]})
                path = folder / f"qpx-{level}.qpx"
                exporter.write(exporter.export(single), path)
                print(path)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
