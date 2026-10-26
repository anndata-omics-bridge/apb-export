"""Run the QPX consumer gate with isolated, pinned downstream tools."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def check(folder: Path) -> None:
    """Validate the examples with QPX and the released, unchanged pMultiQC reader."""
    subprocess.run(
        [
            "uv",
            "run",
            "--no-project",
            "--python",
            "3.13",
            "--with",
            "qpx==1.1.6",
            "--with",
            "pmultiqc==0.0.48",
            "--with",
            "multiqc==1.35",
            "python",
            str(ROOT / "consumers/check_qpx.py"),
            str(folder.resolve()),
        ],
        check=True,
        cwd=ROOT,
    )


def main() -> None:
    """Use CI's existing examples or generate synthetic examples in a temporary directory."""
    if len(sys.argv) == 2:
        check(Path(sys.argv[1]))
        return
    with tempfile.TemporaryDirectory(prefix="apb-qpx-consumer-") as temporary:
        folder = Path(temporary)
        subprocess.run(
            [sys.executable, str(ROOT / "consumers/export_examples.py"), str(folder)],
            check=True,
            cwd=ROOT,
        )
        check(folder)


if __name__ == "__main__":
    main()
