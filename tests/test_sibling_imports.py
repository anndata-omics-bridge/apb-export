"""Other anndata_bridge packages are imported only through their api module."""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).parents[1]
SIBLINGS = frozenset(
    {
        "apb2",
        "apb_aggregate",
        "apb_catalog",
        "apb_fasta",
        "apb_msmu",
        "apb_proteobench",
        "protein_fasta",
        "prozor",
    }
) - {"apb_msmu"}


def _sibling_modules(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None and node.level == 0:
            modules.extend(
                f"{node.module}.{alias.name}" if node.module in SIBLINGS else node.module
                for alias in node.names
            )
    return [module for module in modules if module.split(".")[0] in SIBLINGS]


def test_siblings_are_imported_only_through_their_api() -> None:
    sources = sorted((ROOT / "src").rglob("*.py")) + sorted((ROOT / "tests").rglob("*.py"))
    offenders = [
        f"{path.relative_to(ROOT)}: {module}"
        for path in sources
        for module in _sibling_modules(path)
        if module != f"{module.split('.')[0]}.api"
    ]
    assert offenders == []
