# APB msmu — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --frozen --group dev` |
| Format | `.venv/bin/ruff format src tests && .venv/bin/ruff check --fix src tests` |
| Lint | `.venv/bin/ruff check src tests` |
| Typecheck | `.venv/bin/pyright` |
| Dependencies | `.venv/bin/deptry .` |
| Tests | `.venv/bin/pytest --cov --cov-branch` |
| Build | `uv build && .venv/bin/twine check dist/*` |
| Full gate | `make check` |

## Architecture

One step, vendor output to msmu's MuData. Modules compose downward only:

- `cli`: the one command; APB2's `ParseRuleCompiler` reads the vendor files, then `MsmuExporter` exports; refuses overwrite, reports unreviewed catalogue answers
- `api`: `MsmuExporter(abundance=...).export(parsed, /)`, the family's shape: configuration bound once, one verb on `ParsedLevels`
- `container`: AnnData/MuData assembly; the only module that imports them
- `psm`: Polars only; unpivots APB2's wide ion level into msmu's one-row-per-cell psm table
- `confidence`: q-value and PEP layers by meaning through apb-catalog

Rules for this package:

- **APB2 owns conversion and persistence.** Never define a `convert`, reader or writer here; use `ParseRuleCompiler`, `read_parsed_levels` and `write_parsed_levels`.

- **Meanings come from apb-catalog, never vendor column names.** A vendor without a catalogued field gets no column (q-value) or NaN (PEP); an unreviewed rule is an error, not a guess.
- **Match msmu's readers, not our preferences.** Layout, column names and order, q-value precedence and protein canonicalisation follow `msmu.read_diann`; the parity tests hold this.
- **No real vendor data** in tests, docs or examples; tests synthesise vendor files.
- **Keep this repository private.**

## Code conventions

- Fully annotate every function and method in `src/` and `tests/`, including
  private functions, callbacks, generators, fixtures, and special methods.
- Standard Pyright strict and Ruff are mandatory. Do not create baselines,
  blanket exclusions, file-wide ignores, or unqualified `# type: ignore`.
- Ruff is the sole formatter and linter. Do not add Black, isort, Flake8, mypy,
  or another overlapping formatter/type checker.
- Keep `__init__.py` empty and import from defining modules inside this package. Other anndata_bridge packages import this one only from `apb_msmu.api`, and it imports them only from theirs. The CLI imports this package only from `apb_msmu.api` too.
- Use Google-style docstrings for public APIs and the configured 100-character
  line length.

## Dependency rules

### MUST

- Declare every imported runtime dependency directly in `[project.dependencies]`.
- Put tests, linting, typing, building, and documentation tools in dependency
  groups; optional user-facing capabilities belong in extras.
- Update `pyproject.toml` and `uv.lock` together and run `make check`.

### SHOULD

- Prefer the standard library, then an existing direct dependency, then a small,
  maintained, typed dependency.
- Keep source independent of test, build, documentation, and CLI-only packages.

### MUST NOT

- Depend on unpinned branches or undeclared transitive dependencies.
- Add parallel manifests, lockfiles, formatters, type checkers, or test runners.
- Silence a dependency or typing defect instead of fixing its source.

## Workflow

1. Preserve unrelated worktree changes.
2. Add or update focused tests with each behavioral change.
3. Run the smallest relevant check while iterating.
4. Run `make check` before handoff and report its actual result.
