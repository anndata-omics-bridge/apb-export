# APB Export — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --group dev` |
| Format | `.venv/bin/ruff format src tests && .venv/bin/ruff check --fix src tests` |
| Lint | `.venv/bin/ruff check src tests` |
| Typecheck | `.venv/bin/pyright` |
| Dependencies | `.venv/bin/deptry .` |
| Tests | `.venv/bin/pytest --cov --cov-branch` |
| Docs | `make docs` |
| Build | `uv build && .venv/bin/twine check dist/*` |
| Full gate | `make check` |

## Architecture

Export rules declare what each target reads; one engine applies them. Packaged targets: msmu, prolfqua, ProteoPy and AlphaPeptTools; installed plugins add more through the `apb_export.targets` entry-point group. Modules compose downward only:

- `cli`: one subcommand per target; APB2's `ParseRuleCompiler` reads the target's levels, apb2's `AnnotationCompiler` applies `--annotation`, then `Exporter` exports; refuses overwrite and the wrong file type, reports unreviewed catalogue answers
- `api`: `Exporter(target, abundance=...).export(parsed)`, the family's shape: configuration bound once, one verb on `ParsedLevels`; `write(exported, path)` hands the result to APB2's `write_parsed_levels` with unprefixed var names; `targets()`, `levels` and `extension` describe the packaged rules
- `engine`: `CompiledExport`, one target's level rules checked once; picks the levels for the result's software, aligns several on one observation axis with APB2's `observation_groups()`, binds every entry's sources, and builds a `ParsedLevels`; a rule's `feature_links` become directed feature relations, fine to coarse, as mulink reads them
- `rows`: Polars only; aligns sources to feature, observation and `uns` rows
- `computed`: the conversions a computed entry's `how` names
- `sources`: APB references resolved against one level, catalogue lookups through apb-catalog, and catalogued protein-level fields read through each variable's protein group, named whole or by its leading protein
- `export_rules`: schema models, loader and the packaged `documents/<target>/<version>/rules.json`; `make schema` rewrites `documents/_schema`

Rules for this package:

- **APB2 owns conversion and persistence.** Never define a `convert`, reader or writer here; use `ParseRuleCompiler`, `read_parsed_levels` and `write_parsed_levels`. `Exporter.write` only chooses APB2's writer settings; keys stay columns, as APB2 writes them.

- **Export rules read apb2's rule schema in reverse.** Keep apb2's key names; an entry's `name` is the target field and its `source` an APB reference every vendor shares.
- **Meanings come from apb-catalog, never vendor column names.** A vendor without a catalogued field, or a field apb2 left out, gets no column (q-value) or NaN (PEP); an unreviewed rule is an error, not a guess. Library and global q-values are var sources only; do not add layer candidates for them. In a candidate list, a catalogue candidate the catalogue reviewed but cannot pick (ambiguous, or a requested qualifier unknown) passes to the next candidate and is recorded as `passed`; the last candidate still raises.
- **A software name only chooses a level.** A level block's `software` picks the level a target reads for that vendor, as prolfqua reads MaxQuant's peptides.txt; its fields still come from roles and the catalogue.
- **msmu: its names, our layout.** Column names, q-value precedence and protein canonicalisation follow `msmu.read_diann`; the layout is APB2's wide ion level, marked experimental (`uns["apb_export_status"]`), and msmu's developers are asked to read it. The consumer image's msmu check compares msmu's results with `read_diann`'s.
- **No real vendor data** in tests, docs or examples; tests synthesise vendor files.
- **The repository is public** since 8 October 2026; its docs site is built from `docs/` by the Pages workflow.

## Code conventions

- Fully annotate every function and method in `src/` and `tests/`, including
  private functions, callbacks, generators, fixtures, and special methods.
- Standard Pyright strict and Ruff are mandatory. Do not create baselines,
  blanket exclusions, file-wide ignores, or unqualified `# type: ignore`.
- Ruff is the sole formatter and linter. Do not add Black, isort, Flake8, mypy,
  or another overlapping formatter/type checker.
- Keep `__init__.py` empty and import from defining modules inside this package. Other anndata_bridge packages import this one only from `apb_export.api`, and it imports them only from theirs. The CLI imports this package only from `apb_export.api` too.
- Use Google-style docstrings for public APIs and the configured 100-character
  line length.

## Dependency rules

### MUST

- Declare every imported runtime dependency directly in `[project.dependencies]`.
- Put tests, linting, typing, building, and documentation tools in dependency
  groups; optional user-facing capabilities belong in extras.
- Update `pyproject.toml` and run `make check`; the repository commits no `uv.lock`.

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

The web application and its image are owned by the separate `apb-export-ui` project. This library exposes no web executable or optional web extra.
