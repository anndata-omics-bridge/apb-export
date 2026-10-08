# Architecture

Export rules declare what each target reads; one engine applies them. APB2 owns conversion and persistence, so apb-export defines no reader, writer or `convert`.

## Modules

Modules compose downward only:

| Module | Role |
| --- | --- |
| `cli` | One subcommand per target: APB2 reads the target's levels, `--annotation` is applied, then `Exporter` exports |
| `api` | `Exporter(target, abundance=...).export(parsed)`; `targets()`, `levels` and `extension` describe the packaged rules |
| `engine` | `CompiledExport`: one target's rules checked once; picks the levels for the result's software and builds one AnnData per level |
| `rows` | Polars only; aligns sources to cell, feature, observation and `uns` rows |
| `computed` | The conversions a computed entry's `how` names |
| `sources` | APB references resolved against one level; catalogue lookups through apb-catalog; protein-level fields read through each protein group |
| `container` | AnnData/MuData assembly; the only module that builds them |
| `web` | `apb-export-web`: FastAPI server and job runner; reaches the package only through `apb_export.api` |
| `export_rules` | Schema models, loader, and the packaged `documents/<target>/<version>/rules.json` |

## Rules

- Export rules keep APB2's key names: an entry's `name` is the target field and its `source` an APB reference every vendor shares
- Meanings come from apb-catalog: a vendor without a catalogued field gets no column (q-value) or NaN (PEP); an unreviewed rule is an error
- A software name only chooses a level: a level block's `software` picks the level a target reads, as prolfquapp reads MaxQuant's peptides.txt
- Layout follows each tool's own readers: msmu's `read_diann` sets column names, order, q-value precedence and protein canonicalisation
- `make schema` rewrites the JSON Schemas in `documents/_schema` from the models
