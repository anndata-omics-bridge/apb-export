# Architecture

Export rules declare what each target reads; one engine applies them. APB2 owns conversion and canonical APB persistence. HDF5 exports use its writer; QPX and the `.csv` long table have target-specific serializers of the exported quantities.

## Modules

Modules compose downward only:

| Module | Role |
| --- | --- |
| `cli` | One subcommand per target: APB2 reads the target's levels, `--annotation` is applied, then `Exporter` exports |
| `api` | `Exporter(target, abundance=...).export(parsed)` and `.write(exported, path)` through APB2's writer; `targets()`, `levels` and `extension` describe the packaged rules |
| `engine` | `CompiledExport`: one target's rules checked once; picks the levels for the result's software, aligns several on one observation axis, and builds a `ParsedLevels` for APB2's `write_parsed_levels`, a rule's feature links as one of its directed feature relations |
| `rows` | Polars only; aligns sources to cell, feature, observation and `uns` rows |
| `computed` | The conversions a computed entry's `how` names |
| `sources` | APB references resolved against one level; catalogue lookups through apb-catalog; protein-level fields read through each protein group |
| `outputs` | Selects the serializer once: HDF5 through APB2, target-specific QPX Parquet, or the `.csv` long table; each format says whether it holds one level |
| `long_table` | One level as a CSV: a row per observation and feature some layer holds, obs columns, var columns and layers side by side, categories as their labels, missing as `NA` |
| `qpx` | QPX 1.1 nested types, sparse quantities, run/sample metadata and export manifest; RT conversion uses catalogued original units |
| `export_rules` | Schema models, loader, and the packaged `documents/<target>/<version>/rules.json` |

## Rules

- Export rules keep APB2's key names: an entry's `name` is the target field and its `source` an APB reference every vendor shares
- Meanings come from apb-catalog: a vendor without a catalogued field, or a field apb2 left out, gets no column (q-value) or NaN (PEP); an unreviewed rule is an error
- A software name only chooses a level: a level block's `software` picks the level a target reads, as prolfquapp reads MaxQuant's peptides.txt
- Names follow each tool's own readers: msmu's `read_diann` sets column names, q-value precedence and protein canonicalisation; the layout is APB2's wide level, for msmu marked experimental
- `make schema` rewrites the JSON Schemas in `documents/_schema` from the models
