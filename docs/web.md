# Web app

`apb-export-web` serves a page that converts one vendor result and shows QC before download.

## Inputs and outputs

- Inputs: a vendor result, an optional parameter file, an optional SDRF or prolfquapp sample table, an optional FASTA
- Outputs: any target on [Targets and layouts](targets.md), or an APB2 result as h5ad/h5mu, Parquet or DuckDB
- QC: intensity density per sample, CV per feature within each annotation group, detected and missing counts

```bash
uv pip install 'apb-export[web]'
apb-export-web --store ./web-store --port 8770
```

The page asks for the software first, then shows the files it expects: result and parameter file names from ProteoBench's module documentation, packaged as `src/apb_export/web/hints.json`, and whether the parameter file is required.

## Examples

`--examples corpus.csv --examples-root DIR` offers one row per software of an APB Studio corpus CSV, three ways:

- Result file only
- With its parameter file
- With parameters and SDRF

The SDRF is an `sdrf.tsv` beside the input, else `<module>.sdrf.tsv` from `--examples-sdrf DIR`. The server links the chosen files into the job instead of uploading them; the data stays where the corpus keeps it. At startup it writes the first 20 lines of every example file to the store, and the page shows them when a file name is clicked.

## Jobs

- Each upload becomes a job folder
- One worker runs the `apb-export` or `apb2` command in its own process
- The worker writes `status.json`, `qc.json` and the downloads; the server only serves them
- Finished jobs are deleted after `--ttl-hours`

The page is a Lit/Vite project in `web/`; `make web` rebuilds its bundle into `src/apb_export/web/static`, which the wheel ships, so no Node runs in production.

## Deploy

- `make image` builds `ghcr.io/anndata-omics-bridge/apb-export:local` from this checkout and its sibling packages, with build context `..`
- The `publish` workflow pushes `:vX.Y.Z` for each version tag
- `make web-examples` collects the example corpus, its dataset folders and module SDRFs into `build/web-examples`, one folder a server mounts read-only
