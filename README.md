# APB Export

**[Online documentation](https://anndata-omics-bridge.github.io/apb-export/)** or its [source index](docs/index.md).

APB2 results exported in the AnnData and MuData flavours downstream tools read: one export rule per target, one subcommand per target. Each file gets APB2's `<file>.apb.json` sidecar, which APB Studio's viewer reads. Each target fixes the APB2 levels it reads; `--abundance` picks which layer of that level becomes X, and `--annotation` attaches an SDRF or prolfquapp-style sample table, as `apb2 annotate` does.

```bash
apb-export msmu report.tsv result.h5mu --params report.log.txt
apb-export msmu report.tsv result.h5mu --params report.log.txt --fasta proteins.fasta
apb-export prolfqua report.tsv result.h5ad --params report.log.txt --annotation dataset.csv
apb-export proteopy report.tsv proteins.h5ad --params report.log.txt
apb-export alphapepttools report.tsv linked.h5mu --params report.log.txt
```

The msmu file is the MuData [msmu](https://github.com/bertis-informatics/msmu) builds itself: apb-catalog names the confidence fields, and msmu takes the result with `msmu.read_h5mu` and runs unchanged.

```python
import msmu as mm

mdata = mm.read_h5mu("result.h5mu")
mdata = mm.pp.add_filter(mdata, modality="psm", column="q_value", keep="lt", value=0.01)
mdata = mm.pp.apply_filter(mdata, modality="psm", on="var")
mdata = mm.pp.to_peptide(mdata, calculate_q=False)
```

From Python, APB2 converts and `Exporter` exports, the shape of the other APB tools:

```python
from pathlib import Path

from apb2.api import ParseRuleCompiler
from apb_export.api import Exporter

exporter = Exporter("msmu")
parsed = ParseRuleCompiler(
    Path("report.tsv"), Path("report.log.txt"), requested_levels=exporter.levels
).compile().parse()
exporter.export(parsed).write_h5mu("result.h5mu")
```

An APB2 result already on disk works the same way through `apb2.api.read_parsed_levels`. `Exporter.targets()` lists the packaged targets; `extension` names the file each opens.

## What it writes

msmu's own layout, as `msmu.read_diann` writes it:

- one `psm` modality; runs are observations
- one feature per observed (run, precursor) cell, id `<run>.<ProForma precursor>`
- `X`: block-diagonal sparse float32, the chosen abundance layer (`--abundance`, APB2's primary layer by default)
- `var`: `proteins`, `peptide`, `stripped_peptide`, `filename`, `charge`, `peptide_length`, `decoy`, `contaminant`, `PEP`, and `q_value` where reported
- `varm["search_result"]`: every APB2 feature column and every layer value at that cell
- `uns`: msmu's reader settings, plus APB metadata as in an APB2 result: the MuData's `uns["apb"]` holds the result's root records and the `export` record (rule, target, versions); each modality's `uns["apb"]` holds that level's apb2 records unchanged and the `export` record with each field's APB source and a summary

`peptide` is ProForma for every vendor (`AC[UNIMOD:4]K`), not each vendor's own notation. msmu's modification parser reads it, and `to_ptm` takes the tag as written, for example `"[UNIMOD:21]"`.

## Export rules

The layout above is declared in [msmu's export rule](src/apb_export/export_rules/documents/msmu/v0_4/rules.json), not coded. Rules follow APB2's rule schema read in reverse: an entry's `name` is the target's field, its `source` an APB field every vendor shares, such as a ProForma column, a role, the run key or a catalogued meaning. A catalogued source with `"level": "protein"` reads the protein level through each variable's protein group, named whole or by the protein that leads it: a cell takes its group's value in its run. A level block may name the `software` it is for; an `.h5ad` rule lists its levels in order and writes the first one the result's software and levels allow. One engine checks every source before it builds anything; the published [JSON Schemas](src/apb_export/export_rules/documents/_schema) describe the documents.

| Target | File | APB2 levels | Layout |
| --- | --- | --- | --- |
| [msmu](src/apb_export/export_rules/documents/msmu/v0_4/rules.json) | `.h5mu` | ion | long `psm` modality, one feature per run and precursor |
| [prolfquapp](src/apb_export/export_rules/documents/prolfqua/v2_11/rules.json) | `.h5ad` | ion, or peptide for MaxQuant; reading protein | wide; `uns["prolfquapp"]` for `LFQData_from_anndata`, factors from `apb2 annotate` columns; confidence layers `qValue` (precursor, per run), `pg_qValue` (protein group, per run), `pg_qValue_experiment` (protein group, library or experiment-wide) and `pep`, each where the vendor reports it; a feature without a protein is left out and counted in the export record; with `--fasta`, the leading protein's `fasta.id`, `description`, `IDcolumn` (accession), `gene_name`, `protein_length` and `nr_tryptic_peptides`, as prolfquapp's own FASTA annotation writes them; without one, `fasta.id`, `description` and `IDcolumn` repeat `protein_Id` |
| [ProteoPy](src/apb_export/export_rules/documents/proteopy/v0_1/rules.json) | `.h5ad` | protein | wide; `sample_id` and `protein_id` beside their indexes; with `--fasta`, apb-fasta annotates the protein groups and adds `in_fasta`, `gene_id`, `description` and `fasta_organisms`, which ProteoPy does not read |
| [AlphaPeptTools](src/apb_export/export_rules/documents/alphapepttools/v0_4/rules.json) | `.h5mu` | ion, peptide, protein | one wide modality per level, linked by id columns for `mulink_from_anndatas`; ion and peptide modalities carry `contaminant` and, with `--fasta`, `in_fasta`, `fasta_proteins` and `fasta_organisms`, which AlphaPeptTools does not read |

## Confidence, by meaning

- `q_value`: the library q-value when match-between-runs filled it, else the experiment-wide q-value, else the per-run q-value; the order msmu's DIA-NN reader uses
- `PEP`: the vendor's posterior error probability
- Lookup: apb-catalog's `identification_confidence` and `proteobench_entrapment` catalogues, never vendor column names
- Not reported: no `q_value` column, so msmu's filter cannot act on an invented value; `PEP` is NaN, because msmu's `to_peptide` requires the column
- Unreviewed rule: apb-catalog answers "unknown" and the command stops; unknown is not missing

## Protein groups

Written as msmu's readers write them: members split on `;`, the accession taken from a UniProt `db|ACC|NAME` entry, and members spelled `contam_`, `Cont_` or `CON__` respelled `Cont_ACC`. `contaminant` flags a row that apb2 marked `apb_Contaminant`, a contaminant the software flags or adds itself, or that apb-fasta matched to a FASTA contaminant. `apb-export msmu --fasta` runs that FASTA check after conversion; without a FASTA only the vendor's own flags are known.

## Web app

`apb-export-web` serves a page that takes a vendor result, an optional parameter file, an optional SDRF or prolfquapp sample table and an optional FASTA, converts them to any target above or to an APB2 result (h5ad/h5mu, Parquet, DuckDB), and shows QC before download: intensity density per sample, CV per feature within each annotation group (across all samples without one), and detected and missing counts.

```bash
uv pip install 'apb-export[web]'
apb-export-web --store ./web-store --port 8770
```

The page asks for the software first, then shows which files it writes (result and parameter file names from ProteoBench's module documentation, packaged as `src/apb_export/web/hints.json`) and whether the parameter file is required. `--examples corpus.csv --examples-root DIR` offers one row per software of an APB Studio corpus CSV, such as `apb_studio/corpuses/routine.csv`, three ways: result file only, with its parameter file, and with parameters and SDRF. The SDRF is an `sdrf.tsv` beside the input, else `<module>.sdrf.tsv` from `--examples-sdrf DIR`, such as apb-proteobench's `src/apb_proteobench/data/modules`. The server links the chosen files into the job instead of uploading them; the data stays where the corpus keeps it, and this repository holds none. A row whose software's hints describe no export for it, such as AlphaDIA 1.10's ProteoBench-joined table, is not offered. At startup the server writes the first 20 lines of every example file, including each file of a folder result (a Parquet file's first rows as text), to the store, and the page shows them when a file name is clicked.

Each upload becomes a job folder; one worker runs the `apb-export` or `apb2` command in its own process and writes `status.json`, `qc.json` and the downloads, which the server only serves. Finished jobs are deleted after `--ttl-hours`. The page is a Lit/Vite project in [web](web); `make web` rebuilds its bundle into `src/apb_export/web/static`, which the wheel ships, so no Node runs in production.

To deploy, `make image` builds `ghcr.io/anndata-omics-bridge/apb-export:local` from this checkout and its apb2 and apb-catalog siblings ([Dockerfile](Dockerfile), build context `..`), and the `publish` workflow pushes `:vX.Y.Z` for each version tag. `make web-examples` copies the routine corpus, the dataset folders its rows read, the extra rows of [scripts/web_examples.csv](scripts/web_examples.csv) (ProteoBench submissions for the software routine lacks, apb2's committed Custom format and MetaMorpheus samples, and SDRFs from [scripts/web_sdrf](scripts/web_sdrf) for tools that name samples otherwise than raw files) and the module SDRFs into `build/web-examples`, one folder with one `examples.csv` that a server mounts read-only. The public instance runs behind the FGCZ portal; its descriptor is `portal/apb-export` in the [web-apps repository](https://gitlab.bfabric.org/proteomics/web-apps).

## Limits

- A result without the level a target reads is refused: msmu needs an `ion` level (Sage's `lfq.tsv` is peptide-level), ProteoPy a `protein` level
- Label-free only: a result keyed by more than one observation column, as multiplexed designs are, is refused
- `decoy` is apb2's `apb_Decoy`; most exports carry no decoys, and msmu's target-decoy q-values need them, so use `calculate_q=False` unless the vendor kept them
- Group separators other than `;` pass through unchanged; no APB2 rule declares its separator
- Retention time is not mapped to `rt`, since no catalogue names it; it stays in `varm["search_result"]`

## Development

The development release expects APB2 and apb-catalog checked out beside this repository, on `main`, as CI checks them out:

```text
anndata_bridge/
├── apb2/
├── apb-catalog/
└── apb-export/
```

```bash
uv sync --group dev
make check
```

Tests use synthetic vendor files only and install none of the target tools. Whether the tools read the files is checked in [consumers/](consumers): the [Consumers workflow](.github/workflows/consumers.yml) exports the synthetic examples, builds one image with msmu, ProteoPy and AlphaPeptTools, runs prolfquapp's check in prolfquapp's own image, and runs each tool's check, msmu's comparison with its own `read_diann` among them.
