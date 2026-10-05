# APB Export

APB2 results exported in the AnnData and MuData flavours downstream tools read: one export rule per target, one subcommand per target. Each target fixes the APB2 levels it reads; `--abundance` picks which layer of that level becomes X, and `--annotation` attaches an SDRF or prolfquapp-style sample table, as `apb2 annotate` does.

```bash
apb-export msmu report.tsv result.h5mu --params report.log.txt
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
- `uns`: msmu's reader settings, plus `uns["apb"]` with the export rule, each field's APB source, the parse plan and versions

`peptide` is ProForma for every vendor (`AC[UNIMOD:4]K`), not each vendor's own notation. msmu's modification parser reads it, and `to_ptm` takes the tag as written, for example `"[UNIMOD:21]"`.

## Export rules

The layout above is declared in [msmu's export rule](src/apb_export/export_rules/documents/msmu/v0_4/rules.json), not coded. Rules follow APB2's rule schema read in reverse: an entry's `name` is the target's field, its `source` an APB field every vendor shares, such as a ProForma column, a role, the run key or a catalogued meaning. A catalogued source with `"level": "protein"` reads the protein level through each variable's protein group: a cell takes its group's value in its run. One engine checks every source before it builds anything; the published [JSON Schemas](src/apb_export/export_rules/documents/_schema) describe the documents.

| Target | File | APB2 levels | Layout |
| --- | --- | --- | --- |
| [msmu](src/apb_export/export_rules/documents/msmu/v0_4/rules.json) | `.h5mu` | ion | long `psm` modality, one feature per run and precursor |
| [prolfquapp](src/apb_export/export_rules/documents/prolfqua/v2_11/rules.json) | `.h5ad` | ion, reading protein | wide; `uns["prolfquapp"]` for `LFQData_from_anndata`, factors from `apb2 annotate` columns; q-value layers `qValue` (precursor, per run), `pg_qValue` (protein group, per run) and `pg_qValue_experiment` (protein group, library or experiment-wide), each where the vendor reports it |
| [ProteoPy](src/apb_export/export_rules/documents/proteopy/v0_1/rules.json) | `.h5ad` | protein | wide; `sample_id` and `protein_id` beside their indexes |
| [AlphaPeptTools](src/apb_export/export_rules/documents/alphapepttools/v0_4/rules.json) | `.h5mu` | ion, peptide, protein | one wide modality per level, linked by id columns for `mulink_from_anndatas` |

## Confidence, by meaning

- `q_value`: the library q-value when match-between-runs filled it, else the experiment-wide q-value, else the per-run q-value; the order msmu's DIA-NN reader uses
- `PEP`: the vendor's posterior error probability
- Lookup: apb-catalog's `identification_confidence` and `proteobench_entrapment` catalogues, never vendor column names
- Not reported: no `q_value` column, so msmu's filter cannot act on an invented value; `PEP` is NaN, because msmu's `to_peptide` requires the column
- Unreviewed rule: apb-catalog answers "unknown" and the command stops; unknown is not missing

## Protein groups

Written as msmu's readers write them: members split on `;`, the accession taken from a UniProt `db|ACC|NAME` entry, contaminants marked `contam_`, `Cont_` or `CON__` written `Cont_ACC` and flagged in `contaminant`.

## Limits

- Needs APB2's `ion` level; a vendor export without one (Sage's `lfq.tsv` is peptide-level) is refused
- Label-free only: a result keyed by more than one observation column, as multiplexed designs are, is refused
- No decoys: `decoy` is 0; msmu's target-decoy q-values need decoys, so use `calculate_q=False`
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
