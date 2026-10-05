# APB Export

APB2 results exported in the formats downstream tools read; msmu is the first target, more follow through export rules.

One command from vendor output to the MuData [msmu](https://github.com/bertis-informatics/msmu) builds itself: APB2 reads the vendor files, apb-catalog names their confidence fields, and msmu takes the result with `msmu.read_h5mu` and runs unchanged.

```bash
apb-export report.tsv result.h5mu --params report.log.txt
```

```python
import msmu as mm

mdata = mm.read_h5mu("result.h5mu")
mdata = mm.pp.add_filter(mdata, modality="psm", column="q_value", keep="lt", value=0.01)
mdata = mm.pp.apply_filter(mdata, modality="psm", on="var")
mdata = mm.pp.to_peptide(mdata, calculate_q=False)
```

From Python, APB2 converts and `MsmuExporter` exports, the shape of the other APB tools:

```python
from pathlib import Path

from apb2.api import ParseRuleCompiler
from apb_export.api import MsmuExporter

parsed = ParseRuleCompiler(
    Path("report.tsv"), Path("report.log.txt"), requested_levels=("ion",)
).compile().parse()
MsmuExporter().export(parsed).write_h5mu("result.h5mu")
```

An APB2 result already on disk works the same way through `apb2.api.read_parsed_levels`.

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

The layout above is declared in [msmu's export rule](src/apb_export/export_rules/documents/msmu/v0_4/rules.json), not coded. Rules follow APB2's rule schema read in reverse: an entry's `name` is the target's field, its `source` an APB field every vendor shares, such as a ProForma column, a role, the run key or a catalogued meaning. One engine checks every source before it builds anything; the published [JSON Schemas](src/apb_export/export_rules/documents/_schema) describe the documents.

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

The development release expects APB2 and apb-catalog checked out beside this repository, at the revisions pinned in [quality.yml](.github/workflows/quality.yml):

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

`uv sync --group dev --group parity` adds msmu, and with it the tests that compare an export against `msmu.read_diann` through msmu's own pipeline. Tests use synthetic vendor files only.
