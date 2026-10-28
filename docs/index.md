# APB Export

[APB2](https://anndata-omics-bridge.github.io/apb2/) results exported as the AnnData and MuData files downstream tools open: one export rule per target, one subcommand per target.

APB2 converts the vendor output; apb-export applies the target's packaged export rule and writes the format that tool reads. HDF5 files get APB2's `<file>.apb.json` sidecar; QPX Parquet tables have an export manifest. Field meanings come from [apb-catalog](https://anndata-omics-bridge.github.io/apb-catalog/), never from vendor column names.

## Targets

| Target | File | APB2 levels | Opened with |
| --- | --- | --- | --- |
| [msmu](https://github.com/bertis-informatics/msmu) | `.h5mu` | peptidoform_ion | `msmu.read_h5mu` |
| prolfqua | `.h5ad` | peptidoform_ion, or peptide for MaxQuant | prolfquapp's `LFQData_from_anndata` |
| ProteoPy | `.h5ad` | protein | ProteoPy |
| AlphaPeptTools | `.h5mu` | peptidoform_ion, peptide, protein | `mulink_from_anndatas` |
| QFeatures | `.h5mu` | peptidoform_ion, peptide, protein | `MuData::readH5MU`, then `QFeatures` |
| QPX | `.qpx` directory | peptidoform_ion and/or protein | unchanged pMultiQC `--qpx-plugin` |

[Targets and layouts](targets.md) describes what each file holds.

## Quick start

```bash
apb-export msmu report.tsv result.h5mu --params report.log.txt
apb-export msmu report.tsv result.h5mu --params report.log.txt --fasta proteins.fasta
apb-export prolfqua report.tsv result.h5ad --params report.log.txt --annotation dataset.csv
apb-export proteopy report.tsv proteins.h5ad --params report.log.txt
apb-export alphapepttools report.tsv linked.h5mu --params report.log.txt
apb-export qfeatures report.tsv qfeatures.h5mu --params report.log.txt
apb-export qpx report.tsv result.qpx --params report.log.txt
```

msmu then runs unchanged on the result:

```python
import msmu as mm

mdata = mm.read_h5mu("result.h5mu")
mdata = mm.pp.add_filter(mdata, modality="psm", column="q_value", keep="lt", value=0.01)
mdata = mm.pp.apply_filter(mdata, modality="psm", on="var")
mdata = mm.pp.to_peptide(mdata, calculate_q=False)
```

## From Python

APB2 converts and `Exporter` exports, the shape of the other APB tools:

```python
from pathlib import Path

from apb2.api import ParseRuleCompiler
from apb_export.api import Exporter

exporter = Exporter("msmu")
parsed = ParseRuleCompiler.detect(
    Path("report.tsv"), Path("report.log.txt"), requested_levels=exporter.levels
).compile().parse()
exporter.write(exporter.export(parsed), Path("result.h5mu"))
```

`write` uses APB2's writer for HDF5 and a target-specific serializer for QPX. Stored APB results work through `apb2.api.read_parsed_levels` and the same exporter API. `Exporter.targets()` lists the packaged targets; `extension` names the suffix or `.qpx` directory.

## Install

apb-export is not on PyPI yet. Install from a checkout with its sibling packages beside it, as CI checks them out:

```bash
git clone https://github.com/anndata-omics-bridge/apb2.git
git clone https://github.com/anndata-omics-bridge/apb-catalog.git
git clone https://github.com/anndata-omics-bridge/apb-fasta.git
git clone https://github.com/anndata-omics-bridge/protein-fasta.git protein_fasta
git clone https://github.com/anndata-omics-bridge/prozor.git
git clone https://github.com/anndata-omics-bridge/apb-export.git
cd apb-export
uv sync --group dev
make check
```

`uv tool install --editable .` puts `apb-export` on `PATH`. The separate [APB Export UI](https://github.com/anndata-omics-bridge/apb-export-ui) owns the web application and its locked environment.

## Testing

Tests use synthetic vendor files only and install none of the target tools. The [Consumers workflow](https://github.com/anndata-omics-bridge/apb-export/blob/main/.github/workflows/consumers.yml) exports the synthetic examples and checks that msmu, prolfquapp, ProteoPy, AlphaPeptTools and QFeatures read them, msmu against its own `read_diann`.
