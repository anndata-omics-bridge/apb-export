# Targets and layouts

Each target's layout is declared in its export rule, not coded. One engine checks every source before it builds anything.

## Export rules

Rules follow APB2's rule schema read in reverse:

- `name`: the target's field
- `source`: an APB field every vendor shares: the run key, the primary layer, or a catalogued meaning `{"catalog": <set>, "kind": ..., "location": ...}`. APB's own columns come from apb-catalog's `apb` set, such as `{"catalog": "apb", "kind": "peptidoform_ion", "location": "var"}` for `ProForma_ion` or kind `protein_assignment` for the column a rule marks with that role, so no rule names an APB column. `"enforce": "MUST"` makes the set's MUST fields binding and logs a warning for a missing SHOULD field; prolfqua's intensity sets it
- `"level": "protein"` on a catalogued source: read through each variable's protein group, named whole or by its leading protein; a cell takes its group's value in its run; a var column's value, the group's, repeats in every run
- `software` on a level block: the vendor that block is for
- `.h5ad` rules: list their levels in order and write the first one the result's software and levels allow
- `.csv` rules: choose their level as `.h5ad` rules do and write it as a long table: one row per run and feature that some layer holds, obs columns, var columns and one column per layer, category codes written as the labels the parsing rule declares, missing as `NA`, no sidecar
- `.txt` rules: choose their level as `.h5ad` rules do and write it as a Perseus table, described [below](#perseus)
- `transform` on a layer: `log2`, applied after `missing_values`; a value that is not positive becomes missing
- `how: "short_name"`: per row, the first input holding a value, cut to its first group member and, for `db|ACC|NAME`, to `NAME`; inputs the result lacks, such as FASTA columns without `--fasta`, are passed over; `append` follows after `_`
- `rename` on an `all_values` entry: annotation columns whose name matches a `pattern` are written as text, categories sorted, under `re.sub(pattern, replacement, name)`, the first matching rename winning; `identifier` lower-cases that name and turns every character outside `a-z0-9_` into `_`
- `name_by_kind` on a layer: the layer's name per `abundance` kind of the layer it resolved to, such as `{"maxlfq": "LFQ intensity"}`; any other kind keeps `name`

The published [JSON Schemas](https://github.com/anndata-omics-bridge/apb-export/tree/main/src/apb_export/export_rules/documents/_schema) describe the documents.

HDF5 target files are written by APB2's `write_parsed_levels`, `.apb.json` sidecar included:

- each key stays a column and names its index, a String key as written
- var names carry no level prefix
- several levels share one observation axis; MaxQuant's experiments are aligned onto its raw files when each experiment is one raw file
- a rule's `uns` entries are written beside APB's `uns["apb"]`
- a rule's `feature_links` become one directed relation in the MuData's `varp`, under the link's name: a nonzero (u, v), 1.0, maps feature u to the target level's feature v, as mulink reads it

| Target | Rule | File | APB2 levels | Layout |
| --- | --- | --- | --- | --- |
| msmu | [v0_4](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/msmu/v0_4/rules.json) | `.h5mu` | peptidoform_ion | `psm` modality, precursors by runs, experimental |
| MSstats | [v4_20](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/msstats/v4_20/rules.json) | `.csv` | peptidoform_ion | long table, one row per precursor and run, MSstats' column names |
| Perseus | [v2](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/perseus/v2/rules.json) | `.txt` | protein | wide, proteinGroups.txt's names, `#!{Type}` and `#!{C:…}` header rows |
| prolfquapp | [v2_11](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/prolfqua/v2_11/rules.json) | `.h5ad` | peptidoform_ion, or peptide for MaxQuant; reads protein | wide; `uns["prolfquapp"]` for `LFQData_from_anndata` |
| ProteoPy | [v0_1](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/proteopy/v0_1/rules.json) | `.h5ad` | protein | wide; `sample_id` and `protein_id` beside their indexes |
| exploreDE | [v1_0](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/explorede/v1_0/rules.json) | `.h5ads` directory | protein, peptide, peptidoform_ion, those present | `exploreDE/1.0` without differential expression, one file per level the result holds: `protein.h5ad`, `peptide.h5ad`, `precursor.h5ad`; log2 `intensities` with vendor zeros missing, linear `intensities_raw`, X a copy of `intensities`; SDRF `factor value[…]` and prolfquapp `group`/`bait` columns become `factor_<name>` with sorted categories, other annotation columns pass unprefixed; `feature_id` is the protein group, or `<protein>_<peptide>` and `<protein>_<ProForma ion>`; `label_gene_symbol` is the leading protein's FASTA gene name, else its entry name or accession, with the peptide or precursor appended; `uns` holds exploreDE's identity keys and empty `de`, `enrichment` and `enrichment_support` |
| AlphaPeptTools | [v0_4](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/alphapepttools/v0_4/rules.json) | `.h5mu` | peptidoform_ion, peptide, protein | one wide modality per level, linked by id columns |
| QFeatures | [v1_22](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/qfeatures/v1_22/rules.json) | `.h5mu` | peptidoform_ion, peptide, protein | AlphaPeptTools' layout; assays linked by `addAssayLink` on id columns |
| QPX | [v1_1](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/qpx/v1_1/rules.json) | `.qpx` directory | peptidoform_ion and/or protein | QPX 1.1 Parquet tables, metadata and export manifest |

## QPX

```bash
apb-export qpx report.tsv result.qpx --params report.log.txt
multiqc --qpx-plugin result.qpx -o report
```

- Existing primary precursor quantities become `apb.feature.parquet`; existing primary protein quantities become `apb.pg.parquet`. Either level suffices; no quantities are aggregated or inferred
- Feature quantities use nested `intensities`; protein quantities use flat `label`, `intensity`, and `grouped_runs`, as pMultiQC 0.0.48 reads. Identifiers are opaque integers per table; complete protein-group membership is retained even when distinct groups share a leading protein
- Missing measurements produce no quantity row; measured zeroes remain zero. Float32 intensities follow QPX 1.1; infinite or overflowing values are refused
- `rt` is read only through APB Catalog's `chromatography` set, kind `retention_time`, with original units `minute` or `second`. Minutes become seconds; missing RT remains null. An unknown unit is an unresolved catalogue answer, so the exporter never assumes one
- Unknown spectral values remain typed nulls; unknown scan references and peptide-count lists are empty lists. No spectral measurements or peptide counts are fabricated
- `apb.run.parquet` and `apb.sample.parquet` provide technical identities keyed by run, with optional supplied observation annotations. QPX requires organism and tissue strings, so missing values use the explicit marker `not available`; no organism or grouping column is required from the user
- `export.json` records versions, source mappings and table sizes. HDF5 targets continue to use APB2's sidecar; the QPX directory has its own manifest
- The schema snapshot comes from [qpx 1.1.6](https://github.com/bigbio/qpx/tree/v1.1.6), with its Apache-2.0 license included. pMultiQC's QPX reader does not provide organism-based feature grouping

Report generation belongs to the [apb-pmultiqc](https://github.com/anndata-omics-bridge/apb-pmultiqc) application. [pMultiQC PR #751](https://github.com/bigbio/pmultiqc/pull/751), merged on 29 October 2026, fixes upstream RT plots for missing/all-null RT; until a release contains it, the released reader (0.0.48) can log a plot exception for such valid inputs. Consumer checks require populated reports and no plot errors.

## msmu

msmu's `psm` modality, written wide. Experimental: msmu's own readers write one feature per run and precursor, and msmu's tools may expect that.

- One `psm` modality; runs are observations, precursors are features, id `<ProForma precursor>`; the keys `filename` and `precursor` stay columns and name the indices
- `X`: dense float32, the chosen abundance layer; APB2's primary layer by default; a precursor absent from a run is NaN
- `layers`: `q_value`, the library or experiment-wide q-value repeated in every run where the vendor reports one, else each run's; each run's `PEP`, where the vendor reports it
- `var`: `proteins`, `peptide`, `stripped_peptide`, `charge`, `peptide_length`, `decoy`, `contaminant`, and `PEP` and `q_value` where reported, each the precursor's lowest value across runs (`how: "layer_minimum"`)
- `uns`: `apb_export_status: "experimental"`, msmu's reader settings plus APB metadata; `uns["apb"]` holds the result's root records and the `export` record
- Features without a protein: left out and counted in the export record, since msmu's protein inference would group them as one unnamed protein

`peptide` is ProForma for every vendor (`AC[UNIMOD:4]K`). msmu's modification parser reads it, and `to_ptm` takes the tag as written, for example `"[UNIMOD:21]"`.

## MSstats

```bash
apb-export msstats report.tsv msstats.csv --params report.log.txt --annotation design.tsv
```

```r
input <- read.csv("msstats.csv")
processed <- MSstats::dataProcess(input)
```

- One row per precursor and run that some layer holds; MSstats reads the columns by name, so their order does not matter
- `Run`: the run key; `Condition` and `BioReplicate`: from the sample annotation, through apb-catalog: a prolfquapp-style table's `group…` (or `bait…`, `Experiment…`) and `BioReplicate…` (or `subject…`) columns, an SDRF's `factor value[…]` and `characteristics[biological replicate]`. Without them the export stops
- `ProteinName`: the protein group; `PeptideSequence`: the ProForma peptidoform, modifications included; `PrecursorCharge`: the ion's charge; `PrecursorId`: the ProForma ion, the table's feature key
- `FragmentIon`, `ProductCharge`: `NA`, as MSstatsConvert writes precursor-level input; `IsotopeLabelType`: `L`
- `Intensity`: the `abundance` catalog's intensity without normalization, since MSstats normalizes itself, else APB2's primary layer; zeros written as `NA`; `--abundance` names another layer
- `DetectionQValue`, `PEP`: run-level identification q-value and PEP where the tool reports them; no global or library q-values, and no rows are filtered
- Features without a protein: left out and counted in the export record

## Perseus

```bash
apb-export perseus txt perseus.txt --params mqpar.xml --annotation samples.tsv
```

Perseus: Generic matrix upload, or read with `perseuspy.pd.read_perseus("perseus.txt")`.

- One row per protein group; names as MaxQuant's `proteinGroups.txt`, so Perseus' usual steps run without renaming
- `Protein IDs`: the protein group; `Gene names`: with `--fasta`; `Q-value`: the group's experiment-wide, else library, q-value, where reported
- `Reverse`, `Potential contaminant`: `+` or empty, from APB2's decoy and contaminant flags, contaminants also from the FASTA check
- `<prefix> <sample>`: APB2's primary layer, or `--abundance`'s; prefix by its `abundance` kind: `LFQ intensity` for MaxLFQ, `iBAQ` for iBAQ, else `Intensity`; zeros and missing values written as `NaN`
- `#!{Type}`: sample columns `E` (main), numbers `N`, flags `C`, text `T`
- `#!{C:<column>}`: one categorical row per annotation column, values under the sample columns
- No sidecar

## prolfquapp

- Intensity (`X`): the `abundance` catalogue's intensity without normalization, since prolfquapp normalizes itself, such as DIA-NN's `Precursor_Quantity`; APB2's primary layer where the catalogue has none or cannot pick one; `--abundance` names another layer
- Factors: `apb2 annotate` columns, from `--annotation`
- Confidence layers, each where the vendor reports it: `qValue` (precursor, per run), `pg_qValue` (protein group, per run), `pg_qValue_experiment` (protein group, library or experiment-wide, the same in every run), `pep`
- Features without a protein: left out and counted in the export record
- With `--fasta`: the leading protein's `fasta.id`, `description`, `IDcolumn`, `gene_name`, `protein_length` and `nr_tryptic_peptides`, as prolfquapp's FASTA annotation writes them
- Without a FASTA: `fasta.id`, `description` and `IDcolumn` repeat `protein_Id`

## ProteoPy and AlphaPeptTools

- ProteoPy: `sample_id` and `protein_id` name the indices and stay columns
- ProteoPy with `--fasta`: apb-fasta annotates the protein groups and adds `in_fasta`, `gene_id`, `description` and `fasta_organisms`
- AlphaPeptTools: peptidoform_ion and peptide modalities carry `contaminant`; with `--fasta` also `in_fasta`, `fasta_proteins` and `fasta_organisms`
- AlphaPeptTools on MaxQuant: runs named by raw file; with fractions, several raw files per experiment, only precursors are written and the export record lists the other levels as skipped
- Neither tool reads these FASTA columns; they are kept for inspection

## QFeatures

The modalities `precursors`, `peptides` (where the vendor reports a peptide level) and `proteins` become QFeatures assays. Each finer assay links to the next coarser one on that one's key: `sequence` for peptides, `proteins` for protein groups. MaxQuant peptides name their protein group when proteinGroups.txt is converted with them.

MuData reads the file from GitHub (`BiocManager::install("ilia-kats/MuData")`) until a Bioconductor release carries [MuData#11](https://github.com/ilia-kats/MuData/pull/11).

```r
library(QFeatures)
mae <- MuData::readH5MU("qfeatures.h5mu")
# as(x, "SummarizedExperiment") drops the rowData; X arrives unnamed; categorical columns as factors
as_assay <- function(x) {
  features <- rowData(x)
  for (column in names(features)) {
    if (is.factor(features[[column]])) features[[column]] <- as.character(features[[column]])
  }
  SummarizedExperiment(assays = list(intensity = assay(x)), rowData = features)
}
# Annotations sit in each modality's colData; the MuData's own colData is empty
qf <- QFeatures(lapply(experiments(mae), as_assay), colData = colData(mae[[1]]))
qf <- addAssayLink(qf, from = "precursors", to = "peptides", varFrom = "sequence", varTo = "sequence")
qf <- addAssayLink(qf, from = "peptides", to = "proteins", varFrom = "proteins", varTo = "proteins")
qf <- aggregateFeatures(qf, "precursors", fcol = "proteins", name = "aggregated", fun = colSums, na.rm = TRUE)
```

Without a peptide level, as for DIA-NN, link `precursors` to `proteins` on `proteins`. The Consumers workflow runs [check_qfeatures.R](https://github.com/anndata-omics-bridge/apb-export/blob/main/consumers/r/check_qfeatures.R) on a DIA-NN and a MaxQuant example.

## Plugin targets

An installed package adds a target without changing apb-export: it names a package holding the target's export rule under the `apb_export.targets` entry-point group.

```toml
[project.entry-points."apb_export.targets"]
toy = "toy_plugin.export_rule"
```

- Layout: `<package>/<version>/rules.json`, one version per target, the same schema as the packaged documents
- `Exporter.targets()` lists the target and `Exporter("toy")` applies it; a target named like a packaged one, or by two plugins, is refused
- Command: `apb-export toy DATA OUTPUT` with `--params`, `--software`, `--fasta`, `--sdrf` and `--strict`; no `--abundance`, so the rule's primary layers become X
- Rule features no packaged target uses yet: `obs` sources naming SDRF or annotation columns, and `feature_links`
- The web app offers only the packaged targets it describes

## Confidence, by meaning

- `q_value`: library q-value (var) when it holds a nonzero value, else experiment-wide (var), else per run; the first two repeat in every run; the order msmu's DIA-NN reader uses
- `PEP`: the vendor's posterior error probability
- Lookup: apb-catalog's `confidence` and `proteobench_entrapment` catalogues
- Not reported: no `q_value` column, so msmu's filter cannot act on an invented value; `PEP` is NaN, because msmu's `to_peptide` requires the column
- Unreviewed rule: apb-catalog answers "unknown" and the command stops; unknown is not missing

## Protein groups

Written as msmu's readers write them:

- Members split on `;`
- Accession taken from a UniProt `db|ACC|NAME` entry
- Members spelled `contam_`, `Cont_` or `CON__` respelled `Cont_ACC`
- `contaminant`: set when apb2 marked `apb_Contaminant`, the software flags or adds the contaminant, or apb-fasta matched a FASTA contaminant

`--fasta` runs that FASTA check after conversion; without a FASTA only the vendor's own flags are known.

## Limits

- A result without the level a target reads is refused: msmu needs an `peptidoform_ion` level (Sage's `lfq.tsv` is peptide-level), ProteoPy a `protein` level
- Label-free only: a result keyed by more than one observation column, as multiplexed designs are, is refused
- `decoy` is apb2's `apb_Decoy`; most exports carry no decoys, so use msmu's `calculate_q=False` unless the vendor kept them
- Group separators other than `;` pass through unchanged
- Retention time is mapped to QPX `rt` through the chromatography catalogue; HDF5 targets do not map it
- msmu's file has no per-run `filename` column: msmu tools that group features by file see precursors, not identifications
