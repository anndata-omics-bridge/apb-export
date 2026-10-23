# Targets and layouts

Each target's layout is declared in its export rule, not coded. One engine checks every source before it builds anything.

## Export rules

Rules follow APB2's rule schema read in reverse:

- `name`: the target's field
- `source`: an APB field every vendor shares, such as a ProForma column, a role, the run key or a catalogued meaning
- `"level": "protein"` on a catalogued source: read through each variable's protein group, named whole or by its leading protein; a cell takes its group's value in its run
- `software` on a level block: the vendor that block is for
- `.h5ad` rules: list their levels in order and write the first one the result's software and levels allow

The published [JSON Schemas](https://github.com/anndata-omics-bridge/apb-export/tree/main/src/apb_export/export_rules/documents/_schema) describe the documents.

Every packaged target's file is written by APB2's `write_parsed_levels`, `.apb.json` sidecar included:

- each key stays a column and names its index, a String key as written
- var names carry no level prefix
- several levels share one observation axis; MaxQuant's experiments are aligned onto its raw files when each experiment is one raw file
- a rule's `uns` entries are written beside APB's `uns["apb"]`
- a rule's `feature_links` become one directed relation in the MuData's `varp`, under the link's name: a nonzero (u, v), 1.0, maps feature u to the target level's feature v, as mulink reads it

| Target | Rule | File | APB2 levels | Layout |
| --- | --- | --- | --- | --- |
| msmu | [v0_4](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/msmu/v0_4/rules.json) | `.h5mu` | ion | `psm` modality, precursors by runs, experimental |
| prolfquapp | [v2_11](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/prolfqua/v2_11/rules.json) | `.h5ad` | ion, or peptide for MaxQuant; reads protein | wide; `uns["prolfquapp"]` for `LFQData_from_anndata` |
| ProteoPy | [v0_1](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/proteopy/v0_1/rules.json) | `.h5ad` | protein | wide; `sample_id` and `protein_id` beside their indexes |
| AlphaPeptTools | [v0_4](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/alphapepttools/v0_4/rules.json) | `.h5mu` | ion, peptide, protein | one wide modality per level, linked by id columns |

## msmu

msmu's `psm` modality, written wide. Experimental: msmu's own readers write one feature per run and precursor, and msmu's tools may expect that.

- One `psm` modality; runs are observations, precursors are features, id `<ProForma precursor>`; the keys `filename` and `precursor` stay columns and name the indices
- `X`: dense float32, the chosen abundance layer; APB2's primary layer by default; a precursor absent from a run is NaN
- `layers`: each run's `q_value` and `PEP`, where the vendor reports them
- `var`: `proteins`, `peptide`, `stripped_peptide`, `charge`, `peptide_length`, `decoy`, `contaminant`, and `PEP` and `q_value` where reported, each the precursor's lowest value across runs (`how: "layer_minimum"`)
- `uns`: `apb_export_status: "experimental"`, msmu's reader settings plus APB metadata; `uns["apb"]` holds the result's root records and the `export` record
- Features without a protein: left out and counted in the export record, since msmu's protein inference would group them as one unnamed protein

`peptide` is ProForma for every vendor (`AC[UNIMOD:4]K`). msmu's modification parser reads it, and `to_ptm` takes the tag as written, for example `"[UNIMOD:21]"`.

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
- AlphaPeptTools: ion and peptide modalities carry `contaminant`; with `--fasta` also `in_fasta`, `fasta_proteins` and `fasta_organisms`
- AlphaPeptTools on MaxQuant: runs named by raw file; with fractions, several raw files per experiment, only precursors are written and the export record lists the other levels as skipped
- Neither tool reads these FASTA columns; they are kept for inspection

## Plugin targets

An installed package adds a target without changing apb-export: it names a package holding the target's export rule under the `apb_export.targets` entry-point group.

```toml
[project.entry-points."apb_export.targets"]
toy = "toy_plugin.export_rule"
```

- Layout: `<package>/<version>/rules.json`, one version per target, the same schema as the packaged documents
- `Exporter.targets()` lists the target and `Exporter("toy")` applies it; a target named like a packaged one, or by two plugins, is refused
- Command: `apb-export toy DATA OUTPUT` with `--params`, `--software`, `--fasta`, `--sdrf` and `--strict`; no `--abundance`, so the rule's primary layers become X
- Rule features no packaged target uses yet: `qualifiers` on a catalogued source, `obs` sources naming SDRF or annotation columns, the table's own `uns` entries, and `feature_links`
- The web app offers only the packaged targets it describes

## Confidence, by meaning

- `q_value`: library q-value when match-between-runs filled it, else experiment-wide, else per run; the order msmu's DIA-NN reader uses
- `PEP`: the vendor's posterior error probability
- Lookup: apb-catalog's `identification_confidence` and `proteobench_entrapment` catalogues
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

- A result without the level a target reads is refused: msmu needs an `ion` level (Sage's `lfq.tsv` is peptide-level), ProteoPy a `protein` level
- Label-free only: a result keyed by more than one observation column, as multiplexed designs are, is refused
- `decoy` is apb2's `apb_Decoy`; most exports carry no decoys, so use msmu's `calculate_q=False` unless the vendor kept them
- Group separators other than `;` pass through unchanged
- Retention time is not mapped to `rt`
- msmu's file has no per-run `filename` column: msmu tools that group features by file see precursors, not identifications
