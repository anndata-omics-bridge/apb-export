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

| Target | Rule | File | APB2 levels | Layout |
| --- | --- | --- | --- | --- |
| msmu | [v0_4](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/msmu/v0_4/rules.json) | `.h5mu` | ion | long `psm` modality, one feature per run and precursor |
| prolfquapp | [v2_11](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/prolfqua/v2_11/rules.json) | `.h5ad` | ion, or peptide for MaxQuant; reads protein | wide; `uns["prolfquapp"]` for `LFQData_from_anndata` |
| ProteoPy | [v0_1](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/proteopy/v0_1/rules.json) | `.h5ad` | protein | wide; `sample_id` and `protein_id` beside their indexes |
| AlphaPeptTools | [v0_4](https://github.com/anndata-omics-bridge/apb-export/blob/main/src/apb_export/export_rules/documents/alphapepttools/v0_4/rules.json) | `.h5mu` | ion, peptide, protein | one wide modality per level, linked by id columns |

## msmu

msmu's own layout, as `msmu.read_diann` writes it:

- One `psm` modality; runs are observations
- One feature per observed (run, precursor) cell, id `<run>.<ProForma precursor>`
- `X`: block-diagonal sparse float32, the chosen abundance layer; APB2's primary layer by default
- `var`: `proteins`, `peptide`, `stripped_peptide`, `filename`, `charge`, `peptide_length`, `decoy`, `contaminant`, `PEP`, and `q_value` where reported
- `varm["search_result"]`: every APB2 feature column and every layer value at that cell
- `uns`: msmu's reader settings plus APB metadata; `uns["apb"]` holds the result's root records and the `export` record
- Features without a protein: left out and counted in the export record, since msmu's protein inference would group them as one unnamed protein

`peptide` is ProForma for every vendor (`AC[UNIMOD:4]K`). msmu's modification parser reads it, and `to_ptm` takes the tag as written, for example `"[UNIMOD:21]"`.

## prolfquapp

- Factors: `apb2 annotate` columns, from `--annotation`
- Confidence layers, each where the vendor reports it: `qValue` (precursor, per run), `pg_qValue` (protein group, per run), `pg_qValue_experiment` (protein group, library or experiment-wide), `pep`
- Features without a protein: left out and counted in the export record
- With `--fasta`: the leading protein's `fasta.id`, `description`, `IDcolumn`, `gene_name`, `protein_length` and `nr_tryptic_peptides`, as prolfquapp's FASTA annotation writes them
- Without a FASTA: `fasta.id`, `description` and `IDcolumn` repeat `protein_Id`

## ProteoPy and AlphaPeptTools

- ProteoPy with `--fasta`: apb-fasta annotates the protein groups and adds `in_fasta`, `gene_id`, `description` and `fasta_organisms`
- AlphaPeptTools: ion and peptide modalities carry `contaminant`; with `--fasta` also `in_fasta`, `fasta_proteins` and `fasta_organisms`
- Neither tool reads these FASTA columns; they are kept for inspection

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
- Retention time is not mapped to `rt`; it stays in `varm["search_result"]`
