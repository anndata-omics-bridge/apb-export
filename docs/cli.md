# Command-line interface

```bash
apb-export --help
```

One subcommand per target. Each takes the vendor output and a new output file, converts with APB2, and writes the file the target opens.

```bash
apb-export TARGET DATA OUTPUT [OPTIONS]
```

| Command | Writes |
| --- | --- |
| `msmu` | APB2's peptidoform_ion level as msmu's `psm` modality: one feature per run and precursor |
| `prolfqua` | APB2's peptidoform_ion level for prolfquapp; the annotation's columns become its factors |
| `proteopy` | APB2's protein level for ProteoPy, with `sample_id` and `protein_id` |
| `alphapepttools` | APB2's peptidoform_ion, peptide and protein levels as AlphaPeptTools modalities, each primary layer |
| `qfeatures` | APB2's peptidoform_ion, peptide and protein levels as QFeatures assays, linked by id columns |
| `qpx` | Existing peptidoform_ion and/or protein primary quantities as QPX 1.1 Parquet in a new `.qpx` directory |
| plugin target | Each target an installed plugin adds; see [plugin targets](targets.md#plugin-targets) |

## Arguments

- `DATA`: vendor table or vendor-result directory
- `OUTPUT`: new `.h5ad`/`.h5mu` file or `.qpx` directory; the command refuses to overwrite and refuses the wrong file type

## Options

| Option | Commands | Meaning |
| --- | --- | --- |
| `--params` | all | Vendor parameter file |
| `--software` | all | Restrict recognition to one vendor, as in `apb2 convert` |
| `--abundance` | `msmu`, `prolfqua`, `proteopy` | Layer of the target's level to become X; by default the rule's choice: APB2's primary layer, or for prolfqua the catalogued intensity without normalization |
| `--annotation` | `prolfqua`, `proteopy`, `alphapepttools`, `qfeatures`, `qpx` | SDRF or prolfquapp-style CSV/TSV sample annotation, as `apb2 annotate` reads |
| `--strict` | all | An effectively empty layer is an error, not left out with a finding |
| `--fasta` | all | FASTA files, or their protein-fasta database: check the peptides and annotate the protein groups first |
| `--sdrf` | plugin targets | SDRF-Proteomics sample table: every run must match one of its rows, whose columns fill the sample metadata |

`--abundance` applies only to targets that write one level, so `alphapepttools` and `qfeatures` have none. With `--fasta`, msmu also flags the FASTA's contaminants.

The command reports catalogue answers that are not yet reviewed and stops on them.
