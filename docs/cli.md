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
| `msmu` | APB2's ion level as msmu's `psm` modality: one feature per run and precursor |
| `prolfqua` | APB2's ion level for prolfquapp; the annotation's columns become its factors |
| `proteopy` | APB2's protein level for ProteoPy, with `sample_id` and `protein_id` |
| `alphapepttools` | APB2's ion, peptide and protein levels as AlphaPeptTools modalities, each primary layer |

## Arguments

- `DATA`: vendor table or vendor-result directory
- `OUTPUT`: new `.h5ad` or `.h5mu` file; the command refuses to overwrite and refuses the wrong file type

## Options

| Option | Commands | Meaning |
| --- | --- | --- |
| `--params` | all | Vendor parameter file |
| `--software` | all | Restrict recognition to one vendor, as in `apb2 convert` |
| `--abundance` | `msmu`, `prolfqua`, `proteopy` | Layer of the target's level to become X; APB2's primary layer by default |
| `--annotation` | `prolfqua`, `proteopy`, `alphapepttools` | SDRF or prolfquapp-style CSV/TSV sample annotation, as `apb2 annotate` reads |
| `--strict` | all | Promote APB2 layer-contract warnings to errors |
| `--fasta` | all | FASTA files, or their protein-fasta database: check the peptides and annotate the protein groups first |

`--abundance` applies only to targets that write one level, so `alphapepttools` has none. With `--fasta`, msmu also flags the FASTA's contaminants.

The command reports catalogue answers that are not yet reviewed and stops on them.
