"""ProteoPy accepts apb-export's protein and ion files, and its own DIA-NN reader agrees.

Inside the consumer image: ``python check_proteopy.py EXAMPLES``.

The ion file is compared with what ProteoPy's own reader makes of the same synthetic report at
``aggr_level="Precursor.Id"``: the same ions per sample, each with the same value and protein_id.
ProteoPy's reader for DIA-NN >= 1.9.1 reads Parquet with a ``Decoy`` column, so the check writes
the report once as Parquet, with no decoys; its TSV reader for 1.x cannot aggregate by
``Precursor.Id``. Ions are matched by sample, protein_id and value rather than by id, as the two
spell an ion differently (``PEPTIDEK2`` and ``PEPTIDEK/2``); the synthetic values are unique
within a sample.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from proteopy.read import diann
from proteopy.utils.anndata import check_proteodata


def _cells(adata: ad.AnnData) -> dict[str, list[tuple[str, float]]]:
    """Per sample, the sorted (protein_id, value) pairs of the cells holding a value."""
    values = np.asarray(adata.X, dtype=float)
    proteins = [str(protein) for protein in adata.var["protein_id"]]
    return {
        str(sample): sorted(
            (protein, float(value))
            for protein, value in zip(proteins, row, strict=True)
            if not np.isnan(value)
        )
        for sample, row in zip(adata.obs_names, values, strict=True)
    }


def check_ions(examples: Path) -> None:
    """The ion file is ProteoPy peptide data and equals ProteoPy's own reading of the report."""
    ours = ad.read_h5ad(examples / "proteopy_ion.h5ad")
    verdict = check_proteodata(ours)
    assert verdict == (True, "peptide"), verdict
    report = pd.read_csv(examples / "report.tsv", sep="\t")
    report["Decoy"] = 0
    with tempfile.TemporaryDirectory() as folder:
        parquet = Path(folder) / "report.parquet"
        report.to_parquet(parquet)
        theirs = diann(parquet, aggr_level="Precursor.Id", version="1.9.1")
    assert ours.n_vars == theirs.n_vars, (ours.n_vars, theirs.n_vars)
    assert sorted(map(str, ours.obs_names)) == sorted(map(str, theirs.obs_names))
    assert _cells(ours) == _cells(theirs), "values or protein_id differ from ProteoPy's reader"
    print(
        f"ProteoPy reads {ours.n_obs} samples x {ours.n_vars} ions, "
        "equal to its own DIA-NN reader at Precursor.Id"
    )


if __name__ == "__main__":
    examples = Path(sys.argv[1])
    adata = ad.read_h5ad(examples / "proteopy.h5ad")
    verdict = check_proteodata(adata)
    assert verdict == (True, "protein"), verdict
    print(f"ProteoPy reads {adata.n_obs} samples x {adata.n_vars} proteins")
    check_ions(examples)
