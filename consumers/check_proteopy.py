"""ProteoPy accepts apb-export's protein-level file as proteomics data.

Inside the consumer image: ``python check_proteopy.py EXAMPLES``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
from proteopy.utils.anndata import check_proteodata

if __name__ == "__main__":
    adata = ad.read_h5ad(Path(sys.argv[1]) / "proteopy.h5ad")
    verdict = check_proteodata(adata)
    assert verdict == (True, "protein"), verdict
    print(f"ProteoPy reads {adata.n_obs} samples x {adata.n_vars} proteins")
