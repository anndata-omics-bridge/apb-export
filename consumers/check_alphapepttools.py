"""AlphaPeptTools links apb-export's levels into one MuLink MuData.

Inside the consumer image: ``python check_alphapepttools.py EXAMPLES``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import alphapepttools as apt
import mudata as md

if __name__ == "__main__":
    mdata = md.read_h5mu(Path(sys.argv[1]) / "alphapepttools.h5mu")
    linked = apt.io.mulink_from_anndatas({name: mdata[name] for name in mdata.mod})
    assert "feature_mapping" in linked.varp, "no cross-level links"
    sizes = ", ".join(f"{name} {linked[name].n_vars}" for name in linked.mod)
    print(f"AlphaPeptTools links {sizes}")
