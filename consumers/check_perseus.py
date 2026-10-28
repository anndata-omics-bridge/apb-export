"""perseuspy reads apb-export's Perseus table as Perseus' own files: sample columns typed main.

Inside the consumer image: ``python check_perseus.py EXAMPLES``.
"""

from __future__ import annotations

import sys
from pathlib import Path

from perseuspy import pd

if __name__ == "__main__":
    table = pd.read_perseus(str(Path(sys.argv[1]) / "perseus.txt"))
    names = table.columns.get_level_values("Column Name")
    samples = [name for name in names if name.startswith(("Intensity ", "LFQ intensity ", "iBAQ "))]
    assert samples, f"no sample columns in {list(names)}"
    assert all(table[name].dtypes.iloc[0].kind == "f" for name in samples), table.dtypes
    assert "condition" in table.columns.names, table.columns.names
    print(f"perseuspy reads {len(table)} proteins x {len(samples)} samples, by {table.columns.names}")
