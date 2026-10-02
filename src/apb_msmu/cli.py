"""``apb-msmu``: vendor output in, msmu-ready ``.h5mu`` out."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated

from apb_catalog.resolver import UnresolvedField
from cyclopts import App, Parameter
from loguru import logger

from apb_msmu.api import convert

app = App(
    name="apb-msmu",
    help="Convert vendor output into the MuData msmu's readers build, in one step.",
)


@app.default
def run(
    data: Annotated[Path, Parameter(help="Vendor table or vendor-result directory")],
    output: Annotated[Path, Parameter(help="New .h5mu file for msmu.read_h5mu")],
    /,
    *,
    params: Annotated[Path | None, Parameter(help="Vendor parameter file")] = None,
    software: Annotated[
        str | None, Parameter(help="Restrict recognition to one vendor, as in apb2 convert")
    ] = None,
    abundance: Annotated[
        str | None, Parameter(help="Layer to become msmu's X; APB2's primary layer by default")
    ] = None,
    strict: Annotated[
        bool, Parameter(negative=False, help="Promote APB2 layer-contract warnings to errors")
    ] = False,
) -> int:
    """Write OUTPUT, the msmu MuData of the vendor output DATA."""
    if output.suffix != ".h5mu":
        logger.error(f"output must end in .h5mu: {output}")
        return 1
    if output.exists():
        logger.error(f"refusing to overwrite {output}")
        return 1
    try:
        mdata = convert(data, params, software=software, abundance=abundance, strict=strict)
    except (OSError, ValueError) as error:
        logger.error(str(error))
        return 1
    except UnresolvedField as error:
        # Unknown is not missing: the catalogue has not reviewed this exact vendor rule, so it
        # cannot say whether a q-value or PEP exists. Guessing would export a silent gap.
        logger.error(
            f"apb-catalog cannot answer for this vendor rule ({error}); re-review the "
            "catalogues against this APB2 revision"
        )
        return 1
    output.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(output)
    psm = mdata["psm"]
    apb = psm.uns["apb"]
    logger.info(
        f"wrote {output}: psm {psm.n_obs} runs x {psm.n_vars} features; "
        f"X={apb['abundance_layer']} q_value={apb['q_value_layer']} PEP={apb['pep_layer']}"
    )
    return 0


def main() -> None:
    """Run the command line and exit with its status."""
    sys.exit(app(sys.argv[1:]))
