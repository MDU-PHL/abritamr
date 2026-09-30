"""Assign AMR reporting status to scan results."""

import pandas as pd

from abritamr.utils import (
    abritamr_status_columns,
)

from abritamr.logger import log
from abritamr.parse_reportable import add_abritamr_results


def generate_output(amr: dict, catalog: str) -> dict:
    """Add reportable AMR classes and typing criteria to detected hits."""
    amr = add_abritamr_results(amr=amr, catalog=catalog)
    return amr


def do_typing(amr: pd.DataFrame, reference_catalog: str) -> pd.DataFrame:
    """Validate scan results and annotate each hit with its AMR status."""
    log.info(f"Will now determine status of the AMR genes detected.")

    if (
        "sample_id" in amr.columns.tolist()
        and "species" in amr.columns.tolist()
        and "abritamr_subclass" in amr.columns.tolist()
    ):
        amr = amr.to_dict(orient="records")
        amr = generate_output(amr=amr, catalog=reference_catalog)
    else:
        log.critical(
            f"It looks like your input file is not correctly configured. Please run abritamr scan to generate the appropriate inut file."
        )
        raise SystemExit(1)
    abritamr_cols = abritamr_status_columns()
    amr = pd.DataFrame(amr)
    amr = amr[abritamr_cols]

    return amr


def amr_status(args) -> dict:
    """Read scan output, determine AMR status, and return the annotated table."""
    log.info("Going to try to open input.")
    try:
        amr = pd.read_csv(args.amr)
    except Exception as e:
        amrlist = args.amr.read().split("\n")
        dlm = "," if "," in amrlist[0] else "\t"
        amr = []
        for a in amrlist:
            amr.append(a.split(dlm))
        amr = pd.DataFrame(amr[1:], columns=amr[0])

    amr = do_typing(amr=amr, reference_catalog=args.reference_catalog)

    return amr
