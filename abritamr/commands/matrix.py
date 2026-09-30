"""Generate AMR presence/absence matrices from typed results."""

import pathlib
import pandas as pd

from abritamr.amr_matrix import summary
from abritamr.logger import log
from abritamr.utils import abritamr_matrix_columns


def make_matrix(
    amr: pd.DataFrame,
    facet: str,
    min_coverage: float,
    min_identity: float,
    reference_catalog: str,
) -> pd.DataFrame:
    """Create one matrix per sample using the requested AMR facet."""
    if (
        "sample_id" in amr.columns.tolist()
        and "species" in amr.columns.tolist()
        and "abritamr_subclass" in amr.columns.tolist()
    ):
        lines = []
        for sid in amr["sample_id"].unique().tolist():
            tmp = amr[amr["sample_id"] == sid]

            line = summary(
                results=tmp,
                facet=facet,
                sid=sid,
                minidentity=min_identity,
                mincoverage=min_coverage,
                refgenes=reference_catalog,
            )
            lines.append(line)
        try:
            linelist = pd.concat(lines)
        except ValueError:
            mcols = abritamr_matrix_columns(refgenes=reference_catalog, group=facet)
            linelist = pd.DataFrame(columns=mcols, data=[])
        return linelist
    else:
        log.critical(
            f"It looks like your input file is not correctly configured. Please run abritamr scan and amr_status to generate the appropriate inut file."
        )
        raise SystemExit(1)


def matrix(args) -> dict:
    """Read typed AMR input and generate a matrix from command options."""
    try:
        amr = pd.read_csv(args.amr)
        # amr = amr.to_dict(orient = "records")
        # log.info("Opened input file.")
    except:
        amrlist = args.amr.read().split("\n")
        dlm = "," if "," in amrlist[0] else "\t"
        amr = []
        for a in amrlist:
            amr.append(a.split(dlm))
        amr = pd.DataFrame(amr[1:], columns=amr[0])
    if amr.empty:
        log.critical(f"No results available for generating matrix output")
        raise SystemExit(1)
    matrix = make_matrix(
        amr=amr,
        facet=args.facet,
        min_coverage=args.min_coverage,
        min_identity=args.min_identity,
        reference_catalog=args.reference_catalog,
    )

    return matrix
    # results:pd.DataFrame,
    # _format:str="csv",
    # species:str="",
    # genus:str="",
    # simple:bool=False,
    # sid:str="abritamr",
    # genesonly:bool = False,
    # minidentity:float = 90,
    # mincoverage:float = 90,
    # outname:str = "abritamr_report"
