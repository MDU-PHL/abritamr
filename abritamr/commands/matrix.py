import pathlib
import pandas as pd

from abritamr.amr_matrix import summary
from abritamr.logger import log
from abritamr.utils import abritamr_matrix_columns


def matrix(args) -> dict:

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
                facet=args.facet,
                sid=sid,
                minidentity=args.min_identity,
                mincoverage=args.min_coverage,
                refgenes=args.reference_catalog,
            )
            lines.append(line)
        try:
            linelist = pd.concat(lines)
        except ValueError:
            mcols = abritamr_matrix_columns(
                refgenes=args.reference_catalog, group=args.facet
            )
            linelist = pd.DataFrame(columns=mcols, data=[])
        return linelist
    else:
        log.critical(
            f"It looks like your input file is not correctly configured. Please run abritamr scan and amr_status to generate the appropriate inut file."
        )
        raise SystemExit(1)

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
