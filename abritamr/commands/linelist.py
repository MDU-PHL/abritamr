import pandas as pd

from abritamr.utils import (
    check_assembly,
    check_amrfinder,
    check_any2fasta,
    wrangle_species,
    output_results,
    abritamr_amrtype_columns,
)
from abritamr.amr_report import summary
from abritamr.logger import log


def linelist(args) -> dict:

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
        simple = True if args.viewtype == "compact" else False
        lines = []
        print(amr)
        for sid in amr["sample_id"].unique().tolist():
            tmp = amr[amr["sample_id"] == sid]

            line = summary(
                results=tmp,
                _format=args.format,
                sid=sid,
                simple=simple,
                genesonly=args.genesonly,
                minidentity=args.min_identity,
                mincoverage=args.min_coverage,
            )
            lines.append(line)
        if lines != []:
            linelist = pd.concat(lines)
        else:
            cl = abritamr_amrtype_columns()
            linelist = pd.DataFrame(columns=cl, data=[])
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
