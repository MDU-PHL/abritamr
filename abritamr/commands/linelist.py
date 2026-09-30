import pandas as pd

from abritamr.utils import (
    abritamr_amrtype_columns,
)
from abritamr.amr_report import summary
from abritamr.logger import log


def generate_linelist(
    amr: pd.DataFrame,
    _format: str,
    viewtype: str,
    genesonly: bool,
    min_identity: float,
    min_coverage: float,
) -> pd.DataFrame:

    if (
        "sample_id" in amr.columns.tolist()
        and "species" in amr.columns.tolist()
        and "abritamr_subclass" in amr.columns.tolist()
    ):
        simple = True if viewtype == "compact" else False
        lines = []
        for sid in amr["sample_id"].unique().tolist():
            tmp = amr[amr["sample_id"] == sid]

            line = summary(
                results=tmp,
                _format=_format,
                sid=sid,
                simple=simple,
                genesonly=genesonly,
                minidentity=min_identity,
                mincoverage=min_coverage,
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

        linelist = generate_linelist(
            amr=amr,
            _format=args.format,
            viewtype=args.viewtype,
            genesonly=args.genesonly,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
        )

        return linelist
