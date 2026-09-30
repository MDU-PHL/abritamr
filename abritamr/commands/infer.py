import click
import pathlib
import json
import logging
import pandas as pd
import sys
from io import StringIO

from abritamr.utils import wrangle_species
from abritamr.amr_infer import gdst, gdst_results_to_df_long, gdst_results_to_df_wide
from abritamr.logger import log


def do_gdst(
    amr: pd.DataFrame, reference_folder: str, dflt_result: str, reporttype: str
) -> pd.DataFrame:

    if "sample_id" in amr.columns.tolist() and "species" in amr.columns.tolist():
        # simple = True if args.viewtype == 'compact' else False
        lines = []
        for sid in amr["sample_id"].unique().tolist():
            tmp = amr[amr["sample_id"] == sid]
            sp = tmp["species"].iloc[0]
            line = gdst(
                results=tmp,
                species=sp,
                reference_folder=reference_folder,
                dflt_result=dflt_result,
            )
            lines.extend(line)
        if lines != []:
            if reporttype == "long":
                linelist = gdst_results_to_df_long(lines)
            elif reporttype == "wide":
                linelist = gdst_results_to_df_wide(lines)
            else:
                log.critical(
                    f"Report type {reporttype} is not recognized. Please use 'long' or 'wide'."
                )
                raise SystemExit(1)

            return linelist
        else:
            log.warning(f"There were no results possible")
            return pd.DataFrame()
    else:
        log.critical(
            f"It looks like your input file is not correctly configured. Please run abritamr scan to generate the appropriate inut file."
        )
        raise SystemExit(1)


def abritamr_gdst(args) -> dict:

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
    inferred = do_gdst(
        amr=amr,
        reference_folder=args.reference_folder,
        dflt_result=args.dflt_result,
        reporttype=args.reporttype,
    )
    return inferred
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

