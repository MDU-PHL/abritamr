import pandas as pd

from abritamr.utils import (
    abritamr_status_columns,
)

from abritamr.logger import log
from abritamr.parse_reportable import add_abritamr_results


def generate_output(amr: dict, catalog: str) -> dict:
    amr = add_abritamr_results(amr=amr, catalog=catalog)
    # add_abritamr_results(amr:dict, species: str= "", sid : str = "", cfgpath:str="") -> pd.DataFrame
    return amr


def amr_status(args) -> dict:
    log.info("Going to try to open input.")
    try:
        amr = pd.read_csv(args.amr)
        # # print(amr)
        # amr = amr.to_dict(orient = "records")
        # log.info("Opened input file.")
    except Exception as e:
        amrlist = args.amr.read().split("\n")
        dlm = "," if "," in amrlist[0] else "\t"
        amr = []
        for a in amrlist:
            amr.append(a.split(dlm))
        # # print(amr[0])
        amr = pd.DataFrame(amr[1:], columns=amr[0])
    # # print(amr)
    # amr = amr.to_dict(orient = "records")
    log.info(f"Will now determine status of the AMR genes detected.")

    if (
        "sample_id" in amr.columns.tolist()
        and "species" in amr.columns.tolist()
        and "abritamr_subclass" in amr.columns.tolist()
    ):
        amr = amr.to_dict(orient="records")
        amr = generate_output(amr=amr, catalog=args.reference_catalog)
    else:
        log.critical(
            f"It looks like your input file is not correctly configured. Please run abritamr scan to generate the appropriate inut file."
        )
        raise SystemExit(1)
    # # print(amr)
    abritamr_cols = abritamr_status_columns()
    amr = pd.DataFrame(amr)
    # amr['amrfinderplus_db_version'] = dbv
    amr = amr[abritamr_cols]

    return amr
