import click
import pathlib
import json
import logging
import pandas as pd
import sys
from io import StringIO

from abritamr.utils import check_assembly, check_amrfinder, check_any2fasta, wrangle_species, output_results,abritamr_status_columns
# from abritamr.run_finder import run_amrf
# from abritamr.parse_finder import amrf2dict
# from abritamr.parse_reportable import add_abritamr_results
# from abritamr.parse_amrtype import get_amr_type
# from abritamr.amr_report import summary
from abritamr.logger import log
from abritamr.parse_reportable import add_abritamr_results



def generate_output(amr:dict,cfgpath:str)-> dict:
    amr = add_abritamr_results(amr = amr,cfgpath=cfgpath)
    # add_abritamr_results(amr:dict, species: str= "", sid : str = "", cfgpath:str="") -> pd.DataFrame
    return amr

def amr_status(
    args
        ) -> dict:

    try:
        amr = pd.read_csv(args.amr)
        # amr = amr.to_dict(orient = "records")
        # log.info("Opened input file.")
    except:
        
        amrlist = args.amr.read().split('\n')
        dlm = "," if "," in amrlist[0] else "\t"
        amr = []
        for a in amrlist:
            amr.append(a.split(dlm))
        amr = pd.DataFrame(amr[1:], columns=amr[0])
        # amr = amr.to_dict(orient = "records")
    log.info(f"Will now determine status of the AMR genes detected.")

    if 'sample_id' in amr.columns.tolist() and 'species' in amr.columns.tolist()  and 'abritamr_subclass' in amr.columns.tolist():
        amr = amr.to_dict(orient = "records")
        amr = generate_output(amr = amr, cfgpath = args.reportable_config)
    else:

        log.critical(f"It looks like your input file is not correctly configured. Please run abritamr scan to generate the appropriate inut file.")
        raise SystemExit(1)
    # print(amr)
    abritamr_cols = abritamr_status_columns()
    amr = pd.DataFrame(amr)
    # amr['amrfinderplus_db_version'] = dbv
    amr = amr[abritamr_cols]
    

    return amr
