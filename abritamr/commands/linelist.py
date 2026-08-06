import click
import pathlib
import json
import logging
import pandas as pd
import sys
from io import StringIO

from abritamr.utils import check_assembly, check_amrfinder, check_any2fasta, wrangle_species, output_results
# from abritamr.run_finder import run_amrf
# from abritamr.parse_finder import amrf2dict
# from abritamr.parse_reportable import add_abritamr_results
# from abritamr.parse_amrtype import get_amr_type
from abritamr.amr_report import summary
from abritamr.abritamr_logging import log



def linelist(
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

    
    if 'sample_id' in amr.columns.tolist() and 'species' in amr.columns.tolist() and 'abritamr_subclass' in amr.columns.tolist():
        simple = True if args.viewtype == 'compact' else False
        lines = []
        for sid in amr['sample_id'].unique().tolist():
            tmp = amr[amr['sample_id'] == sid]

            line = summary(
                results =tmp, 
                _format = args.format, 
                sid = sid,
                simple = simple, 
                genesonly = args.genesonly,
                minidentity = args.min_identity,
                mincoverage = args.min_coverage,
                )
            lines.append(line)
        linelist = pd.concat(lines)

        return linelist
    else:

        log.critical(f"It looks like your input file is not correctly configured. Please run abritamr scan and amr_status to generate the appropriate inut file.")
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
    
    