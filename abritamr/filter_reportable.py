
from dataclasses import dataclass, asdict
from cel import evaluate
import logging
import sys
import pandas as pd
from abritamr.logger import log
from abritamr.utils import get_refgenes


def construct_filter(result:dict, refgenes:pd.DataFrame):
    result['abritamr_priority_status'] = "not-reportable"
    result['criteria_id'] = ""
    result['criteria_version'] = ""
    result['abritamr_AMR_type'] = ""
    # refgenes = get_refgenes()
    refgenes = refgenes.fillna("")
    # log.info(f"Will extract criteria")
    # listofreportable = get_abritamr_configs(cfgtype = "reportable", cfgpath = cfgpath)

    row = refgenes[refgenes[result['amrfinder_accession_field']] == result['Closest reference accession']]
    if not row.empty:
        rw = row.to_dict(orient = "records")[0]
        sts = rw['priority_status']
        tpe = rw["amrtype"]
        cid = rw['criteria_id']
        civ = rw['criteria_version']
        status_criteria = rw["additional_status_criteria"]
        tpe_criteria = rw["additional_type_criteria"]
        if status_criteria != "":
            rpt = evaluate(status_criteria, result)
            if not rpt:
                sts = f"not-{sts}"
        if tpe_criteria != "":
            trpt = evaluate(tpe_criteria, result)
            if not trpt:
                tpe = ""
        result['abritamr_priority_status'] = sts
        result['criteria_id'] = cid
        result['criteria_version'] = civ
        result['abritamr_AMR_type'] = tpe

    else:
        log.warning(f"{result['Element symbol']} with accession {result['Closest reference accession']} could not be found.")
           
    return result


        