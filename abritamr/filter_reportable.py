
from dataclasses import dataclass, asdict
from cel import evaluate
import logging
import sys
import pandas as pd
from abritamr.logger import log
from abritamr.utils import get_refgenes
from abritamr.cel_functions import create_cel_context, evaluate_rule

def construct_filter(result:dict, refgenes:pd.DataFrame):
    # # print(result)
    result['abritamr_priority_status'] = "not-reportable"
    result['criteria_id'] = ""
    result['criteria_version'] = ""
    result['abritamr_AMR_type'] = ""
    # refgenes = get_refgenes()
    refgenes = refgenes.fillna("")
    # log.info(f"Will extract criteria")
    # listofreportable = get_abritamr_configs(cfgtype = "reportable", cfgpath = cfgpath)
    # # print(result['Closest reference accession'])
    row = refgenes[refgenes['abritamr_accession_key'].str.contains(result['Closest reference accession'], na=False)]
    if not row.empty:
        rw = row.to_dict(orient = "records")[0]
        sts = rw['priority_status']
        tpe = rw["amrtype"]
        cid = rw['criteria_id']
        civ = rw['criteria_version']
        status_criteria = rw["additional_status_criteria"]
        tpe_criteria = rw["additional_type_criteria"]
        # # print({'row': result})
        data = {'row': result}
        ctx = create_cel_context(data = result, name = 'row')
        # print(data)
        if status_criteria != "":
            # # print(ctx)
            # # print(f"Evaluating status criteria: {status_criteria}")
            # # print(f"Context: {ctx}")
            
            rpt = evaluate_rule(status_criteria, ctx)
            # # print(f"Result of status criteria evaluation: {rpt}")
            if not rpt:
                sts = f"not-{sts}"
        if tpe_criteria != "":
            # print(f"Evaluating type criteria: {tpe_criteria}")
            trpt = evaluate_rule(tpe_criteria, ctx)
            if not trpt:
                tpe = ""
        result['abritamr_priority_status'] = sts
        result['criteria_id'] = cid
        result['criteria_version'] = civ
        result['abritamr_AMR_type'] = tpe

    else:
        log.warning(f"{result['Element symbol']} with accession {result['Closest reference accession']} could not be found.")
           
    return result


        