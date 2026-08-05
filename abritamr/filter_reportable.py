from abritamr.criteria import get_abritamr_configs
from dataclasses import dataclass, asdict
from cel import evaluate
import logging
import sys
from abritamr.abritamr_logging import log

# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
# handler = logging.StreamHandler(sys.stderr)

# log = logging.getLogger(__name__)
# log.addHandler(handler)
# log.setLevel(logging.DEBUG)

def construct_filter(result:dict, cfgpath:str):
    result['abritamr_AMR_reporting_status'] = "not-reportable"
    result['criteria_id'] = ""
    result['criteria_version'] = ""
    result['abritamr_AMR_type'] = ""
    
    # log.info(f"Will extract criteria")
    listofreportable = get_abritamr_configs(cfgtype = "reportable", cfgpath = cfgpath)
    for criteria in listofreportable:
        
        crt = asdict(criteria)
        res = evaluate(crt['criteria'],result)
        if res:
            result['abritamr_AMR_reporting_status'] = crt['status']
            result['abritamr_AMR_type'] = crt["amrtype"] if crt["amrtype"] else ""
            result['criteria_id'] = crt['criteria_id']
            result['criteria_version'] = crt['criteria_version']
            
    return result


        