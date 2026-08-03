from abritamr.criteria import get_abritamr_configs
from dataclasses import dataclass, asdict
from cel import evaluate
import logging

logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def filter_string(crt:dict, dlm:str = " && ")-> str:

    fl = []

    for k in crt:
        if f"{crt[k]}" != "None" and k != "exception" and k != "action":
            if k not in ["abritamr_class","abritamr_subclass","gene"]:
                f = f"{k} in {crt[k]}"
            else:
                lst = eval(crt[k]) if not isinstance(crt[k],list) else crt[k]
                tmp = []
                for l in lst:
                    tmp.append(f"'{l}' in {k}")
                f = ' || '.join(tmp)
                
            fl.append(f)
    return f'{dlm}'.join(fl)

def extract_criteria(criteria):

    crt = {k: str(v) for k, v in asdict(criteria).items() if k != "action"}
    return crt

def construct_filter(result:dict, cfgpath:str):
    rpt = "not-reportable"
    # log.info(f"Will extract criteria")
    listofreportable = get_abritamr_configs(cfgtype = "reportable", cfgpath = cfgpath)
    for criteria in listofreportable:
        crt = extract_criteria(criteria = criteria)
        action = {k: str(v) for k, v in asdict(criteria).items() if k == "action"}
        primary_filter = filter_string(crt)
        res = evaluate(primary_filter,result)
        if res:
            rpt = action['action']
        if 'exception' in crt and crt['exception'] != "None":
            sfilterres = set()
            for xcptn in eval(crt['exception']):
                flt = filter_string(xcptn,dlm = " && ")
                
                resxptn = evaluate(flt, result)
                if resxptn:
                    sfilterres.add(xcptn["action"])
            if sfilterres:
                rpt = "|".join(list(sfilterres))
    
    return rpt


        