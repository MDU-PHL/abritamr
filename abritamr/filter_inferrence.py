 # for now - will need to be cleverer/cleaner/better 
from dataclasses import dataclass, asdict
from cel import evaluate
import pandas as pd
# print(criterias)


def filter_string(ifr:dict)-> str:

    fl = []
    # print(ifr)
    for k in ifr:
        # print(k)
        # print(ifr[k])
        if f"{ifr[k]}" != "None" and k != "exception" and k != "result" and k != "drugname":
            # print(ifr[k])
            # if k  in ["genus","species"]:
            #     f = f"{k} == {ifr[k]}"
            # else:
                # print(type(ifr[k])) 
            dlm = " || " 
            key = k.split("_")
            # print(ifr[k])
            cond = '_'.join(key[:2]) if len(key) > 3 else key[0]
            # print(cond)
            comp = key[3] if len(key) == 4 else "in"
            # print(comp)
            # val = ifr[k].split(",")
            lst = eval(ifr[k]) if not isinstance(ifr[k],list) else ifr[k]
            # print(lst)
            tmp = []
            for l in lst:
                if comp == "equals":
                    tmp.append(f"'{l}' == {cond}")
                
                else:
                    tmp.append(f"'{l}' in {cond}")
            if comp == "and":
                dlm = " && "
            f = f'{dlm}'.join(tmp)
            f =f"({f})"
                # print(f)
            fl.append(f)
    # # print(f'{dlm}'.join(fl))
    return f' && '.join(fl)

def extract_criteria(criteria:dict) -> dict:

    crt = {k: str(v) for k, v in asdict(criteria).items() if k not in ["result"]}
    # print(crt)
    return crt

def combine_results( result: list, species:str="", genus:str="" ) -> dict:
    res = {}
    to_test = {
        'abritamr_class':[],
        'abritamr_subclass':[],
        'gene':[],
       
    }
    # print(result)
    for row in result:
        for col in to_test:
            if col in row:
                to_test[col].append(row[col])
        
    
    for k in to_test:
        # print(to_test[k])
        val = ','.join(to_test[k]) if to_test[k] != [] else "None"
        
        res[k] = val
    # print(res)
    res['species']= species
    res['genus']=genus
    return res

def get_mechs(crt:dict,result:dict, primary_filter:str="", mechs:list=[]) -> str:

    conds = [i for i in crt if ('class' in i  or 'gene' in i) and (crt[i] != "None" or not crt[i])]
   
    cols = eval(crt[conds[0]]) if not isinstance(crt[conds[0]], list) else crt[conds[0]]
   
    for row in result:
        for k in row:
            if k in conds[0] and row[k] in primary_filter:
                mechs.append(row['gene'])
            elif k in conds[0]:
                for c in cols:
                    if c in row[k]:
                        mechs.append(row['gene'])
   
    return mechs

def construct_filter(result:list, species:str="",genus:str="",sid:str="",default:str = "Susceptible"):
  
    listofinfer = get_abritamr_configs(cfgtype = "infer")
    
    to_test = combine_results(result = result, species = species, genus = genus)

    inferred = {}
    for criteria in listofinfer:
        rpt = f"{default}"
        # print(criteria)
        mechs = []
        crt = extract_criteria(criteria = criteria)
        dr = crt['drugname']
        resistance = {k: str(v) for k, v in asdict(criteria).items() if k == "result"}
        primary_filter = filter_string(crt)
        
        res = evaluate(primary_filter,to_test)
        if res:
            rpt = resistance['result']
            mechs = get_mechs(crt = crt, result = result, primary_filter = primary_filter, mechs = mechs)
            
        if 'exception' in crt and crt['exception'] != "None":
            print(f"Found an exception rule")
            sfilterres = set()
            mechs = []
            for xcptn in eval(crt['exception']):
                flt = filter_string(xcptn)
                
                resxptn = evaluate(flt, to_test)
                
                if resxptn:
                    mechs = get_mechs(crt = xcptn, result = result, primary_filter = flt, mechs=mechs)
                    sfilterres.add(xcptn["result"])
            if sfilterres:
                rpt = "|".join(list(sfilterres))
            
        mechs = ";".join(sorted(mechs)) if mechs != [] else "No mechanisms identified"
        
        inferred[f"{dr} - mechanisms"] = mechs
        inferred[f"{dr} - interpretation"] = rpt

    
    return inferred


        

