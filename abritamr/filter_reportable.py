from criteria import listofreportable # for now - will need to be cleverer/cleaner/better 
from dataclasses import dataclass, asdict
from cel import evaluate

# print(criterias)


def filter_string(crt:dict, dlm:str = " && ")-> str:

    fl = []

    for k in crt:
        # print(k)
        # print(crt[k])
        if f"{crt[k]}" != "None" and k != "exception" and k != "action":
            if k not in ["abritamr_class","abritamr_subclass","gene"]:
                f = f"{k} in {crt[k]}"
            else:
                # print(type(crt[k])) 
                lst = eval(crt[k]) if not isinstance(crt[k],list) else crt[k]
                tmp = []
                for l in lst:
                    tmp.append(f"'{l}' in {k}")
                f = ' || '.join(tmp)
                # print(f)
            fl.append(f)
    # print(f'{dlm}'.join(fl))
    return f'{dlm}'.join(fl)

def extract_criteria(criteria):

    crt = {k: str(v) for k, v in asdict(criteria).items() if k != "action"}
    # print(crt)
    return crt

def construct_filter(result:dict):
    rpt = "not-reportable"
    
    for criteria in listofreportable:
       
        
        crt = extract_criteria(criteria = criteria)
        # print(crt)
        action = {k: str(v) for k, v in asdict(criteria).items() if k == "action"}
        primary_filter = filter_string(crt)
        # print(action)
        # print(primary_filter)
        # print(result['abritamr_subclass'])
        res = evaluate(primary_filter,result)
        if res:
            rpt = action['action']
        if 'exception' in crt and crt['exception'] != "None":
            # print(crt['exception'])
            sfilterres = set()
            for xcptn in eval(crt['exception']):
                # print(xcptn)
                flt = filter_string(xcptn,dlm = " && ")
                # print(flt)
                # print(result)
                resxptn = evaluate(flt, result)
                if resxptn:
                    sfilterres.add(xcptn["action"])
            if sfilterres:
                rpt = "|".join(list(sfilterres))

    return rpt


        



results = [
    {"data":{"abritamr_subclass": "Carbapenemase","gene":"blabla", "species":"Salmonella enterica", "genus": "Salmonella"}, "res": "reportable"}, #this should be report
    {"data":{"abritamr_subclass": "Carbapen","gene":"blabla", "species":"Salmonella enterica","genus": "Salmonella"}, "res":"not-reportable"}, # this will not report}
    {"data":{"abritamr_subclass": "Carbapenemase (MBL)", "gene": "bla1", "species":"Salmonella enterica", "genus": "Salmonella"},"res":"reportable"}, # this should report
    {"data":{"abritamr_subclass": "Carbapenemase (MBL)", "gene": "bla1", "species":"Stenotrophomonas maltophilia", "genus":"Stenotrophomonas"},"res":"not-reportable"}, # this not should report
    {"data": {"abritamr_subclass": "Carbapenemase (MBL)", "gene": "bla2", "species":"Stenotrophomonas maltophilia","genus":"Stenotrophomonas"}, "res":"reportable"}, # this  should report
    {"data": {"abritamr_subclass": "Linezolid/Phenicol", "gene": "bla2", "species":"Staphylococcus aureus","genus":"Staphylococcus"}, "res":"reportable"}, # this  should report

]

for result in results:
    res = construct_filter(result["data"])
    if res == result["res"]:
        print("Success")
    else:
        print(f"Something is wrong, {result['data']} should return {result['res']}")

    # break