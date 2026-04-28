from criteria import listofcriteria # for now - will need to be cleverer/cleaner/better 
from dataclasses import dataclass, asdict
from cel import evaluate

# print(criterias)


def filter_string(crt:dict, dlm:str = " && ")-> str:

    fl = []

    for k in crt:
        if crt[k] and crt[k] != "None" and k != "exception" and k != "action":
            f = f"{k} in {crt[k]}"
            fl.append(f)
    return f'{dlm}'.join(fl)

def extract_criteria(criteria):

    crt = {k: str(v) for k, v in asdict(criteria).items() if k != "action"}
    return crt

def construct_filter(result:dict):
    rpt = "not-reportable"
    for criteria in listofcriteria:
       
        # print(crt)
        crt = extract_criteria(criteria = criteria)
        action = {k: str(v) for k, v in asdict(criteria).items() if k == "action"}
        primary_filter = filter_string(crt)
        print(primary_filter)
        res = evaluate(primary_filter,result)
        if res:
            rpt = action['action']
        if 'exception' in crt and crt['exception'] != "None":
            # print(crt['exception'])
            sfilterres = set()
            for xcptn in eval(crt['exception']):
                
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
    {"data":{"drugsubclass": "Carbapenemase","gene":"blabla", "species":"Salmonella enterica"}, "res": "reportable"}, #this should be report
    {"data":{"drugsubclass": "Carbapenemases","gene":"blabla", "species":"Salmonella enterica"}, "res":"not-reportable"}, # this will not report}
    {"data":{"drugsubclass": "Carbapenemase (MBL)", "gene": "bla1", "species":"Salmonella enterica"},"res":"reportable"}, # this should report
    {"data":{"drugsubclass": "Carbapenemase (MBL)", "gene": "bla1", "species":"Stenotrophomonas maltophilia"},"res":"not-reportable"}, # this not should report
    {"data": {"drugsubclass": "Carbapenemase (MBL)", "gene": "bla2", "species":"Stenotrophomonas maltophilia"}, "res":"reportable"}, # this  should report

]

for result in results:
    res = construct_filter(result["data"])
    if res == result["res"]:
        print("Success")
    else:
        print(f"something is wrong, {result['data']} should return {result['res']}")

    # break