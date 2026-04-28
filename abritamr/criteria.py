from dataclasses import dataclass
from typing import Optional
import pandas as pd
import json

CRITERIA_COLS = {
    'required': ['action'],
    'one_of': ['drugclass','drugsubclass','gene'],
    'optional':['species','drug','abx','exception']
}

@dataclass(frozen=True)
class BaseCriteria:
    action: str   # "report" or "not_report"
    drugclass: Optional[list] = None
    drugsubclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    drug: Optional[str]  = None
    abx: Optional[str]  = None
    gene: Optional[list]  = None
    # exception: Optional[list]  = None
    
@dataclass(frozen=True)
class Criteria:
    action: str   # "report" or "not_report"
    drugclass: Optional[list] = None
    drugsubclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    drug: Optional[str]  = None
    abx: Optional[str]  = None
    gene: Optional[list]  = None
    exception: BaseCriteria  = None
    



def check_required(criteria:dict)-> bool:
    print(f"Checking that required fields are supplied for criteria.")
    for req in CRITERIA_COLS['required']:
        if req not in criteria:
            print(f"You must have {req} in your criteria. Please try again.")
            raise SystemExit(1)
        if criteria[req] == "":
            print(f"You must supply a valid value for {req}. Please try again.")
            raise SystemExit(1)
    print(f"Looks like all the required fields are present and not empty.")
    return True

def check_one_ofs(criteria:dict) -> bool:
    print(f"Checking that the criteria has at least one of {' | '.join(CRITERIA_COLS['one_of'])}")

    ones = []
    for o in CRITERIA_COLS['one_of']:
        if o in criteria and criteria[o] != "":
            ones.append(o)
    if ones == []:
        print(f"You must supply at least one of {' | '.join(CRITERIA_COLS['one_of'])}. Please try again.")
        raise SystemExit(1)
    print(f"Looks like you have correctly supplied {' | '.join(ones)}")
    return True

def wrangle_lists(criteria:dict) -> dict:

    listkeys = ["gene","species","genus","drugclass","drugsubclass"]

    for l in listkeys:
        if l in criteria:
            lst = []
            # print(l)
            if criteria[l] != "":
                lst = criteria[l].split(",")
            criteria[l] = lst
    return criteria

def construct_criteria(criteria:dict) -> Criteria:
    check_required(criteria = criteria)
    check_one_ofs(criteria = criteria)
    criteria = wrangle_lists(criteria)
    # print(criteria)
    if 'exception' in criteria and criteria['exception'] != []:
        exceptioncriteria = []
        for cr in criteria['exception']:
            ct = wrangle_lists(cr)
            bcr = BaseCriteria(**ct)
            exceptioncriteria.append(bcr)
        criteria['exception'] = exceptioncriteria
    
    ctr = Criteria(**criteria)

    return ctr
# elif i == "Carbapenemase (MBL)" and species == "Stenotrophomonas maltophilia":
#                             # if species is "Stenotrophomonas maltophilia" don't report blaL1
#                             genes_reported.extend([g for g in genes if not g.startswith("blaL1")])
#                             genes_not_reported.extend([g for g in genes if g.startswith("blaL1")])

criterias= [
    {
        "drugsubclass": "Carbapenemase",
        "action":"reportable"
    },
    {
        "drugsubclass": "Carbapenemase (MBL)",
        "action":"reportable",
        "exception":
            [
                {
                    "species": "Stenotrophomonas maltophilia",
                    "drugsubclass": "Carbapenemase (MBL)",
                    "gene":"bla1",
                    "action":"not-reportable"
                }
            ]
    }
            ]

listofcriteria = []
for c in criterias:

    c1 = construct_criteria(c)
    # print(c1)
    listofcriteria.append(c1)

# abritamr = Abritamr(RULES)
# abritamr.run(data, output, etc)