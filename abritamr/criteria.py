from dataclasses import dataclass
from typing import Optional,Literal
import pandas as pd
import json

CRITERIA_COLS = {
    'required': ['action'],
    'one_of': ['drugclass','drugsubclass','gene'],
    'optional':['species','drug','abx','exception']
}

@dataclass(frozen=True)
class BaseCriteria:
    action: Literal['reportable','not-reportable']
    drugclass: Optional[list] = None
    drugsubclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    drug: Optional[str]  = None
    abx: Optional[str]  = None
    gene: Optional[list]  = None
    # exception: Optional[list]  = None
    def __post_init__(self):
        if self.action not in ['reportable','not-reportable']:
            raise ValueError(f"Action must be one of reportable or not-reportable")
        if not self.drugclass and not self.drugsubclass and not self.gene:
            raise ValueError(f"You must supply one of 'drugclass', 'drugsubclass' or 'gene'")
        
@dataclass(frozen=True)
class Criteria:
    action: Literal['reportable','not-reportable']
    drugclass: Optional[list] = None
    drugsubclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    drug: Optional[str]  = None
    abx: Optional[str]  = None
    gene: Optional[list]  = None
    exception: BaseCriteria  = None
    def __post_init__(self):
        if self.action not in ['reportable','not-reportable']:
            raise ValueError(f"Action must be one of reportable or not-reportable")
        if not self.drugclass and not self.drugsubclass and not self.gene:
            raise ValueError(f"You must supply one of 'drugclass', 'drugsubclass' or 'gene'")
    



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
    # check_required(criteria = criteria)
    # check_one_ofs(criteria = criteria)
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
