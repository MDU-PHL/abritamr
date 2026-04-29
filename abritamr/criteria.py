from dataclasses import dataclass
from typing import Optional,Literal
import pandas as pd
import json



@dataclass(frozen=True)
class BaseInfer:
    result: Literal['Susceptible','Resistant','Intermediate', 'Not-reported'] = "Not-reported"
    drug: Optional[list] = None
    gene_present: Optional[list] = None
    gene_missing: Optional[list] = None
    mech_present: Optional[list] = None
    mech_missing: Optional[list] = None
    drugsubclass_present: Optional[list] = None
    drugsubclass_missing: Optional[list] = None
    species: Optional[list] = None

    def __post_init__(self):
        if self.result not in ['Susceptible','Resistant','Intermediate','Not-reported']:
            raise ValueError(f"result must be one of 'Susceptible', 'Resistant', 'Intermediate' or 'Not-reported'")
    

@dataclass(frozen=True)
class BaseCriteria:
    action: Literal['reportable','not-reportable']
    drugclass: Optional[list] = None
    drugsubclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    gene: Optional[list]  = None
    # exception: Optional[list]  = None
    def __post_init__(self):
        if self.action not in ['reportable','not-reportable']:
            raise ValueError(f"action must be one of reportable or not-reportable")
        if not self.drugclass and not self.drugsubclass and not self.gene:
            raise ValueError(f"You must supply one of 'drugclass', 'drugsubclass' or 'gene'")
        
@dataclass(frozen=True)
class Criteria:
    action: Literal['reportable','not-reportable']
    drugclass: Optional[list] = None
    drugsubclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
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
                lst = [i.strip() for i in lst]
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

criterias= [

    {
        "drugsubclass": "Carbapenemase",
        "action":"reportable"
    },
    {
        "drugsubclass":"ESBL (KPC variant)",
        "action":"reportable"
    },
    {
        "drugsubclass":"Aminoglycosides (Ribosomal methyltransferase)",
        "action":"reportable"
    },
    {
        "drugsubclass":"Colistin",
        "action":"reportable",
        "exception":[
            {
                "drugsubclass":"Colistin",
                "action":"not-reportable",
                "gene":"mcr9"
            }
        ]
    },
    {
        "drugsubclass":"Carbapenemase (OXA-51 family)",
        "action":"reportable",
        "exception":[
            {
                "action":"not-reportable",
                "drugsubclass":"Carbapenemase (OXA-51 family)",
                "species":"Acinetobacter baumannii,Acinetobacter calcoaceticus,Acinetobacter nosocomialis,Acinetobacter pittii,Acinetobacter baumannii complex"
            }
        ]
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
    },
    {
        "drugsubclass": "ESBL",
        "genus": "Salmonella,Shigella",
        "action":"reportable"
    },
    {
        "drugsubclass":"AmpC",
        "genus":"Salmonella,Shigella",
        "action":"reportable",
        "exception":[
            {
                "drugsubclass":"AmpC",
                "genus":"Shigella",
                "action":"not-reportable",
                "gene":"blaEC"
            }
        ]
    },
    {
        "drugsubclass":"Vancomycin",
        "gene" : "vanA,vanB,vanC,vanD,vanE,vanE,vanG,vanL,vanM,vanN",
        "action": "reportable"
    },
    {
        "drugsubclass":"Methicillin",
        "gene" : "mecI,mecR",
        "action": "reportable"
    },
    {
        "drugsubclass":"Oxazolidinone",
        "species" : "Staphylococcus aureus,Staphylococcus argenteus",
        "action": "reportable"
    },
    {
        "drugsubclass":"Oxazolidinone",
        "genus" : "Enterococcus",
        "action": "reportable"
    },
    {
        "drugsubclass":"Linezolid",
        "species" : "Staphylococcus aureus,Staphylococcus argenteus",
        "action": "reportable"
    },
    {
        "drugsubclass":"Linezolid",
        "genus" : "Enterococcus",
        "action": "reportable"
    }
    ]

listofreportable = []
for c in criterias:

    c1 = construct_criteria(c)
    # print(c1)
    listofreportable.append(c1)
