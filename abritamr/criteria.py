from dataclasses import dataclass
from typing import Optional,Literal
import pandas as pd
import json



@dataclass(frozen=True)
class BaseInfer:
    result: Literal['Susceptible','Resistant','Intermediate', 'Not-reported'] = "Not-reported"
    drugname: Optional[list] = None
    gene_present: Optional[list] = None
    gene_missing: Optional[list] = None
    mech_present: Optional[list] = None
    mech_missing: Optional[list] = None
    drugsubclass_present_contains: Optional[list] = None
    drugsubclass_present_equals: Optional[list] = None
    drugsubclass_present_and: Optional[list] = None
    drugsubclass_present_or: Optional[list] = None
    species: Optional[list] = None
    genus: Optional[list] = None

    def __post_init__(self):
        if self.result not in ['Susceptible','Resistant','Intermediate','Not-reported']:
            raise ValueError(f"result must be one of 'Susceptible', 'Resistant', 'Intermediate' or 'Not-reported'")
        if not self.drugname :
            raise ValueError(f"You must supply a 'drugname'")
        if not self.drugsubclass_present_contains and not self.drugsubclass_present_equals and not self.drugsubclass_present_and and not self.drugsubclass_present_or and not self.gene_present and not self.mech_present:
            raise ValueError(f"You must supply one a drugsubclass, drugclass, gene or other mechanism")
        
@dataclass(frozen=True)
class Infer:
    result: Literal['Susceptible','Resistant','Intermediate', 'Not-reported'] = "Not-reported"
    drugname: Optional[list] = None
    gene_present: Optional[list] = None
    gene_missing: Optional[list] = None
    mech_present: Optional[list] = None
    mech_missing: Optional[list] = None
    drugsubclass_present_contains: Optional[list] = None
    drugsubclass_present_equals: Optional[list] = None
    drugsubclass_present_and: Optional[list] = None
    drugsubclass_present_or: Optional[list] = None
    species: Optional[list] = None
    genus: Optional[list] = None
    exception: BaseInfer  = None
    def __post_init__(self):
        if self.result not in ['Susceptible','Resistant','Intermediate','Not-reported']:
            raise ValueError(f"result must be one of 'Susceptible', 'Resistant', 'Intermediate' or 'Not-reported'")
        if not self.drugname :
            raise ValueError(f"You must supply a 'drugname'")
        if not self.drugsubclass_present_contains and not self.drugsubclass_present_equals and not self.drugsubclass_present_and and not self.drugsubclass_present_or and not self.gene_present and not self.mech_present:
            raise ValueError(f"You must supply one a drugsubclass, drugclass, gene or other mechanism")
        
        

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
    

def wrangle_lists_infer(infer:dict) -> dict:

    listkeys = [
        "gene",
        "species",
        "genus","drugsubclass_present_contains", 
        "drugsubclass_present_equals",
        "drugsubclass_present_and",
        "drugsubclass_present_or",
        "gene_present",
        "gene_missing",
        "mech_present",
        "mech_missing"
        ]

    for l in listkeys:
        if l in infer:
            lst = []
            # print(infer[l])
            if infer[l] != "":
                lst = infer[l].split(",")
                lst = [i.strip() for i in lst]
            infer[l] = lst
    return infer


def construct_inference(resistance:dict) -> Infer:
    # check_required(criteria = criteria)
    # check_one_ofs(criteria = criteria)
    infer = wrangle_lists_infer(resistance)
    # print(criteria)
    if 'exception' in infer and infer['exception'] != []:
        exceptioninfer = []
        for cr in infer['exception']:
            ct = wrangle_lists_infer(cr)
            bif = BaseInfer(**ct)
            
            exceptioninfer.append(bif)
        infer['exception'] = exceptioninfer
    
    ifr = Infer(**infer)

    return ifr


def wrangle_lists_criteria(criteria:dict) -> dict:

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
    criteria = wrangle_lists_criteria(criteria)
    # print(criteria)
    if 'exception' in criteria and criteria['exception'] != []:
        exceptioncriteria = []
        for cr in criteria['exception']:
            ct = wrangle_lists_criteria(cr)
            bcr = BaseCriteria(**ct)
            
            exceptioncriteria.append(bcr)
        criteria['exception'] = exceptioncriteria
    
    ctr = Criteria(**criteria)

    return ctr



resistance = [
    {
        "drugname":"Colistin",
        "species":"Salmonella enterica",
        "drugsubclass_present_contains":"Colistin",
        "result":"Resistant"
       
    },
    {
        "drugname":"Aminoglycosides (RMT)",
        "species":"Salmonella enterica",
        "drugsubclass_present_contains":"Aminoglycosides (Ribosomal methyltransferase)",
        "result":"Resistant",
    },
    {
        "drugname":"Chloramphenicol",
        "species":"Salmonella enterica",
        "drugsubclass_present_contains":"phenicol",
        "result":"Resistant",
    },
    {
        "drugname":"Trim-Sulpha",
        "species":"Salmonella enterica",
        "drugsubclass_present_and":"Trimethoprim,Sulfathiazole",
        "result":"Resistant",
    },
    {
        "drugname":"Trimethoprim",
        "species":"Salmonella enterica",
        "drugsubclass_present_equals":"Trimethoprim",
        "result":"Resistant",
    },
    {
        "drugname":"Sulfathiazole",
        "drugsubclass_present_equals":"Sulfathiazole",
        "result":"Resistant",
    },
    {
        "drugname":"Ciprofloxacin",
        "species":"Salmonella enterica",
        "drugsubclass_present_contains":"quinolone",
        "result":"Resistant",
    },
    {
        "drugname":"Tetracycline",
        "species":"Salmonella enterica",
        "drugsubclass_present_equals":"Tetracycline",
        "result":"Resistant",
    },
    {
        "drugname":"Streptomycin",
        "species":"Salmonella enterica",
        "drugsubclass_present_or":"Streptomycin,Spectinomycin",
        "result":"Resistant",
    },
    {
        "drugname":"Kanamycin",
        "species":"Salmonella enterica",
        "drugsubclass_present_or":"Kanamycin,Aminoglycosides (Ribosomal methyltransferase)",
        "result":"Resistant",
    },
    {
        "drugname":"Gentamicin",
        "species":"Salmonella enterica",
        "drugsubclass_present_or":"Gentamicin,Aminoglycosides (Ribosomal methyltransferase)",
        "result":"Resistant",
    },
    {
        "drugname":"Azithromycin",
        "species":"Salmonella enterica",
        "drugsubclass_present_or":"Macrolide,Lincosamide/Macrolide/Streptogramin,Azithromycin",
        "result":"Resistant",
    },
    {
        "drugname":"Meropenem",
        "species":"Salmonella enterica",
        "drugsubclass_present_contains":"Carbapenemase",
        "result":"Resistant",
    },
    {
        "drugname":"Cefotaxime (AmpC)",
        "species":"Salmonella enterica",
        "drugsubclass_present_contains":"AmpC",
        "result":"Resistant",
    },
    {
        "drugname":"Cefotaxime (ESBL)",
        "species":"Salmonella enterica",
        "drugsubclass_present_equals":"ESBL",
        "result":"Resistant",
    },
    {
        "drugname":"Ampicillin",
        "species":"Salmonella enterica",
        "drugsubclass_present_or":"Beta-lactam,ESBL,AmpC,Ampicillin",
        "result":"Resistant",
    },

]

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
listofinfer = []
for c in criterias:

    c1 = construct_criteria(c)
    # print(c1)
    listofreportable.append(c1)


for r in resistance:

    r1 = construct_inference(r)
    print(r1)
    listofinfer.append(r1)
