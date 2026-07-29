from dataclasses import dataclass
from typing import Optional,Literal
import pandas as pd
import json
import pathlib


@dataclass(frozen=True)
class BaseInfer:
    result: Literal['Susceptible','Resistant','Intermediate', 'Not-reported'] = "Not-reported"
    drugname: Optional[list] = None
    gene_present: Optional[list] = None
    gene_missing: Optional[list] = None
    mech_present: Optional[list] = None
    mech_missing: Optional[list] = None
    abritamr_subclass_present_contains: Optional[list] = None
    abritamr_subclass_present_equals: Optional[list] = None
    abritamr_subclass_present_and: Optional[list] = None
    abritamr_subclass_present_or: Optional[list] = None
    abritamr_subclass_present_in: Optional[list] = None
    species: Optional[list] = None
    genus: Optional[list] = None

    def __post_init__(self):
        if self.result not in ['Susceptible','Resistant','Intermediate','Not-reported']:
            raise ValueError(f"result must be one of 'Susceptible', 'Resistant', 'Intermediate' or 'Not-reported'")
        if not self.drugname :
            raise ValueError(f"You must supply a 'drugname'")
        if not self.abritamr_subclass_present_contains and not self.abritamr_subclass_present_equals and not self.abritamr_subclass_present_in and not self.abritamr_subclass_present_or and not self.gene_present and not self.mech_present:
            raise ValueError(f"You must supply one a abritamr_subclass, abritamr_class, gene or other mechanism")
        
@dataclass(frozen=True)
class Infer:
    result: Literal['Susceptible','Resistant','Intermediate', 'Not-reported'] = "Not-reported"
    drugname: Optional[list] = None
    gene_present: Optional[list] = None
    gene_missing: Optional[list] = None
    mech_present: Optional[list] = None
    mech_missing: Optional[list] = None
    abritamr_subclass_present_contains: Optional[list] = None
    abritamr_subclass_present_equals: Optional[list] = None
    abritamr_subclass_present_and: Optional[list] = None
    abritamr_subclass_present_or: Optional[list] = None
    abritamr_subclass_present_in: Optional[list] = None
    species: Optional[list] = None
    genus: Optional[list] = None
    exception: BaseInfer  = None
    def __post_init__(self):
        if self.result not in ['Susceptible','Resistant','Intermediate','Not-reported']:
            raise ValueError(f"result must be one of 'Susceptible', 'Resistant', 'Intermediate' or 'Not-reported'")
        if not self.drugname :
            raise ValueError(f"You must supply a 'drugname'")
        if not self.abritamr_subclass_present_contains and not self.abritamr_subclass_present_equals and not self.abritamr_subclass_present_in and not self.abritamr_subclass_present_and and not self.abritamr_subclass_present_or and not self.gene_present and not self.mech_present:
            raise ValueError(f"You must supply one a abritamr_subclass, abritamr_class, gene or other mechanism")
        
        

@dataclass(frozen=True)
class BaseCriteria:
    action: Literal['reportable','not-reportable']
    abritamr_class: Optional[list] = None
    abritamr_subclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    gene: Optional[list]  = None
    # exception: Optional[list]  = None
    def __post_init__(self):
        if self.action not in ['reportable','not-reportable']:
            raise ValueError(f"action must be one of reportable or not-reportable")
        if not self.abritamr_class and not self.abritamr_subclass and not self.gene:
            raise ValueError(f"You must supply one of 'abritamr_class', 'abritamr_subclass' or 'gene'")
        
@dataclass(frozen=True)
class Criteria:
    action: Literal['reportable','not-reportable']
    abritamr_class: Optional[list] = None
    abritamr_subclass: Optional[list]  = None
    species: Optional[list]  = None
    genus: Optional[list]  = None
    gene: Optional[list]  = None
    exception: BaseCriteria  = None
    def __post_init__(self):
        if self.action not in ['reportable','not-reportable']:
            raise ValueError(f"Action must be one of reportable or not-reportable")
        if not self.abritamr_class and not self.abritamr_subclass and not self.gene:
            raise ValueError(f"You must supply one of 'abritamr_class', 'abritamr_subclass' or 'gene'")

@dataclass(frozen=True)
class AMRtype:
    amrtype: Optional[list] = None
    abritamr_subclass: Optional[list] = None
    gene: Optional[list] = None



def wrangle_lists_infer(infer:dict) -> dict:

    listkeys = [
        "gene",
        "species",
        "genus","abritamr_subclass_present_contains", 
        "abritamr_subclass_present_equals",
        "abritamr_subclass_present_and",
        "abritamr_subclass_present_or",
        "abritamr_subclass_present_in",
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
    # print(infer)
    if 'exception' in infer and infer['exception'] != []:
        print(f"Found an exception")
        exceptioninfer = []
        for cr in infer['exception']:
            ct = wrangle_lists_infer(cr)
            bif = BaseInfer(**ct)
            
            exceptioninfer.append(bif)
        infer['exception'] = exceptioninfer
    
    ifr = Infer(**infer)

    return ifr


def wrangle_lists_criteria(criteria:dict) -> dict:

    listkeys = ["gene","species","genus","abritamr_class","abritamr_subclass"]

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


def wrangle_lists_amrtype(amrtype:dict) -> dict:

    listkeys = ["gene","abritamr_subclass"]

    for l in listkeys:
        if l in amrtype:
            lst = []
            # print(l)
            if amrtype[l] != "":
                lst = amrtype[l].split(",")
                lst = [i.strip() for i in lst]
            amrtype[l] = lst
    return amrtype

def construct_amrtype(amrtype:dict) -> Criteria:
    # check_required(criteria = criteria)
    # check_one_ofs(criteria = criteria)
    amrtype = wrangle_lists_amrtype(amrtype)
    
    
    ctr = AMRtype(**amrtype)

    return ctr

def get_abritamr_configs(cfgtype:str):
    
    _funcs= {
        'infer': construct_inference,
        'amr_type': construct_amrtype,
        'reportable': construct_criteria
    }

    try:
        with open(f"{cfgtype}_criteria.json", "r") as j:

            config = json.load(j)
    except Exception as e:
        print(f"An error has occured: {e}")
        raise SystemExit(1)
    

    listcriteria = []

    for c in config:
        c1 = _funcs[cfgtype](c)
        listcriteria.append(c1)

    return listcriteria
