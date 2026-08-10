from dataclasses import dataclass
from dataclasses import asdict
from cel import evaluate
from typing import Optional,Literal
import pandas as pd
import json
import pathlib


@dataclass(frozen=True)
class ClassDefintion:
    class_id:int
    class_curation_id:str
    defintion:str
    _class:str = None
    _subclass:str = None

@dataclass(frozen=True)
class Criteria:
    
    criteria_id: str
    criteria_version:str
    criteria: str
    additional_status_criteria:str = None
    additional_type_criteria: str = None
    status: str = None
    amrtype: str = None
    species: str = None
    drugname:str = None
    inferred: str = None
    

    def __post_init__(self):
        
        if not self.drugname   and self.inferred:
            raise ValueError(f"You must supply a 'drugname'")
        if not self.status and not self.inferred:
            raise ValueError(f"You must supply one of status (reportable or not-reportable) or inferred (example resistant, R, susceptible)")
      

def get_abritamr_reporting(cfgpath:str= "") -> list:
    reps = pd.read_csv(f"{cfgpath}")

    abritamr_rep = []
    for criteria in reps.to_dict(orient="records"):
        # print(criteria)
        rep = Criteria(**criteria)
        abritamr_rep.append(rep)

    return abritamr_rep

def get_abritamr_defs(cfgpath:str= "") -> list:

    defs = pd.read_csv(cfgpath)
    abritamr_defs = []
    for d in defs.to_dict(orient = "records"):
        dfs = ClassDefintion(**d)
        abritamr_defs.append(d)
    
    return abritamr_defs