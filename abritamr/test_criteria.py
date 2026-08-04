from dataclasses import dataclass
from dataclasses import asdict
from cel import evaluate
from typing import Optional,Literal
import pandas as pd
import json
import pathlib


@dataclass(frozen=True)
class Criteria:
    
    criteria_id: str
    criteria_version:str
    criteria: str
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
      


reps = pd.read_csv(f"{pathlib.Path(__file__).parent / 'configs' / 'abritamr_reporting.csv'}")

test = pd.read_csv(f"{pathlib.Path(__file__).parent / 'configs' / 'abritamr_criteria_tests.csv'}")
abritamr_rep = []
for criteria in reps.to_dict(orient="records"):
    # print(criteria)
    rep = Criteria(**criteria)
    abritamr_rep.append(rep)

test_result = []
for row in test.to_dict(orient = "records"):

    for rep in abritamr_rep:
        res = evaluate(asdict(rep)['criteria'], row)
        if res:
            row['abritamr_AMR_reporting'] = asdict(rep)['status']
            row['criteria_id'] = asdict(rep)['criteria_id']
            row['criteria_version'] = asdict(rep)['criteria_version']
            if row['abritamr_AMR_reporting'] == row['test']:
                tr = 'pass'
            else:
                tr = 'fail'
            row['pass_fail'] = tr
    
            test_result.append(row)

    # break

t = pd.DataFrame(test_result)
print(t)
t.to_csv("test_result.csv", index = False)