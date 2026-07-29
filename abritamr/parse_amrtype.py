import pandas as pd
import numpy as np

import json
import pathlib

from filter_amrtype import construct_filter

def get_amr_type(amr:dict, species: str= "", genus:str="") -> pd.DataFrame:

    for row in amr:
        
        result_tocheck = {
            'abritamr_subclass':_subclass,
            'species':species,
            'genus':genus,
            'gene':row['Gene symbol']
        }
        report = construct_filter( result = result_tocheck )
        row['abritamr AMRtype'] = report

        # print(row)
    
    # result = pd.DataFrame(amr)
    
    return amr
    
