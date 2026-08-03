import pandas as pd
import numpy as np

import json
import pathlib

from abritamr.filter_amrtype import construct_filter

import logging

logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def get_amr_type(amr:dict, species: str= "", genus:str="", cfgpath:str="") -> pd.DataFrame:
    
    for row in amr:
        
        result_tocheck = {
            'abritamr_subclass':row['abritamr_subclass'],
            'species':species,
            'genus':genus,
            'gene':row['Element symbol']
        }
        report = construct_filter( result = result_tocheck, cfgpath = cfgpath )
        row['abritamr AMR type'] = report

        
    return amr
    
