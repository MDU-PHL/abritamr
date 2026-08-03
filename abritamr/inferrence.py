import pandas as pd
import numpy as np

import json
import pathlib

from abritamr.filter_inferrence import construct_filter


    
def infer_resistance(amr:dict, species: str= "", genus:str="", default:str="Susceptible") -> dict:
    
    inferred = construct_filter(result = amr, species=species, genus=genus, default= default)

    return inferred


def infer(amr:dict, species: str= "", genus:str="", default:str="Susceptible") ->  dict:

    for a in amr:
        a['gene'] = a["Element symbol"]

    inferred = infer_resistance(amr = amr, species = species, genus= genus, default = default)
    return inferred
    # if inferred != {}:

