import pandas as pd
import numpy as np

import json
import pathlib

from filter_inferrence import construct_filter


    
def infer_resistance(amr:dict, species: str= "", genus:str="", default:str="Susceptible") -> dict:

    inferred = construct_filter(result = amr, species=species, genus=genus, default= default)


    
    return inferred


def infer(amr:dict, species: str= "", genus:str="", default:str="Susceptible", sid:str="abritamr", output:str="abritamr_inferred_ast.csv") -> bool:

    inferred = infer_resistance(amr = amd, species = species, genus= genus, default = default)

    if inferred != {}:

        