import pandas as pd
import numpy as np
import logging
import json
import pathlib
import sys

from abritamr.utils import get_refgenes

from abritamr.filter_reportable import construct_filter
from abritamr.logger import log

# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
# handler = logging.StreamHandler(sys.stderr)

# log = logging.getLogger(__name__)
# log.addHandler(handler)


def get_class(key:str, val:str, refgenes:pd.DataFrame, _type:str) -> str:

    try:
        return refgenes[refgenes[key] == val][f"{_type}"].values[0]
    except:
        log.critcal("The entry is not present in the reference catalog")
        return "unknown"

def get_classes(key:str, val:str, refgenes:pd.DataFrame) -> tuple:

    _class = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_class")
    _subclass = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_subclass")

    return _class,_subclass

def find_classes(refgenes:pd.DataFrame,accession:str) -> str:

    _class = _subclass =  "unknown"
    
    for key in ['refseq_protein_accession', 'refseq_nucleotide_accession', 'genbank_protein_accession', 'genbank_nucleotide_accession']:
        if accession in refgenes[key].unique().tolist():
            _class,_subclass = get_classes(key = key, val = accession,refgenes = refgenes)
            break

    return _class,_subclass

def add_abritamr_results(amr:dict, species: str= "", sid : str = "", cfgpath:str="") -> pd.DataFrame:
    log.info(f"Adding abritamr classes and determining gene status.")
    refgenes = get_refgenes()
    for row in amr:
        _class,_subclass = find_classes(refgenes = refgenes, accession = row['Closest reference accession'])

        
        row['abritamr_class'] = _class
        row['abritamr_subclass'] = _subclass
        row['species']= species
        row['sample_id'] = sid
        row['gene'] = row['Element symbol']
        row = construct_filter( result = row , cfgpath = cfgpath)
        
        
        
        # print(row)
    return amr



    