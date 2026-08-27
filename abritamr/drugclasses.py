from webbrowser import get

import pandas as pd
import numpy as np
import logging
import json
import pathlib
import sys
from abritamr.logger import log
from abritamr.utils import get_refgenes


def get_class(key:str, val:str, refgenes:pd.DataFrame, _type:str) -> str:

    try:
        return refgenes[refgenes[key] == val][f"{_type}"].values[0]
    except:
        log.critcal("The entry is not present in the reference catalog")
        return "unknown"



def get_classes(key:str, val:str, refgenes:pd.DataFrame) -> tuple:

    _class = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_class")
    _subclass = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_subclass")
    pmid = get_class(key = key, val = val, refgenes = refgenes, _type = "pubmed_reference")
    db_version = get_class(key = key, val = val, refgenes = refgenes, _type = "db_version")
    key = get_class(key = key, val = val, refgenes = refgenes, _type = "abritamr_accession_key")
    amrrules_mut = get_class(key = key, val = val, refgenes = refgenes, _type = "amrrules_mutation")
    
    return _class,_subclass,pmid,db_version,key,amrrules_mut

def find_classes(refgenes:pd.DataFrame,accession:str) -> str:

    _class = _subclass =  "unknown"
    
    for key in ['refseq_protein_accession', 'refseq_nucleotide_accession', 'genbank_protein_accession', 'genbank_nucleotide_accession']:
        if accession in refgenes[key].unique().tolist():
            _class,_subclass,pmid,mech,db_version,akey,amrrules_mut = get_classes(key = key, val = accession,refgenes = refgenes)
            break

    return _class,_subclass,pmid,mech,db_version,akey,amrrules_mut

def apply_classes(amr:dict, species:str, sid:str, catalog:str) -> dict:

    refgenes = get_refgenes(pth = catalog)
    for row in amr:
        _class,_subclass,pmid,db_version,key,amrrules_mut = find_classes(refgenes = refgenes, accession = row['Closest reference accession'])

        
        row['abritamr_class'] = _class
        row['abritamr_subclass'] = _subclass
        row['species']= species
        row['sample_id'] = sid
        row['pmid'] = pmid
        row['abritamr_class_version'] = db_version
        row['abritamr_accession_key'] = key
        row['amrrules_mutation'] = amrrules_mut
    
    return amr