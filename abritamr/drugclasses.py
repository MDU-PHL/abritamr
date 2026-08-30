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
        # # print(f"{val} not found in {key} for {_type}")
        log.critical("The entry is not present in the reference catalog")
        return "unknown"

def get_amrrules_mutation(val:str, refgenes:pd.DataFrame) -> str:

    try:
        return refgenes[refgenes["abritamr_accession_key"] == val][f"amrrules_mutation"].values[0]
    except:
        # # print(f"{val} not found in abritamr_accession_key for amrrules_mutation")
        # log.critical("The entry is not present in the reference catalog")
        return "-"

def get_mechanism(val:str, refgenes:pd.DataFrame) -> str:

    try:
        return refgenes[refgenes["abritamr_accession_key"] == val][f"abritamr_mechanism"].values[0]
    except:
        # # print(f"{val} not found in abritamr_accession_key for abritamr_mechanism")
        log.critical("The entry is not present in the reference catalog")
        return "unknown"

def get_accession_key(key:str, val:str, refgenes:pd.DataFrame) -> str:
    
    try:
        return refgenes[refgenes[key] == val][f"abritamr_accession_key"].values[0]
    except:
        # # print(f"{val} not found in {key} for abritamr_accession_key")
        log.critical("The entry is not present in the reference catalog")
        return "unknown"
    
def get_classes(key:str, val:str, refgenes:pd.DataFrame) -> tuple:
    abacc = get_accession_key(key = key, val = val, refgenes = refgenes)
    _class = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_class")
    _subclass = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_subclass")
    pmid = get_class(key = key, val = val, refgenes = refgenes, _type = "pubmed_reference")
    db_version = get_class(key = key, val = val, refgenes = refgenes, _type = "db_version")
    key = get_class(key = key, val = val, refgenes = refgenes, _type = "abritamr_accession_key")
    amrrules_mut = get_amrrules_mutation(val = abacc, refgenes = refgenes)
    mech = get_mechanism(val = abacc, refgenes = refgenes)
    
    return _class,_subclass,pmid,db_version,key,amrrules_mut,mech

def find_classes(refgenes:pd.DataFrame,accession:str) -> str:

    _class = _subclass =  "unknown"
    
    for key in ['refseq_protein_accession', 'refseq_nucleotide_accession', 'genbank_protein_accession', 'genbank_nucleotide_accession']:
        if accession in refgenes[key].unique().tolist():
            # # print(f"Found {accession} in {key}")
            _class,_subclass,pmid,db_version,akey,amrrules_mut,mech = get_classes(key = key, val = accession,refgenes = refgenes)
            break

    return _class,_subclass,pmid,db_version,akey,amrrules_mut,mech

def apply_classes(amr:dict, species:str, sid:str, catalog:str) -> dict:
    # # print(amr)
    refgenes = get_refgenes(pth = catalog)
    for row in amr:
        # # print(row)
        _class,_subclass,pmid,db_version,key,amrrules_mut,mech = find_classes(refgenes = refgenes, accession = row['Closest reference accession'])

        
        row['abritamr_class'] = _class
        row['abritamr_subclass'] = _subclass
        row['species']= species
        row['sample_id'] = sid
        row['pmid'] = pmid
        row['abritamr_class_version'] = db_version
        row['abritamr_accession_key'] = key
        row['amrrules_mutation'] = amrrules_mut
        row['abritamr_mechanism'] = mech
        # # print(row)
    return amr