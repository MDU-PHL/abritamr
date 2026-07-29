import pandas as pd
import numpy as np

import json
import pathlib

from filter_reportable import construct_filter

def open_amrfinder(pth : pathlib.Path) -> pd.DataFrame:

    if pth.exists():
        df =  pd.read_csv(pth, sep = "\t")
        return df.to_dict(orient='records')
    else:
        print(f"The input file provided {pth} does not exist or is inaccessible. Please try again.")
        raise SystemExit(1)

def get_refgenes()-> pd.DataFrame:

    rf = pathlib.Path(__file__).parent / 'db' /'refgenes_latest.csv'
    if rf.exists():
        return pd.read_csv(rf)
    else:
        print(f"Something is very wrong - the reference catalog is not present.")
        raise SystemExit(1)

def get_class(key:str, val:str, refgenes:pd.DataFrame, _type:str) -> str:

    try:
        return refgenes[refgenes[key] == val][f"{_type}"].values[0]
    except:
        print("The entry is not present in the reference catalog")
        return "unknown"

def get_classes(key:str, val:str, refgenes:pd.DataFrame) -> tuple:

    _class = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_class")
    _subclass = get_class(key = key, val= val, refgenes = refgenes, _type = "abritamr_subclass")

    return _class,_subclass

def find_classes(refgenes:pd.DataFrame,accession:str) -> str:

    _class = _subclass =  "unknown"
    
    for key in ['refseq_protein_accession', 'refseq_nucleotide_accession', 'genbank_protein_accession', 'genbank_nucleotide_accession']:
        # print(key)
        # print(accession)
        if accession in refgenes[key].unique().tolist():
            # print(refgenes[key])
            _class,_subclass = get_classes(key = key, val = accession,refgenes = refgenes)
            break

    return _class,_subclass

def infer_resistance():
    pass
    
def add_abritamr_results(refgenes:pd.DataFrame, amr:dict, species: str= "", genus:str="") -> pd.DataFrame:

    for row in amr:
        _class,_subclass = find_classes(refgenes = refgenes, accession = row['Accession of closest sequence'])
        result_tocheck = {
            'abritamr_class':_class,
            'abritamr_subclass':_subclass,
            'species':species,
            'genus':genus,
            'gene':row['Gene symbol']
        }
        report = construct_filter( result = result_tocheck )
        row['abritamr_class'] = _class
        row['abritamr_subclass'] = _subclass
        row['abritamr_AMR_reporting'] = report
        row['species'] = species
        row['genus'] = genus

        # print(row)
    
    # result = pd.DataFrame(amr)
    
    return amr
    