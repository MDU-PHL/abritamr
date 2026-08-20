import pandas as pd
import pathlib
import datetime
import numpy
import json
import subprocess
from dataclasses import asdict
from abritamr.logger import log
from abritamr.utils import get_refgenes
from abritamr.criteria import get_abritamr_reporting,get_abritamr_defs
from abritamr.cel_functions import contains_any,create_cel_context,evaluate_rule
# from cel import evaluate
import sys


def _get_date():

    return datetime.datetime.today().strftime('%Y-%m-%d')


def _get_new_catalog() -> str:
    log.info(f"Getting reference catalog from ncbi.")
    updated_html = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"
                    
    p = subprocess.run(f"wget -O ReferenceGeneCatalog.txt {updated_html}", shell = True, capture_output = True, encoding = "utf-8")
    if p.returncode == 0:
        log.info(f"Reference catalog downloaded.")
        return "ReferenceGeneCatalog.txt"
    else:
        log.critical(f"Something went wrong with the download of reference catalog. {p.stderr}.")
        raise SystemExit

def _make_key(df):
    if 'abritamr_accession_key' not in list(df.columns):
        df['mut_acc'] = df[['allele','whitelisted_taxa','refseq_nucleotide_accession']].apply( lambda x: '_'.join(x), axis = 1)
        df['abritamr_accession_key'] = numpy.where(df['refseq_protein_accession'] == '',df['genbank_protein_accession'],df['refseq_protein_accession'])
        df['abritamr_accession_key'] = numpy.where(df['subtype'] == 'POINT', df['mut_acc'], df['abritamr_accession_key'])
    return df

def catalog_df(catalog:str,src:str = 'abritamr') -> pd.DataFrame: # will need to generalise for amrrules and others??

    new_catalog = catalog
    if catalog == "":
        new_catalog = _get_new_catalog()
        
    refs = pd.read_csv(new_catalog, sep = "\t")
   
    refs['gene'] = refs['allele'].fillna(refs['gene_family'])
    refs = refs.fillna("")
    refs = _make_key(df = refs)
    
    
    return refs

def _capitalise(x):

    a = x.split('/')
    a = [i.capitalize() for i in a]
    _list = list(map(lambda x: x.replace('Carbapenem','Carbapenemase'), a))
 
    return '/'.join(_list)

def get_current(pth:str) -> pd.DataFrame:

    current = get_refgenes(pth)
    current = current.fillna("")
    return current


def get_reportable_criteria(cfgpath: str) -> dict:

    rcdict = get_abritamr_reporting(cfgpath = cfgpath)

    return dict

def get_class_criteria(cfgpath:str) -> dict:

    ccdict = get_abritamr_defs(cfgpath = cfgpath)

    return ccdict


def _updated_entries(new_catalog, previous_catalog) -> pd.DataFrame:
    previous_catalog = previous_catalog.rename(columns = {"class_new": "class_prev", "subclass_new":"subclass_prev"})
    _tmp= new_catalog.merge(previous_catalog, on = ['abritamr_accession_key'])
    _tmp['changed'] = numpy.where((_tmp['class_prev'] != _tmp['class']) , 'updated', '')
    _tmp['changed'] = numpy.where((_tmp['subclass_prev'] != _tmp['subclass']) , 'updated', _tmp['changed'])
    
    _tmp['Previous_class'] = numpy.where(_tmp['changed'] == 'updated', _tmp['class_prev'], '')
    _tmp['Previous_subclass'] = numpy.where(_tmp['changed'] == 'updated', _tmp['subclass_prev'],'')
    changed = list(_tmp[_tmp['changed']!='']['key'])
    
    new_catalog['Status'] = numpy.where(new_catalog['abritamr_accession_key'].isin(changed), 'updated',new_catalog['Status'])
    new_catalog = new_catalog.merge(_tmp[['abritamr_accession_key','Previous_class','Previous_subclass']], on = ['abritamr_accession_key'], how = 'left')
    return new_catalog


def _new_entries(new_catalog, previous_catalog) -> pd.DataFrame:
    
    log.info(f"Checking for new entries.")
    new_catalog['Status'] = numpy.where(new_catalog['abritamr_accession_key'].isin(list(previous_catalog['abritamr_accession_key'])), 'existing','new')
    log.info(f"Checking for updated entries.")
    new_catalog = _updated_entries(new_catalog=new_catalog,previous_catalog=previous_catalog)
    return new_catalog

def _update_status(new_catalog, previous_catalog) -> pd.DataFrame:

    new_catalog = _new_entries(new_catalog=new_catalog,previous_catalog=previous_catalog)

    return new_catalog

def _compare_to_existing(new_catalog, previous_catalog) -> pd.DataFrame:
    try:
        if isinstance(previous_catalog,pd.DataFrame):
            new_catalog = _update_status(new_catalog=new_catalog,previous_catalog=previous_catalog)
        
            return new_catalog
        else:
            log.info(f"There is in no previous file to compare with.")
    except Exception as e:
        log.warning(f"Cannot compare to previous catalog : {e}")

    return new_catalog
    

def construct_refgenes_starter(catalog:str, previous_catalog:str, src:str= "abritamr") -> list:

    refs = catalog_df(catalog = catalog)
    crnt = get_current(pth = previous_catalog)
    refs = _compare_to_existing(new_catalog = refs, previous_catalog = crnt)
    refs = refs.fillna("")
    refgenes = refs.to_dict(orient= "records")

    return refgenes


def apply_classes(class_definitions: list, refgenes : list, src: str = 'abritamr') -> list:
    # rename = _get_rename()
    rows = []
    for row in refgenes:
        ctx = create_context(data = row, name = 'row')
        for c in class_definitions:
            rs = evaluate_rule(c['definition'], ctx)
            # rs = evaluate(c['definition'], ctx)
            if rs:
                # print(c)
                row['abritamr_class'] = c['abritamr_class'] if c['abritamr_class'] else _capitalise(row['abritamr_class'])
                if f"{c['abritamr_subclass']}"  != "nan" and 'append' not in f"{c['abritamr_subclass']}":
                    # log.info(f"criteria subclass is : {c['abritamr_subclass']}")
                    row['abritamr_subclass'] = c['abritamr_subclass']
                elif f"{c['abritamr_subclass']}"  != "nan"  and 'append' in f"{c['abritamr_subclass']}":
                    row['abritamr_subclass'] = c['abritamr_subclass'].replace('append', row['subclass'].lower())
                else:
                    row['abritamr_subclass'] = _capitalise(row['subclass'])
                row['abritamr_class_definition'] = c['class_curation_id']
            
            # else:
        
        # print(row)
            row['abritamr_class'] = row['abritamr_class'] if 'abritamr_class' in row else _capitalise(row['class'])
            row['abritamr_subclass']= row['abritamr_subclass'] if 'abritamr_subclass' in row else _capitalise(row['subclass'])
        # print(row)
        rows.append(row)
    return rows

def apply_amrtyping(amrtyping_definitions:list,refgenes:list, dflt_rs: str= 'not-reportable') -> list:


    rows = []
    for row in refgenes:
        row['reportable_status'] = dflt_rs
        
        ctx = create_context(data = row, name = 'row')
        for cr in amrtyping_definitions:
            crt = asdict(cr)
            
            rs = evaluate_rule(crt['criteria'], ctx)
            if rs:
                
                row['reportable_status'] = crt['status']
                row['amrtype'] = crt['amrtype']
                row['criteria_id'] = crt['criteria_id']
                row['criteria_version'] = crt['criteria_version']
                row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
                row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
                row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
                row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
                break
                
            
        rows.append(row)
    
    return rows


def wrangle_catalog(catalog:str, previous_catalog:str, amrtyping_definitions:str, class_definitions:str, output:str, src:str = "abritamr") -> bool:

    refgenes = construct_refgenes_starter(catalog = catalog, previous_catalog=previous_catalog, src = src)
   
    ccdict = get_class_criteria(cfgpath = class_definitions)
    refgenes = apply_classes(class_definitions = ccdict, refgenes = refgenes)
    
    rcdict = get_abritamr_reporting(cfgpath =amrtyping_definitions)
    if rcdict == []:
        log.warning(f"No amr typing rules found in {amrtyping_definitions}. No typing will be captured and amrtyping will not be able to be undertaken using this catalogue.")

    refgenes = apply_amrtyping(amrtyping_definitions = rcdict, refgenes= refgenes)

    
    # pd.DataFrame(refgenes).to_csv(f"{output}", index = False)
    
    return refgenes