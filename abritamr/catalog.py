import pandas as pd
import pathlib
import datetime
import json
from dataclasses import asdict
from logger import log
from utils import get_refgenes
from criteria import get_abritamr_reporting,get_abritamr_defs
from cel import evaluate
import sys


def _get_date():

    return datetime.datetime.today().strftime('%Y-%m-%d')


def _get_new_catalog() -> str:
    logger.info(f"Getting reference catalog from ncbi.")
    updated_html = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"
                    
    p = subprocess.run(f"wget -O ReferenceGeneCatalog.txt {updated_html}", shell = True, capture_output = True, encoding = "utf-8")
    if p.returncode == 0:
        log.info(f"Reference catalog downloaded.")
        return "ReferenceGeneCatalog.txt"
    else:
        logger.critical(f"Something went wrong with the download of reference catalog. {p.stderr}.")
        raise SystemExit

def _make_key(df):
    if 'key' not in list(df.columns):
        df['mut_acc'] = df[['allele','whitelisted_taxa','refseq_nucleotide_accession']].apply( lambda x: '_'.join(x), axis = 1)
        df['key'] = numpy.where(df['refseq_protein_accession'] == '',df['genbank_protein_accession'],df['refseq_protein_accession'])
        df['key'] = numpy.where(df['subtype'] == 'POINT', df['mut_acc'], df['key'])

def catalog_df(catalog:str,src:str = 'abritamr') -> pd.DataFrame: # will need to generalise for amrrules and others??

    new_catalog = catalog
    if catalog == "":
        new_catalog = _get_new_catalog()
    
    refs = pd.read_csv(new_catalog, sep = "\t")
    refs = refs.fillna("")
    if src == "abritamr":
        refs['gene'] = refs['allele'].fillna(refgenes['gene_family'])
        refs = _make_key(df = refs)
    
    # refs = refs.to_dict(orient = "records")
    
    return refs

# abritamr specific
def _capitalise(x):

    a = x.split('/')
    a = [i.capitalize() for i in a]
    _list = list(map(lambda x: x.replace('Carbapenem','Carbapenemase'), a))
 
    return '/'.join(_list)
# abritamr specific
def _get_rename()-> dict:

    cfg = pathlib.Path(__file__).parent / "configs" / "update_vars.json"
    if cfg.exists():
        with open(f"{cfg}", "r") as j:
            return json.load(j)

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
    # print(previous_catalog)
    previous_catalog = previous_catalog.rename(columns = {"class_new": "class_prev", "subclass_new":"subclass_prev"})
    new_catalog = new_catalog.rename(columns = {"class":"class_new", "subclass": "subclass_new"})
    _tmp= new_catalog.merge(previous_catalog, on = ['key'])
    # print(_tmp.columns)
    _tmp['changed'] = numpy.where((_tmp['class_prev'] != _tmp['class_new']) , 'updated', '')
    _tmp['changed'] = numpy.where((_tmp['subclass_prev'] != _tmp['subclass_new']) , 'updated', _tmp['changed'])
    
    _tmp['Previous_class'] = numpy.where(_tmp['changed'] == 'updated', _tmp['class_prev'], '')
    _tmp['Previous_subclass'] = numpy.where(_tmp['changed'] == 'updated', _tmp['subclass_prev'],'')
    changed = list(_tmp[_tmp['changed']!='']['key'])
    
    new_catalog['Status'] = numpy.where(new_catalog['key'].isin(changed), 'updated',new_catalog['Status'])
    new_catalog = new_catalog.merge(_tmp[['key','Previous_class','Previous_subclass']], on = ['key'], how = 'left')
    return new_catalog


def _new_entries(new_catalog, previous_catalog) -> pd.DataFrame:
    
    logger.info(f"Checking for new entries.")
    new_catalog['Status'] = numpy.where(new_catalog['key'].isin(list(previous_catalog['key'])), 'existing','new')
    logger.info(f"Checking for updated entries.")
    new_catalog = _updated_entries(new_catalog=new_catalog,previous_catalog=previous_catalog)
    return new_catalog

def _update_status(new_catalog, previous_catalog) -> pd.DataFrame:

    new_catalog = _new_entries(new_catalog=new_catalog,previous_catalog=previous_catalog)

    return new_catalog

def _compare_to_existing(new_catalog, previous_catalog) -> pd.DataFrame:
    try:
        if isinstance(previous_catalog,pandas.DataFrame):
            new_catalog = _update_status(new_catalog=tab,previous_catalog=previous_catalog)
            # return new_catalog
        else:
            logger.info(f"There is in no previous file to compare with.")
    except Exception as e:
        log.warning(f"Cannot compare to previous catalog : {e}")

    return new_catalog
    

def construct_refgenes_starter(catalog:str, previous_catalog:str src:str= "abritamr") -> list:

    refs = catalog_df(catalog = catalog)
    crnt = get_current(pth = previous_catalog)

    refs = _compare_to_existing(new_catalog = refs, previous_catalog = crnt)

    refgenes = refs.to_dict(orient= "records")


    return refgenes


def _rename(rename_key, row):

    _cls = _capitalise(row['class']) if row['class'] not in ['FLUOROQUINOLONE',"NITROFURAN"] else rename_key[row['class']]
    sbcls = rename_key[row['subclass']] if row['subclass'] in rename_key else _capitalise(row['subclass'])

    return _cls,sbcls


def apply_classes(class_defintitions: list, refgenes : list, src: str = 'abritamr') -> list:
    rename = _get_rename()
    rows = []
    for row in refgenes:
        # abritamr specific logic around renaming - may not be used in other future formats
        if src == 'abritamr':
            for k in rename:
                if k in row['class']:
                    row['abritamr_class'], row['abritamr_subclass'] = _rename(rename,row)

        for c in class_defintitions:
            crt = asdict(c)
            rs = evaluate(crt['definition'], row)
            if rs:
                row['abritamr_class'] = crt['class'] if crt['class'] else capitalize(row['class'])
                if crt['subclass'] and 'append' not in crt['subclass']:
                    row['abritamr_subclass'] = crt['subclass']
                elif crt['subclass'] and 'append' in crt['subclass']:
                    row['abritamr_subclass'] = crt['subclass'].replace('append', row['subclass'].lower())
                else:
                    row['abritamr_subclass'] = capitalize(row['subclass'])
                row['abritamr_class_definition'] = crt['class_curation_id']
            
        
        
        row['abritamr_class'] = capitalize(row['abritamr_class']) if 'abritamr_class' in row else capitalize(row['class'])
        row['abritamr_subclass']= capitalize(row['abritamr_subclass']) if 'abritamr_subclass' in row else capitalize(row['subclass'])
        rows.append(row)
    return rows

def apply_amrtyping(amrtyping_definitions:list,refgenes:list, dflt_rs: str= 'not-reportable') -> list:


    rows = []
    for row in refgenes:
        row['reportable_status'] = dflt_rs
        for cr in rcdict:
            crt = asdict(cr)
            rs = evaluate(crt['criteria'], row)
            if rs:
                # print(row)
                # print(crt['status'])
                row['reportable_status'] = crt['status']
                row['amrtype'] = crt['amrtype']
                row['criteria_id'] = crt['criteria_id']
                row['criteria_version'] = crt['criteria_version']
                row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
                row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
                row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
                row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
                break
                
            # else:
            #     row['reportable_status'] = dflt_rs
        rows.append(row)
    
    return rows

def wrangle_catalog(catalog:str, previous_catalog:str, amrtyping_definitions:str, class_defintitions:str output:str, src:str = "abritamr") -> bool:

    refgenes = construct_refgenes_starter(catalog = catalog, previous_catalog=previous_catalog, src = src)

    ccdict = get_class_criteria(cfgpath = class_defintitions)
    refgenes = apply_classes(class_defintitions = ccdict, refgenes = refgenes)

    rcdict = get_abritamr_reporting(cfgpath =amrtyping_definitions)
    refgenes = apply_amrtyping(amrtyping_definitions = rcdict, refgenes= refgenes)

    pd.DataFrame(refgenes).to_csv(f"{output}", index = False)
    
    return True