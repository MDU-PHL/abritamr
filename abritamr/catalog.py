import re

import pandas as pd
import pathlib
import datetime
import numpy
import json
import re
import subprocess
from dataclasses import asdict
from abritamr.logger import log
from abritamr.utils import get_refgenes
from abritamr.criteria import get_abritamr_reporting,get_abritamr_defs
from abritamr.cel_functions import contains_any,create_cel_context,evaluate_rule
from abritamr.parse_gdstrules import get_cfg,get_rules,open_rules,_get_date
from abritamr.utils import _get_date



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
        ctx = create_cel_context(data = row, name = 'row')
        for c in class_definitions:
            cl = asdict(c)
            rs = evaluate_rule(cl['definition'], ctx)
            # rs = evaluate(c['definition'], ctx)
            if rs:
                # print(c)
                row['abritamr_class'] = cl['abritamr_class'] if cl['abritamr_class'] else _capitalise(row['abritamr_class'])
                if f"{cl['abritamr_subclass']}"  != "nan" and 'append' not in f"{cl['abritamr_subclass']}":
                    # log.info(f"criteria subclass is : {cl['abritamr_subclass']}")
                    row['abritamr_subclass'] = cl['abritamr_subclass']
                elif f"{cl['abritamr_subclass']}"  != "nan"  and 'append' in f"{cl['abritamr_subclass']}":
                    row['abritamr_subclass'] = cl['abritamr_subclass'].replace('append', row['subclass'].lower())
                else:
                    row['abritamr_subclass'] = _capitalise(row['subclass'])
                row['abritamr_class_definition'] = cl['class_curation_id']
            
            # else:
        
        # print(row)
            row['abritamr_class'] = row['abritamr_class'] if 'abritamr_class' in row else _capitalise(row['class'])
            row['abritamr_subclass']= row['abritamr_subclass'] if 'abritamr_subclass' in row else _capitalise(row['subclass'])
        # print(row)
        rows.append(row)
    return rows


def mutant_nomenclature_converter() -> dict:
    converter = {'G': 'Gly', 'A': 'Ala', 'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'C': 'Cys', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 
                 'M': 'Met', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys', 'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu', 'W': 'Trp', 
                 'Y': 'Tyr', 'F': 'Phe'}
    return converter

def sub_aa_del(ref:str, pos:int, alt:str, xtr:str) -> str:
    cvt = mutant_nomenclature_converter()
    pos1 = f"{pos}"
    pos2 = f"{pos+len(ref)-1}"
    ref1= ref[0]
    ref2 = ref[-1]
    res = f"{ref1}{pos1}_{ref2}{pos2}{alt}{xtr}"
    for c in cvt:
        res= res.replace(c, cvt[c])
        # alt = alt.replace(c, cvt[c])
    # var = f"p.{ref}{pos}del{alt}{xtr}"
    return f"p.{res}"

def sub_aa_mutations(ref:str,pos:int, alt:str, xtr:str) -> str:
    reqa = len(ref) == len(alt) # check if it is a deletion or mutli substitution
    if 'del' in alt and reqa == False:
        # alt = f"del{alt.replace('del','')}"
        return sub_aa_del(ref = ref, pos = pos, alt = alt, xtr = xtr)
        # ref = f"{ref[0]}"

    cvt = mutant_nomenclature_converter()
    if len(alt) > 1 and 'Ter' not in alt:
        # ref = ""
        pos = f"{pos}_{pos+len(alt)-1}"
    
    for c in cvt:
        ref = ref.replace(c, cvt[c])
        if 'Ter' not in alt:
            alt = alt.replace(c, cvt[c])

    var = f"p.{ref}{pos}{alt}{xtr},p.{pos}{alt}{xtr}"
    return var

def sub_nt_mutations(ref:str, pos:int, alt:str, promoter: bool = False) -> str:
    var = f"c.{pos}{ref}>{alt}" if promoter else f"c.[{pos}{ref}>{alt}]"
    return var

def parse_snps_for_amrrules(refgenes:list) -> list:
    rgx = re.compile(r'(\D+)(-?\d+)(\D+)(.*)')
    for row in refgenes:
        if "POINT" == row['subtype']:
            # print(f"Parsing SNP for {row['allele']}")
            try:
                ref,pos,alt,xtr = rgx.match(row['allele'].split("_")[-1]).groups()
                if "promoter" not in row["product_name"] and "ribosomal RNA" not in row["product_name"]:
                    row['amrrules_mutation'] = sub_aa_mutations(ref = ref, pos = int(pos), alt = alt, xtr = xtr)
                elif "promoter" in row["product_name"]:
                    row['amrrules_mutation'] = sub_nt_mutations(ref = ref, pos = int(pos), alt = alt, promoter = True)
                elif "ribosomal RNA" in row["product_name"]:
                    row['amrrules_mutation'] = sub_nt_mutations(ref = ref, pos = int(pos), alt = alt, promoter = False)
                else:
                    row['amrrules_mutation'] = row['allele'].split("_")[-1]
            except Exception as e:
                log.warning(f"Could not parse SNP for {row['allele']}. The following error occured : {e}")
                row['amrrules_mutation'] = row['allele'].split("_")[-1]
    return refgenes

def get_amrrules(output_dir:str) -> dict:
    cfg = get_cfg()
    species_list = cfg['species']
    rules_dict = {}
    for sp in species_list:
        rule_file = get_rules(species = sp, output_dir = output_dir)
        rules_dict[sp] = open_rules(pth = rule_file)

    return rules_dict

def apply_amrrule_genecontext(output_dir:str, rows:list, dflt_rs: str= '-', dflt_hg : str = 'high') -> list:

    rules_dict = get_amrrules(output_dir = output_dir)
    for row in rows:
        if row['priority_status'] != dflt_hg:
            ar = set()
            for sp,rule in rules_dict.items():
                for rl in rule:
                    acc = rl['protein accession'] if rl['protein accession'] != '-' else ""
                    if acc == row['abritamr_accession_key']:
                        row['priority_status'] = rl['gene context'] if rl['gene context'] != '-' else dflt_rs
                        ar.add(f"('{sp}' in row.species)")
                        row['criteria_id'] = f"AMRrules{rl['ruleID']}"
                        row['criteria_version'] = f"AMRrules-downloaded-{_get_date()}"
            row['additional_status_criteria'] = " | ".join(ar)
    return rows


def apply_amrtyping(amrtyping_definitions:list,refgenes:list, dflt_rs: str= '-') -> list:


    rows = []
    for row in refgenes:
        row['priority_status'] = dflt_rs
        
        ctx = create_cel_context(data = row, name = 'row')
        for cr in amrtyping_definitions:
            crt = asdict(cr)
            
            rs = evaluate_rule(crt['criteria'], ctx)
            if rs:
                
                row['priority_status'] = crt['status']
                row['amrtype'] = crt['amrtype']
                row['criteria_id'] = crt['criteria_id']
                row['criteria_version'] = crt['criteria_version']
                row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
                row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
                
                break
                
            
        rows.append(row)
    
    return rows


def wrangle_catalog(catalog:str, previous_catalog:str, amrtyping_definitions:str, class_definitions:str, output_dir:str, src:str = "abritamr") -> bool:

    refgenes = construct_refgenes_starter(catalog = catalog, previous_catalog=previous_catalog, src = src)
   
    ccdict = get_class_criteria(cfgpath = class_definitions)
    refgenes = apply_classes(class_definitions = ccdict, refgenes = refgenes)
    
    rcdict = get_abritamr_reporting(cfgpath =amrtyping_definitions)
    if rcdict == []:
        log.warning(f"No amr typing rules found in {amrtyping_definitions}. No typing will be captured and amrtyping will not be able to be undertaken using this catalogue.")

    refgenes = apply_amrtyping(amrtyping_definitions = rcdict, refgenes= refgenes)
    refgenes = apply_amrrule_genecontext(output_dir = output_dir, rows = refgenes)
    refgenes = parse_snps_for_amrrules(refgenes = refgenes)
    
    # pd.DataFrame(refgenes).to_csv(f"{output}", index = False)
    
    return refgenes