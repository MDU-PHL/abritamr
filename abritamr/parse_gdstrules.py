import re

import pandas as pd
import pathlib
import datetime
import json
import subprocess
from dataclasses import asdict
from abritamr.logger import log
from abritamr.utils import _get_date
# from abritamr.utils import get_refgenes
# from abritamr.criteria import get_abritamr_reporting,get_abritamr_defs



def get_cfg() -> list:

    try:
        with open(f"{pathlib.Path(__file__).parent / 'configs' / 'amr_rules_config.json'}", "r") as s:
            sp = json.load(s)

        return sp
    except Exception as e:
        print(f"An error occured identifying supported species.")
        # raise SystemExit

def url_stub() ->str:

    return "https://raw.githubusercontent.com/AMRverse/AMRrules/refs/heads/main/rules/"


def get_rules(species:str, output_dir:str) -> dict:

    sfx = "_".join(species.split(" "))
    url = url_stub()
    p = subprocess.run(f"wget -nc -O {output_dir}/amr_rules_{sfx}.tsv {url}{sfx}.tsv", shell = True, capture_output = True, encoding = "utf-8")
    if pathlib.Path(f"{output_dir}/amr_rules_{sfx}.tsv").exists():
        log.info(f"AMR rule downloaded for {species}.")
        return f"{output_dir}/amr_rules_{sfx}.tsv"
    else:
        log.critical(f"Something went wrong with the download of rules. {p.stderr}.")
        raise SystemExit
    
def open_rules(pth:str) -> dict:

    try:
        df = pd.read_csv(pth, sep = "\t")
        return df.fillna("").to_dict(orient = "records")
    except Exception as e:
        print(f"Something went wrong opening the rules file. The following error occured : {e}")
        raise SystemExit(1)



def get_accession_key(row:dict) -> str:
    acc = row['protein accession'] if row['protein accession'] != '-' else row['nucleotide accession']
    acc = acc.split(":")[0] if ":" in acc else acc
    return acc

def special_rule_for_amrrules_rna(mutation : str) -> str:
    if mutation.endswith("]"):
        return f"{mutation.split(']')[0]}]"
    else:
        return mutation

def parse_rule(row:dict, simple_rules:dict) -> str:

    acc = get_accession_key(row = row)
    
    mt = f"contains_any(row.amrrules_mutation, '{special_rule_for_amrrules_rna(row['mutation'])}') && " if row['mutation'] != "-" else ""
    if "|" not in row['gene'] and "&" not in row['gene'] and acc != "-":
        return f"({mt}  contains_any(row.abritamr_accession_key,'{acc}'))"
    elif "|" not in row['gene'] and "&" not in row['gene'] and acc == "-":
        return "-"
    else:
        rl = row['gene']
        for s in simple_rules:
            # print(f"Replacing {s} with {simple_rules[s]}")
            rl = rl.replace(s, f"contains_any(row.abritamr_accession_key,'{simple_rules[s]}')")
        rl = rl.replace("&", " && ")
        rl = rl.replace("|", " || ")
        rl = f"{rl} && {mt}" if mt != "" else rl
        return rl


def get_simple_rules(rules:list) -> list:

    # rows = []
    simple_rules = {}
    for rule in rules:
        if "&" not in rule['gene'] and "|" not in rule['gene']:
            acc = get_accession_key(row = rule)
            if acc != "-":
                simple_rules[rule['ruleID']] = acc
            
    
    return simple_rules
        

def wrangle_the_rules(rules:list, simple_rules:dict, species:str, cfg:dict, evidence_grade:str = 'very low') -> list:

    grades = cfg['grades']
    ccl = cfg['clinical_category']
    row = []
    for rule in rules:
        try:
            ev = int(grades[rule['evidence grade'].strip()])
            yev = int(grades[evidence_grade.strip()])
            if ev >= yev:
                dr = rule['drug'].capitalize() if rule['drug'] != '-' else rule['drug class'].capitalize()
                rule_id = f"AMRrules-{rule['ruleID']}"
                rule_version = f"AMRrules-downloaded-{_get_date()}"
                cc = ccl[rule['clinical category']] if rule['clinical category'] in ccl else "Unknown"
                infrd = f"{cc} ({rule['phenotype']})"
                rl = parse_rule(rule, simple_rules = simple_rules)
                src = f"AMRrules" if rule['PMID'] == '-' else f"AMRrules (PMID: {rule['PMID']})"
                if rl != "-":
                    d = {
                        'drugname':dr,
                        'species': species,
                        'rule_id':rule_id,
                        'rule_version':rule_version,
                        'rule':rl,
                        'original_rule': rule['gene'], # remove after testing
                        'inferred':infrd,
                        'source':src
                        }
                    row.append(d)
        except Exception as e:
            log.warning(f"Something went wrong parsing the rule {rule}. The following error was reported : {e}")
            continue
    return row



def generate_rules(species:str, evidence_grade:int, cfg:dict, output_dir:str) -> list:

    rule_file = get_rules(species = species, output_dir = output_dir)
    rules = open_rules(pth = rule_file)
    simple_rules = get_simple_rules(rules = rules)
    # print(simple_rules)
    results = wrangle_the_rules(rules = rules, simple_rules = simple_rules, species = species, evidence_grade = evidence_grade, cfg = cfg)
    # print(results)
    return results

def get_amrrules_for_species(evidence_grade:int,output_dir:str, species:str="all", rules_dict:dict={}) -> list:
    cfg = get_cfg()
    species_list = cfg['species']
    
    if species == "all":
        
        for sp in species_list:
            # print(sp)
            rules = generate_rules(species = sp, evidence_grade = evidence_grade, cfg = cfg, output_dir = output_dir)
            rules_dict[sp] = rules
    else:
        if species not in species_list:
            log.critical(f"Species {species} is not supported by AMRverse rules. Please check your input.")
            raise SystemExit(1)
        else:
            rules = generate_rules(species = species, evidence_grade = evidence_grade, cfg = cfg, output_dir = output_dir)
            rules_dict[species] = rules

    return rules_dict


def get_additional_rules(pth:str) -> list:
    log.info(f"Opening additional rules file {pth}")
    try:
        rules = pd.read_csv(pth)
        return rules.fillna("").to_dict(orient = "records")
    except Exception as e:
        log.critical(f"Something went wrong opening the additional rules file. The following error occured : {e}")
        raise SystemExit(1)
    

def add_rules_to_existing(rules:dict, additional_rules:str) -> list:
    additional_rules = get_additional_rules(pth = additional_rules)
    for ar in additional_rules:
        sp = ar['species']
        if sp in rules:
            rules[sp].append(ar)
        else:
            rules[sp] = [ar]

    return rules
# @dataclass(frozen=True)
# class InferRules:
#     drugname: str
#     species: str
#     rule_id:str
#     rule_version:str
#     rule:str
#     inferred:str
#     source:str = "not supplied"


