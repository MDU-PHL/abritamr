import pandas as pd
import pathlib
import datetime
import numpy
import json
import subprocess
from dataclasses import asdict
from abritamr.logger import log
# from abritamr.utils import get_refgenes
# from abritamr.criteria import get_abritamr_reporting,get_abritamr_defs
from cel import evaluate
import sys


def get_species_list() -> list:

    try:
        with open(f"{pathlib.Path(__file__).parent / 'configs' / 'amr_rules_species.json'}", "r") as s:
            sp = json.load(s)
        

        return sp
    except Exception as e:
        print(f"An error occured identifying supported species.")
        # raise SystemExit

def url_stub() ->str:

    return "https://raw.githubusercontent.com/AMRverse/AMRrules/refs/heads/main/rules/"


def get_rules(species:str ) -> dict:

    sfx = "_".join(species.split(" "))
    url = url_stub()
    p = subprocess.run(f"wget {url}{sfx}.tsv", shell = True, capture_output = True, encoding = "utf-8")
    if p.returncode == 0:
        log.info(f"AMR rule downloaded for {species}.")
        return f"{sfx}.tsv"
    else:
        log.critical(f"Something went wrong with the download of rules. {p.stderr}.")
        raise SystemExit

def open_rules(pth:str) -> dict:

    try:
        df = pd.read_csv(pth, sep = "\t")
        return df.to_dict(orient = "records")
    except Exception as e:
        print(f"Something went wrong opening the rules file. The following error occured : {e}")
        raise SystemExit(1)



def _get_date():

    return datetime.datetime.today().strftime('%Y-%m-%d')

     

def parse_rule(row:dict, simple_rules:dict) -> str:
    
    if "|" not in row['gene'] and "&" not in row['gene']:
        acc = row['protein accession'] if row['protein accession'] != '-' else row['nucleotide accession']
        rl = f"(abritamr_accession_key == '{acc}')"
    else:
        rl = row['gene']
        for s in simple_rules:
            rl = rl.replace(s, f"abritamr_accession_key == '{simple_rules[s]}'")
        rl = rl.replace("&", " && ")
        rl = rl.replace("|", " || ")
    
    return rl

def get_simple_rules(rules:list) -> list:

    # rows = []
    simple_rules = {}
    for rule in rules:
        if "&" not in rule['gene'] and "|" not in rule['gene']:
            simple_rules[rule['ruleID']] = rule['gene']
    
    return simple_rules
        

def wrangle_the_rules(rules:list, simple_rules:dict, species:str, cfg:dict, evidence_grade:str = 0) -> list:

    grades = cfg['grades']
    row = []
    for rule in rules:
        ev = grade[rule['evidence grade']]
        if ev >= evidence_grade:
            dr = rule['drug'].capitalize()
            rule_id = f"AMRrules-{rule['ruleID']}"
            rule_version = f"AMRrules-downloaded-{_get_date()}"
            inferred = f"{rule['phenotype']} ({rule['clinical category']})"
            rl = parse_rule(rule)
            src = f"AMRrules" if rule['PMID'] == '-' else f"AMRrules (PMID: {rule['PMID']})"
            d = {
                'drugname':dr,
                'species': species,
                'rule_id':rule_id,
                'rule_version':rule_version,
                'rule':rl,
                'original_rule': row['gene'], # remove after testing
                'inferred':inferred,
                'source':src
                }
            row.append(d)
def generate_rules(species:str, evidence_grade:int, cfg:dict) -> list:

    rule_file = get_rules(species = species)
    rules = open_rules(pth = rule_file)
    simple_rules = get_simple_rules(rules = rules)
    results = wrangle_the_rules(rules = rules, simple_rules = simple_rules, species = species, evidence_grade = evidence_grade, cfg = cfg)

    return rows
# @dataclass(frozen=True)
# class InferRules:
#     drugname: str
#     species: str
#     rule_id:str
#     rule_version:str
#     rule:str
#     inferred:str
#     source:str = "not supplied"


