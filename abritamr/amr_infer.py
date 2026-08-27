 # for now - will need to be cleverer/cleaner/better 
from dataclasses import dataclass, asdict
import json

from cel_functions import create_cel_context, evaluate_rule
from criteria import InferRules
import pandas as pd
import pathlib
# print(criterias)

def priority_gdst() -> dict:

    return {
        "S":0,
        "I":1,
        "R":2,
    }

def combine_results( result: list) -> dict:
    res = []
    to_test = {
        'abritamr_class':[],
        'abritamr_subclass':[],
        'abritamr_accession_key':[],
        'amrrules_mutation':[],
        'abritamr_mechanism':[],
        # 'gene':[],
       
    }
    # print(result)
    for row in result:
        for col in to_test:
            if col in row:
                to_test[col].append(row[col])
        
    
    # for k in to_test:
    #     # print(to_test[k])
    #     val = ','.join(to_test[k]) if to_test[k] != [] else "None"
        
    #     res[k] = val

    res.append(to_test)
    
    return res

def find_rules(species:str, reference_folder:str) -> dict:

    
    with open(f"{pathlib.Path(__file__).parent/ 'configs'/ 'amr_rules_config.json'}", "r") as s:
        sp = json.load(s)

    spcs = sp['species']
    # print(spcs)
    rules = []
    for s in spcs:
        if species in s or s in species:
            print(f"Found rules for {species} in {reference_folder}/02_abritamr_{s.replace(' ', '_')}_rules.csv")
            # try:
            ruleset = pd.read_csv(f"{reference_folder}/02_abritamr_{s.replace(' ', '_')}_rules.csv")
            print(ruleset)
            rules.append(ruleset.fillna(""))
            # except Exception as e:
            #     print(f"Could not find rules for {species} in {reference_folder}/02_abritamr_{s}_rules.csv. The following error was reported : {e}. Please check the rules file exists and is formatted correctly.")
    # print(rules)
    if rules == []:
        print(f"Could not find rules for {species} in {reference_folder}. Please check the rules file exists and is formatted correctly.")
        return {}
    else:
        rules = pd.concat(rules)
        return rules.to_dict(orient = "records")

def filter_results(results:pd.DataFrame, min_cov:float = 0.9, min_id:float = 0.9) -> pd.DataFrame:
    results = results[(results['% Coverage of reference'] >= min_cov) & (results['% Identity to reference'] >= min_id)]
    return results

def create_rules(species:str, reference_folder:str) -> list:
    rules = find_rules(species = species, reference_folder = reference_folder)
    try:
        rules = [InferRules(**r) for r in rules]
    except Exception as e:
        print(f"No rules available for {species} in {reference_folder}: {e}")
        rules = []
    return rules

def gdst(results:pd.DataFrame, species:str, reference_folder:str, dflt_result:str = "Susceptible (default)") -> list:
    # results = pd.read_csv(results)
    # print(results)
    sid = results.iloc[0].get("sample_id", "unknown")
    resultsmooshed = combine_results(result = results.to_dict(orient = "records"))
    # print(resultsmooshed)
    rules = create_rules(species = species, reference_folder = reference_folder)
    # print(rules)
    if rules == []:
        print(f"No gDST will be provided for {sid}. There are no rules available for {species}.")
        return []
    gdst_final = []
    for row in  resultsmooshed:
        gdst_results = {'sample_id': sid, 'species': species}
        ctx = create_cel_context(data = row, name = 'row')
        rlt = {}
        gdst_report = []
        for rule in rules:
            if rule.drugname not in rlt:
                rlt[rule.drugname] = {'drugname': rule.drugname, 'mechanisms': [], 'inferred': [], 'rule_id': [], 'rule_version': [], 'source': []}
            # print(f"Evaluating rule {rule.rule_id} for {rule.drugname} with rule {rule.rule}")
            if evaluate_rule(rule = rule.rule, ctx = ctx):
                print(f"Rule {rule.rule_id} evaluated to True for {rule.drugname}. Inferred value: {rule.inferred}")
                mechs = []
                for key in row:
                    if key in rule.rule:
                        print(f"Key {key} triggered {rule.rule}")
                        vals = results[key].unique().tolist()
                        keys = [i for i in vals if i in rule.rule]
                        print(f"Keys {keys} triggered {rule.rule}")
                        mech = results[results[key].isin(keys)]['abritamr_mechanism'].unique().tolist()
                        mechs.extend(mech)
                rlt[rule.drugname]['mechanisms'].extend(mechs)
                rlt[rule.drugname]['inferred'].append(rule.inferred)
                rlt[rule.drugname]['rule_id'].append(rule.rule_id)
                rlt[rule.drugname]['rule_version'].append(rule.rule_version)
                rlt[rule.drugname]['source'].append(rule.source)
        for drug in rlt:
            if rlt[drug]['inferred'] == []:
                rlt[drug]['inferred'] = [dflt_result]
            else:       
                # print(f"Prioritizing {rlt[drug]['inferred']} for {drug}")
                rlt[drug]['inferred'] = [sorted(rlt[drug]['inferred'], key = lambda x: priority_gdst()[x[0].upper()] if x[0].upper() in priority_gdst() else -1, reverse = True)[0]]
            # print(f"Prioritized {rlt[drug]['inferred']} for {drug}")
            for key in ['mechanisms', 'rule_id', 'rule_version', 'source', 'inferred']:
                # print(f"Joining {key} for {drug}: {rlt[drug][key]}")
                rs = ','.join(rlt[drug][key]) if rlt[drug][key] != [] else "-"
                # print(f"Joined {key} for {drug}: {rs}")
                rlt[drug][key] = rs
            # print(rlt[drug])
                
        gdst_results = {'sample_id': sid, 'species': species, 'results': list(rlt.values())}
        gdst_final.append(gdst_results)

    return gdst_final

def gdst_results_to_df_wide(gdst_results: list) -> pd.DataFrame:
    rows = []
    for res in gdst_results:    
        row = {
            'sample_id': res['sample_id'],
            'species': res['species'],
        }
        for drug_res in res['results']:
            row[f"{drug_res['drugname']}_gDST"] = drug_res['inferred']
            row[f"{drug_res['drugname']}_mechanisms"] = drug_res['mechanisms']
            row[f"{drug_res['drugname']}_rule_id"] = drug_res['rule_id']
            row[f"{drug_res['drugname']}_rule_version"] = drug_res['rule_version']
            row[f"{drug_res['drugname']}_source"] = drug_res['source']
        rows.append(row)
    dr_cols = [f"{drug}_{suffix}" for drug in set(drug_res['drugname'] for res in gdst_results for drug_res in res['results']) for suffix in ['gDST', 'mechanisms', 'rule_id', 'rule_version', 'source']]
    cols_wide = ['sample_id', 'species'] + sorted(dr_cols)
    # cols_wide = ['sample_id', 'species'] + [f"{drug}_{suffix}" for drug in set(drug_res['drugname'] for res in gdst_results for drug_res in res['results']) for suffix in ['gDST', 'mechanisms', 'rule_id', 'rule_version', 'source']]
    return pd.DataFrame(rows)[cols_wide].sort_values(by = ['sample_id'])

def gdst_results_to_df_long(gdst_results: list) -> pd.DataFrame:
    rows = []
    for res in gdst_results:    
        for drug_res in res['results']:
            row = {
                'sample_id': res['sample_id'],
                'species': res['species'],
                'drugname': drug_res['drugname'],
                'gDST': drug_res['inferred'],
                'rule_id': drug_res['rule_id'],
                'rule_version': drug_res['rule_version'],
                'source': drug_res['source'],
                'mechanisms': drug_res['mechanisms']
            }
            rows.append(row)
    # for row in rows:
    #     print(row)
    cols_long = ['sample_id', 'species', 'drugname', 'mechanisms', 'gDST', 'rule_id', 'rule_version', 'source']

    return pd.DataFrame(rows)[cols_long].sort_values(by = ['drugname'])

    # print(gdst_results)

