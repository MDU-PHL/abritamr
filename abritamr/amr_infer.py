"""Infer genotypic drug susceptibility from detected AMR mechanisms."""

import json

from abritamr.cel_functions import create_cel_context, evaluate_rule
from abritamr.criteria import InferRules
import pandas as pd
from abritamr.logger import log
import pathlib


def priority_gdst() -> dict:
    """Return the ordering used to prioritize inferred susceptibility results."""
    return {
        "S": 0,
        "I": 1,
        "R": 2,
    }


def combine_results(result: list) -> dict:
    """Collect rule-relevant fields from AMR result records."""
    to_test = {
        "abritamr_class": [],
        "abritamr_subclass": [],
        "abritamr_accession_key": [],
        "amrrules_mutation": [],
        "abritamr_mechanism": [],
    }
    for row in result:
        for col in to_test:
            if col in row:
                to_test[col].append(row[col])

    return to_test


def find_rules(species: str, reference_folder: str) -> dict:
    """Load the rule records configured for a species."""
    with open(
        f"{pathlib.Path(__file__).parent / 'configs' / 'amr_rules_config.json'}", "r"
    ) as s:
        sp = json.load(s)

    spcs = sp["species"]
    rules = []
    for s in spcs:
        if species in s or s in species:
            log.info(
                f"Found rules for {species} in {reference_folder}/02_abritamr_{s.replace(' ', '_')}_rules.csv"
            )
            # try:
            ruleset = pd.read_csv(
                f"{reference_folder}/02_abritamr_{s.replace(' ', '_')}_rules.csv"
            )
            rules.append(ruleset.fillna(""))
    if rules == []:
        log.warning(
            f"Could not find rules for {species} in {reference_folder}. Please check the rules file exists and is formatted correctly."
        )
        return {}
    else:
        rules = pd.concat(rules)
        return rules.to_dict(orient="records")


def filter_results(
    results: pd.DataFrame, min_cov: float = 0.9, min_id: float = 0.9
) -> pd.DataFrame:
    """Keep AMR hits meeting the minimum coverage and identity thresholds."""
    results = results[
        (results["% Coverage of reference"] >= min_cov)
        & (results["% Identity to reference"] >= min_id)
    ]
    return results


def create_rules(species: str, reference_folder: str) -> list:
    """Load species rules and convert them to inference-rule objects."""
    rules = find_rules(species=species, reference_folder=reference_folder)
    try:
        rules = [InferRules(**r) for r in rules]
    except Exception as e:
        log.warning(f"No rules available for {species} in {reference_folder}: {e}")
        rules = []
    return rules


def gdst(
    results: pd.DataFrame,
    species: str,
    reference_folder: str,
    dflt_result: str = "Susceptible (default)",
) -> list:
    """Evaluate species rules and return inferred susceptibility results."""

    sid = results.iloc[0].get("sample_id", "unknown")
    resultsmooshed = combine_results(result=results.to_dict(orient="records"))
    rules = create_rules(species=species, reference_folder=reference_folder)
    if rules == []:
        log.info(
            f"No gDST will be provided for {sid}. There are no rules available for {species}."
        )
        return []
    gdst_final = []
    for row in resultsmooshed:
        gdst_results = {"sample_id": sid, "species": species}
        data = {row: resultsmooshed[row]}
        ctx = create_cel_context(data=data, name="row")
        rlt = {}
        for rule in rules:
            if row in rule.rule:
                if rule.drugname not in rlt:
                    rlt[rule.drugname] = {
                        "drugname": rule.drugname,
                        "mechanisms": [],
                        "inferred": [],
                        "rule_id": [],
                        "rule_version": [],
                        "source": [],
                    }
                if evaluate_rule(rule=rule.rule, ctx=ctx):
                    mechs = []
                    for key in results:
                        if key in rule.rule:
                            vals = (
                                resultsmooshed[key]
                                if isinstance(resultsmooshed[key], list)
                                else [resultsmooshed[key]]
                            )
                            keys = [i for i in vals if i in rule.rule]
                            keys = []
                            for i in vals:
                                for j in i.split("_"):
                                    if j in rule.rule:
                                        keys.append(j)
                            mechs = []
                            for k in keys:
                                tmp = (
                                    results[
                                        results["abritamr_accession_key"].str.contains(
                                            k, na=False
                                        )
                                    ]["abritamr_mechanism"]
                                    .unique()
                                    .tolist()
                                )
                                mechs.extend(tmp)
                    rlt[rule.drugname]["mechanisms"].extend(mechs)
                    rlt[rule.drugname]["inferred"].append(rule.inferred)
                    rlt[rule.drugname]["rule_id"].append(rule.rule_id)
                    rlt[rule.drugname]["rule_version"].append(rule.rule_version)
                    rlt[rule.drugname]["source"].append(rule.source)
        for drug in rlt:
            if rlt[drug]["inferred"] == []:
                rlt[drug]["inferred"] = [dflt_result]
            else:
                rlt[drug]["inferred"] = [
                    sorted(
                        rlt[drug]["inferred"],
                        key=lambda x: (
                            priority_gdst()[x[0].upper()]
                            if x[0].upper() in priority_gdst()
                            else -1
                        ),
                        reverse=True,
                    )[0]
                ]
            for key in ["mechanisms", "rule_id", "rule_version", "source", "inferred"]:
                rs = (
                    ";".join(rlt[drug][key])
                    if rlt[drug][key] != [] or set(rlt[drug][key]) != {"-"}
                    else "-"
                )
                rlt[drug][key] = rs

        gdst_results = {
            "sample_id": sid,
            "species": species,
            "results": list(rlt.values()),
        }
        gdst_final.append(gdst_results)

    return gdst_final


def gdst_results_to_df_wide(gdst_results: list) -> pd.DataFrame:
    """Convert inferred results to one row per sample with drug-specific columns."""
    rows = []
    for res in gdst_results:
        row = {
            "sample_id": res["sample_id"],
            "species": res["species"],
        }
        for drug_res in res["results"]:
            row[f"{drug_res['drugname']}_gDST"] = drug_res["inferred"]
            row[f"{drug_res['drugname']}_mechanisms"] = drug_res["mechanisms"]
            row[f"{drug_res['drugname']}_rule_id"] = drug_res["rule_id"]
            row[f"{drug_res['drugname']}_rule_version"] = drug_res["rule_version"]
            row[f"{drug_res['drugname']}_source"] = drug_res["source"]
        rows.append(row)
    dr_cols = [
        f"{drug}_{suffix}"
        for drug in set(
            drug_res["drugname"] for res in gdst_results for drug_res in res["results"]
        )
        for suffix in ["gDST", "mechanisms", "rule_id", "rule_version", "source"]
    ]
    cols_wide = ["sample_id", "species"] + sorted(dr_cols)
    return pd.DataFrame(rows)[cols_wide].sort_values(by=["sample_id"])


def gdst_results_to_df_long(gdst_results: list) -> pd.DataFrame:
    """Convert inferred results to one row per sample and drug."""
    rows = []
    for res in gdst_results:
        for drug_res in res["results"]:
            row = {
                "sample_id": res["sample_id"],
                "species": res["species"],
                "drugname": drug_res["drugname"],
                "gDST": drug_res["inferred"],
                "rule_id": drug_res["rule_id"],
                "rule_version": drug_res["rule_version"],
                "source": drug_res["source"],
                "mechanisms": drug_res["mechanisms"],
            }
            rows.append(row)
    cols_long = [
        "sample_id",
        "species",
        "drugname",
        "mechanisms",
        "gDST",
        "rule_id",
        "rule_version",
        "source",
    ]

    return pd.DataFrame(rows)[cols_long].sort_values(by=["drugname"])
