import click
import pathlib
import json
import logging
import pandas as pd
import sys

from abritamr.logger import log
from abritamr.catalog import wrangle_catalog
from abritamr.parse_gdstrules import add_rules_to_existing, get_rules_for_species as get_amrrules_for_species


def create_db_folder(pth:str) -> bool:
    try:
        pathlib.Path(pth).absolute().mkdir(parents=True, exist_ok=True)
        return True
    except Exception as e:
        log.critical(f"Cannot create path for database files folder {pth}. {e}")
        raise SystemExit(1)

def update_catalog(args) -> bool:

    if create_db_folder(args.output_dir):
        log.info(f"Database folder created at {args.output_dir}")

        try:

            created = wrangle_catalog(
                catalog = args.catalog,
                previous_catalog = args.previous_catalog,
                amrtyping_definitions = args.amrtyping_definitions,
                class_definitions = args.class_definitions)
            created.to_csv(f"{args.output_dir}/01_abritamr_reference_gene_catalog.csv", index = False)
            log.info(f"The basic reference gene catalogue has been created: {args.output_dir}/01_abritamr_reference_gene_catalog.csv")


            
        except Exception as e:
            log.critical(f"Looks like something has gone wrong with generating the reference catalogue. The following error was reported : {e}. Please try again.")
            raise SystemExit(1)

def update_rules(args) -> bool:
    
    if create_db_folder(args.output_dir):
        log.info(f"Database folder created at {args.output_dir}")

        try:
            rules_dict = {}
            if not args.no_amrrules:
                
                rules_dict = get_amrrules_for_species(
                    evidence_grade = args.evidence_grade,
                    species = args.species,
                    rules_dict = rules_dict)

            rules_dict = add_rules_to_existing(
                rules = rules_dict,
                additional_rules = args.inference_definitions
            )
            for sp in rules_dict:
                fn = sp.replace(" ", "_")
                pd.DataFrame(rules_dict[sp]).to_csv(f"{args.output_dir}/02_abritamr_{fn}_rules.csv", index = False)
                log.info(f"Rules for {sp} have been created: {args.output_dir}/02_abritamr_{fn}_rules.csv")

            
        except Exception as e:
            log.critical(f"Looks like something has gone wrong with generating the rules. The following error was reported : {e}. Please try again.")
            raise SystemExit(1) 