"""Create and update local reference catalogs and inference rules."""

import pathlib
import pandas as pd

from abritamr.logger import log
from abritamr.catalog import wrangle_catalog
from abritamr.parse_gdstrules import add_rules_to_existing, get_amrrules_for_species


def create_db_folder(pth: str) -> bool:
    """Create the destination directory for generated database files."""
    try:
        pathlib.Path(pth).absolute().mkdir(parents=True, exist_ok=True)
        return True
    except Exception as e:
        log.critical(f"Cannot create path for database files folder {pth}. {e}")
        raise SystemExit(1)


def update_catalog(args) -> bool:
    """Build the reference gene catalog from the supplied definitions."""
    if create_db_folder(args.output_dir):
        log.info(f"Database folder created at {args.output_dir}")

        created = wrangle_catalog(
            catalog=args.catalog,
            previous_catalog=args.previous_catalog,
            amrtyping_definitions=args.amrtyping_definitions,
            class_definitions=args.class_definitions,
            output_dir=args.output_dir,
        )
        pd.DataFrame(created).to_csv(
            f"{args.output_dir}/01_abritamr_reference_gene_catalog.csv", index=False
        )
        log.info(
            f"The basic reference gene catalogue has been created: {args.output_dir}/01_abritamr_reference_gene_catalog.csv"
        )


def update_rules(args) -> bool:
    """Generate species inference rules and save them in the output folder."""
    if create_db_folder(args.output_dir):
        log.info(f"Database folder created at {args.output_dir}")
        log.info(
            f"Generating rules for {args.species} with evidence grade {args.evidence_grade}"
        )

        rules_dict = {}
        if not args.no_amrrules:
            rules_dict = get_amrrules_for_species(
                evidence_grade=args.evidence_grade,
                species=args.species,
                rules_dict=rules_dict,
                output_dir=args.output_dir,
            )

        rules_dict = add_rules_to_existing(
            rules=rules_dict, additional_rules=args.inference_definitions
        )
        for sp in rules_dict:
            fn = sp.split()
            tmp = pd.DataFrame(rules_dict[sp])
            tmp = tmp[tmp["drugname"] != "-"]
            if len(fn) <= 2:
                tmp = tmp[tmp["drugname"] != "-"]
                if not tmp.empty:
                    tmp.to_csv(
                        f"{args.output_dir}/02_abritamr_{'_'.join(fn)}_rules.csv",
                        index=False,
                    )
                    log.info(
                        f"Rules for {sp} have been created: {args.output_dir}/02_abritamr_{'_'.join(fn)}_rules.csv"
                    )
            else:
                for s in tmp["amrrules_species"].unique():
                    tmp2 = tmp[
                        (tmp["amrrules_species"] == s) & (tmp["drugname"] != "-")
                    ]
                    if not tmp2.empty:
                        tmp2.to_csv(
                            f"{args.output_dir}/02_abritamr_{'_'.join(s.split())}_rules.csv",
                            index=False,
                        )
                        log.info(
                            f"Rules for {s} have been created: {args.output_dir}/02_abritamr_{'_'.join(s.split())}_rules.csv"
                        )


def generate_database(args) -> bool:
    """Generate both the reference gene catalog and inference rules."""
    update_catalog(args=args)
    update_rules(args=args)
    log.info(
        f"Database generation complete. Please check the output folder {args.output_dir} for the generated files."
    )
