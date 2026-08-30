import pathlib, argparse, sys, os, logging,json


spcs = ["all"]

with open(f"{pathlib.Path(__file__).parent.parent / 'configs'/ 'amr_rules_config.json'}", "r") as s:
    sp = json.load(s)

spcs.extend(sp['species'])


def database_args(subparsers):
    # parser_sub_utils_catalog = subparsers.add_parser('catalog', help='Generate the reference gene catalog for abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_database = subparsers.add_parser('generate_database', help='Create a reference gene catalogue for abritamr and genomic DST rules from the latest sources. This will download latest versions of the amrfinder gene catalog and genomic DST rules from the AMRrules github.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_database.add_argument(
        "--latest",
        "-l",
        action="store_true",
        help="Download the latest reference gene catalogue. Please note that if you have not installed or updated the AMRfinder database you may see unexpected behaviour. Please ensure you have the latest version of the AMRfinder database installed before running abritamr with this catalogue."
    )
    parser_sub_database.add_argument(
        "--catalog",
        default=f"",
        help="Path to a previous downloaded reference gene catalog. If you do not supply this, the latest version will be downloaded and used to create a new catalog."
    )
    parser_sub_database.add_argument(
        "--class-definitions",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_class_definitions.csv'}",
        help="Path to definitions for class curation."
    )
    parser_sub_database.add_argument(
        "--amrtyping-definitions",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_reporting.csv'}",
        help="Path to criteria for reporting and amr typing."
    )
    parser_sub_database.add_argument(
        "--inference-definitions",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_inference.csv'}",
        help="Path to criteria for AMR AST inference. This file may contain rules for multiple species (please see documentation for the correct format)."
    )
    parser_sub_database.add_argument(
        "--replace-amrverse-rules",
        action='store_true',
        help = "Set to if you want to replace the AMRverse rules with your own rules. Please see documentation for the correct format. Default behaviour is to append your rules to the AMRverse rules."
    )
    parser_sub_database.add_argument(
        "--previous-catalog",
        default=f"{pathlib.Path(__file__).parent.parent / 'db' / 'abritamr_reference_catalog.csv'}",
        help="Path to where previous reference gene catalog. Used to identify updates or changes for checking."
    )
    parser_sub_database.add_argument(
            "--no-amrrules",
            action='store_true',
            help = "Set to if you DO NOT want to generate rules from AMR rules (https://github.com/AMRverse/AMRrules)"
        )
    parser_sub_database.add_argument(
        "--species",
        "-s",
        choices = spcs,
        help="Species to generate rules for (only required if --no-amrrules is not set).",
        default = "all"
    )
    parser_sub_database.add_argument(
        "--evidence_grade",
        choices = ['very low', 'low', 'moderate', 'high'],
        help="Minimum evidence grade to include (only required if --no-amrrules is not set).",
        default = "very low"
    )
    parser_sub_database.add_argument(
        "--output-dir",
        default=f"{pathlib.Path.cwd() / 'abritamr_db'}",
        help="Path to save the abritamr configuration files."
    )
    
    return parser_sub_database