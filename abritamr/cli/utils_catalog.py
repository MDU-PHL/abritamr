import pathlib, argparse, sys, os, logging,json


spcs = ["all"]

with open(f"{pathlib.Path(__file__).parent.parent / 'configs'/ 'amr_rules_config.json'}", "r") as s:
    sp = json.load(s)

spcs.extend(sp['species'])


def catalog_args(subparsers):
    # parser_sub_utils_catalog = subparsers.add_parser('catalog', help='Generate the reference gene catalog for abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_catalog = subparsers.add_parser('generate_gene_catalog', help='Create a reference gene catalogue ONLY for abritamr (this will not update any rules for genomic DST)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_catalog.add_argument(
        "--latest",
        "-l",
        action="store_true",
        help="Download the latest reference gene catalogue. Please note that if you have not installed or updated the AMRfinder database you may see unexpected behaviour. Please ensure you have the latest version of the AMRfinder database installed before running abritamr with this catalogue."
    )
    parser_sub_catalog.add_argument(
        "--catalog",
        default=f"",
        help="Path to a previous downloaded reference gene catalog. If you do not supply this, the latest version will be downloaded and used to create a new catalog."
    )
    parser_sub_catalog.add_argument(
        "--class-definitions",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_class_definitions.csv'}",
        help="Path to definitions for class curation."
    )
    parser_sub_catalog.add_argument(
        "--amrtyping-definitions",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_reporting.csv'}",
        help="Path to criteria for reporting and amr typing."
    )
    parser_sub_catalog.add_argument(
        "--previous-catalog",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_reference_catalog.csv'}",
        help="Path to where previous reference gene catalog. Used to identify updates or changes for checking."
    )
    parser_sub_catalog.add_argument(
        "--output-dir",
        default=f"{pathlib.Path.cwd() / 'abritamr_db'}",
        help="Path to save the abritamr configuration files."
    )
    
    return parser_sub_catalog