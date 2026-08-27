import pathlib, argparse, sys, os, logging,json


spcs = ["all"]

with open(f"{pathlib.Path(__file__).parent.parent / 'configs'/ 'amr_rules_config.json'}", "r") as s:
    sp = json.load(s)

spcs.extend(sp['species'])


def rules_args(subparsers):
    # parser_sub_utils_catalog = subparsers.add_parser('catalog', help='Generate the reference gene catalog for abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_rules = subparsers.add_parser('generate_rules', help='Generate ONLY abritamr compatible rules for reporting genomic DST.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_rules.add_argument(
        "--catalog",
        "-c",
        default=f"",
        help="Path to reference gene catalogue to use. If none supplied the version compatible with your version of abritamr will be downloaded."
    )
    parser_sub_rules.add_argument(
            "--inference-definitions",
            default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_inference.csv'}",
            help="Path to criteria for AMR AST inference. This file may contain rules for multiple species (please see documentation for the correct format)."
    )
    parser_sub_rules.add_argument(
        "--no-amrrules",
        action='store_true',
        help = "Set to if you DO NOT want to generate rules from AMR rules (https://github.com/AMRverse/AMRrules)"
    )
    parser_sub_rules.add_argument(
        "--species",
        "-s",
        choices = spcs,
        help="Species to generate rules for.",
        default = "all"
    )
    parser_sub_rules.add_argument(
        "--evidence_grade",
        choices = ['very low', 'low', 'moderate', 'high'],
        help="Minimum evidence grade to include.",
        default = "very low"
    )
    parser_sub_rules.add_argument(
            "--output-dir",
            default=f"{pathlib.Path.cwd() / 'abritamr_db'}",
            help="Path to save the abritamr configuration files."
        )
    return parser_sub_rules