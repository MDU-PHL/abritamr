import pathlib, argparse, sys, os, logging,json




def catalog_args(subparsers):
    # parser_sub_utils_catalog = subparsers.add_parser('catalog', help='Generate the reference gene catalog for abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_catalog = subparsers.add_parser('catalog', help='Create a reference gene catalogue for abritamr. This will generate subclass and amr status columns to an AMRfinder reference gene catalogue for use in abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_catalog.add_argument(
        "--catalog",
        "-c",
        default="",
        help="Path to reference gene catalogue to use. If none supplied the version compatible with your version of abritamr will be downloaded."
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
        "--inferrence-definitions",
        default=f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_inferrence.csv'}",
        help="Path to criteria for AMR AST inferrence."
    )
    parser_sub_catalog.add_argument(
        "--previous-catalog",
        default=f"{pathlib.Path(__file__).parent.parent / 'db' / 'refgenes_latest.csv'}",
        help="Path to were previous reference gene catalog. Used to identify updates or changes for checking."
    )
    parser_sub_catalog.add_argument(
        "--output",
        default=f"{pathlib.Path.cwd() / 'refgenes.csv'}",
        help="Path to save you catalog. Please check this file for correctness."
    )