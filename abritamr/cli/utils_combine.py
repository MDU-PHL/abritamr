"""Argument parser for combining output files into a single file"""

import argparse, pathlib


def combine_args(subparsers):
    """Add the combine utility and its options to the parser."""
    # parser_sub_utils_catalog = subparsers.add_parser('catalog', help='Generate the reference gene catalog for abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_combine = subparsers.add_parser(
        "combine",
        help="Combine abritamr outputs for each stage into single files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_sub_combine.add_argument(
        "--input-directory",
        default=f"{pathlib.Path.cwd()}",
        help="Path to look in for abritamr generated outputs.",
    )
    parser_sub_combine.add_argument(
        "--output-directory",
        default=f"{pathlib.Path.cwd()}",
        help="Directory to save output files to.",
    )
    parser_sub_combine.add_argument(
        "--prefix", default="abritamr", help="Prefix to use for naming files."
    )

    return parser_sub_combine
