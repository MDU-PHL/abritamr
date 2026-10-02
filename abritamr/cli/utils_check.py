"""Argument parser configuration for checking installation."""

import argparse


def check_args(subparsers):
    """Add the check utility and its options to the parser."""
    # parser_sub_utils_catalog = subparsers.add_parser('catalog', help='Generate the reference gene catalog for abritamr.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_check = subparsers.add_parser(
        "check",
        help="Check that abritamr dependencies are installed.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    return parser_sub_check
