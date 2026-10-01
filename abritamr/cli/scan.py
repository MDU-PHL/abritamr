"""Argument parser configuration for AMRFinder scanning."""

import pathlib, argparse, sys, os, logging, json

from abritamr.cli.basic_args import (
    inputs,
    references,
    rules,
    detection_args,
    basic_output,
)


def scan_args(subparsers):
    """Add the scan command and its options to the parser."""
    parser_sub_scan = subparsers.add_parser(
        "scan",
        help="Run amrfinder and apply abritamr drug classes. Outputs a file with rows populated with per-gene information.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_sub_scan = inputs(parser=parser_sub_scan)
    parser_sub_scan = basic_output(parser=parser_sub_scan)
    parser_sub_scan = rules(parser=parser_sub_scan)
    parser_sub_scan = detection_args(parser=parser_sub_scan)
    parser_sub_scan = references(parser=parser_sub_scan)
    return parser_sub_scan
