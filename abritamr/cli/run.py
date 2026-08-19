import pathlib, argparse, sys, os, logging,json

from abritamr.cli.basic_args import inputs, references, detection_args


def run_args(subparsers):
    parser_sub_complete = subparsers.add_parser('complete', help='Run the complete suite of abritamr functions. This will create a folder for each sample and generate a linelist report and inferred antibiogram (if supported for your species).', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_complete = inputs(parser = parser_sub_complete)
    parser_sub_complete.add_argument(
        '--prefix',
        '-p',
        help = "Prefix for collated output files.",
        default = "abritamr"

    )
    parser_sub_complete.add_argument(
        '--multi',
        help = "Tab-delimited file with columns sample_id, path to assembly (and/or amrfinder outputs) and species (optional).",
        default = ""
    )
    parser_sub_complete.add_argument(
        '--format',
        '-f',
        help = "Output format",
        choices=['csv', 'tab'],
        default = "csv"
    )
    parser_sub_complete.add_argument(
        '--viewtype',
        '-vt',
        help = "For line list, full will print out all drug classes available in the reference set, compact will only print out hat is detected.",
        choices=['full', 'compact'],
        default = 'compact'
    )
    parser_sub_complete.add_argument(
        '--facet',
        help = "Feature to facet by for generating a matrix",
        choices=['abritamr_subclass', 'abritamr_class', 'gene'],
        default = "abritamr_subclass"
    )
    parser_sub_complete.add_argument(
        '--workdir',
        '-w',
        help = "Working directory where output folders and results will be created.",
        default = f"{pathlib.Path.cwd()}"
    )
    parser_sub_complete.add_argument(
        '--no-keep',
        action='store_true',
        help = "Set to if you DO NOT want to keep intermediate files."
    )
    parser_sub_complete.add_argument(
        '--genesonly',
        action='store_true',
        help = "Set to if you want to limit reportable to only be genes (excludes reporting of any SNPs)."
    )
    
    parser_sub_complete = detection_args(parser = parser_sub_complete)
    parser_sub_catalog = references(parser = parser_sub_complete)
    
    
    return parser_sub_complete