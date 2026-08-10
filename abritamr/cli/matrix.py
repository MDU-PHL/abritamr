import pathlib, argparse, sys, os, logging,json
from abritamr.cli.basic_args import inputs, references, detection_args, basic_output



def matrix_args(subparsers):


    parser_sub_matrix = subparsers.add_parser('matrix', help='Generate a presence/absence table based on genes or drugclasses', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_matrix.add_argument(
        '--amr',
        '-a',
        help = "result file(s) from 'abritamr scan'. Default stdin",
        metavar='FILE', 
        nargs='?',
        default = sys.stdin
    )
    parser_sub_matrix = basic_output(parser = parser_sub_matrix)
    parser_sub_matrix.add_argument(
        '--facet',
        help = "Feature to facet by for generating a matrix",
        choices=['abritamr_subclass', 'abritamr_class', 'gene'],
        default = "abritamr_subclass"
    )
    parser_sub_matrix = detection_args(parser = parser_sub_matrix)
    parser_sub_matrix = references(parser_sub_matrix)

    return parser_sub_matrix