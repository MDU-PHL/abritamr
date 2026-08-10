import pathlib, argparse, sys, os, logging,json




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
    parser_sub_matrix.add_argument(
        '--output',
        '-o',
        help = "Filename to save output - default stdout.",

    )
    parser_sub_matrix.add_argument(
        '--facet',
        help = "Feature to facet by for generating a matrix",
        choices=['abritamr_subclass', 'abritamr_class', 'gene'],
        default = "abritamr_subclass"
    )
    parser_sub_matrix.add_argument(
        '--format',
        '-f',
        help = "Output format",
        choices=['csv', 'tab'],
        default = "csv"
    )
    parser_sub_matrix.add_argument(
        "--min-identity",
        default = 0.9,
        help ="Minimum identity to reference gene for reporting a match."
    )
    parser_sub_matrix.add_argument(
        "--min-coverage",
        default = 0.5,
        help ="Minimum coverage of the reference gene for reporting a complete match."
    )

    return parser_sub_matrix