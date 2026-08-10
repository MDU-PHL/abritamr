import pathlib, argparse, sys, os, logging,json

from abritamr.cli.basic_args import inputs, references, detection_args, basic_output


def status_args(subparsers):

    parser_sub_status = subparsers.add_parser('amr_status', help='Determine status of genes recovered from abritamr scan.  Outputs a file with rows populated with per-gene information.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_status.add_argument(
        '--amr',
        '-a',
        help = "result file(s) from 'abritamr scan'. Default stdin",
        # metavar='FILE', 
        nargs='?',
        default = sys.stdin
    )
    parser_sub_status.add_argument(
        '--species',
        '-sp',
        help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence.",
       
    )
    parser_sub_status = basic_output(parser = parser_sub_status)
    parser_sub_status = detection_args(parser = parser_sub_status)
    parser_sub_status = references(parser = parser_sub_status)

    return parser_sub_status
