import pathlib, argparse, sys, os, logging,json
from abritamr.cli.basic_args import inputs, references, detection_args, basic_output



def llist_args(subparsers):


    parser_sub_llist = subparsers.add_parser('linelist', help='Generate a linelist report summarising genes observed by abritamr drugclass or criteria supplied.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_llist.add_argument(
        '--amr',
        '-a',
        help = "result file(s) from 'abritamr scan'. Default stdin",
        metavar='FILE', 
        nargs='?',
        default = sys.stdin
    )
    parser_sub_llist = basic_output(parser = parser_sub_llist)
    parser_sub_llist.add_argument(
        '--viewtype',
        '-vt',
        help = "For line list, full will print out all drug classes available in the reference set, compact will only print out hat is detected.",
        choices=['full', 'compact'],
        default = 'compact'
    )
    parser_sub_llist.add_argument(
        '--genesonly',
        action='store_true',
        help = "Set to if you want to limit reportable to only be genes (excludes reporting of any SNPs)."
    )
    parser_sub_llist = detection_args(parser = parser_sub_llist)

    return parser_sub_llist