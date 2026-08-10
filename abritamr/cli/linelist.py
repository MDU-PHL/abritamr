import pathlib, argparse, sys, os, logging,json




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
    parser_sub_llist.add_argument(
        '--output',
        '-o',
        help = "Filename to save output - default stdout.",

    )
    parser_sub_llist.add_argument(
        '--format',
        '-f',
        help = "Output format",
        choices=['csv', 'tab'],
        default = "csv"
    )
    parser_sub_llist.add_argument(
        '--viewtype',
        '-vt',
        help = "For line list, full will print out all drug classes available in the reference set, compact will only print out hat is detected.",
        choices=['full', 'compact'],
        default = 'compact'
    )
    parser_sub_llist.add_argument(
        "--min-identity",
        default = 0.9,
        help ="Minimum identity to reference gene for reporting a match."
    )
    parser_sub_llist.add_argument(
        "--min-coverage",
        default = 0.5,
        help ="Minimum coverage of the reference gene for reporting a complete match."
    )
    parser_sub_llist.add_argument(
        '--genesonly',
        action='store_true',
        help = "Set to if you want to limit reportable to only be genes (excludes reporting of any SNPs)."
    )

    return parser_sub_llist