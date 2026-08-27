import pathlib, argparse, sys, os, logging,json

from abritamr.cli.basic_args import inputs, references, detection_args, basic_output




def infer_args(subparsers):

    parser_sub_infer = subparsers.add_parser('infer', help='Infer phenotype from detected AMR mechanisms.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_infer.add_argument(
        '--amr',
        '-a',
        help = "result file(s) from 'abritamr scan'. Default stdin",
        # metavar='FILE', 
        nargs='?',
        default = sys.stdin
    )
    parser_sub_infer = basic_output(parser = parser_sub_infer)
    parser_sub_infer.add_argument(
        '--dflt_result',
        help = "Default result to report if no rules are matched. Default is 'Susceptible (default)'.",
        default = "Susceptible (default)"
    )
    parser_sub_infer.add_argument(
            '--reporttype',
            '-rt',
            help = "For gDST report, long will print out the drugs available for the species as rows, wide will print out the drugs available for the species as columns.",
            choices=['long', 'wide'],
            default = 'long'
    )
    parser_sub_infer = detection_args(parser = parser_sub_infer)

