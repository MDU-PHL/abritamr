import pathlib, argparse, sys, os, logging,json




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
        '--output',
        '-o',
        help = "Filename to save output - default stdout.",

    )
    parser_sub_status.add_argument(
        '--format',
        '-f',
        help = "Output format",
        choices=['csv', 'tab'],
        default = "csv"
    )
    parser_sub_status.add_argument(
        '--species',
        '-sp',
        help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence.",
       
    )
    parser_sub_status.add_argument(
        "--reportable-config",
        "-rc",
        default = f"{pathlib.Path(__file__).parent / 'configs' / 'abritamr_reporting.csv'}",
        help = "Path to config that defines the criteria for highlighting relevant genes and mechanisms for reporting. Supplying a new config will override the default abritamr criteria."
    )

    return parser_sub_status
