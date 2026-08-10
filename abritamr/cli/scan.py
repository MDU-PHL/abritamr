import pathlib, argparse, sys, os, logging,json




def scan_args(subparsers):
    parser_sub_scan = subparsers.add_parser('scan', help='Run amrfinder and apply abritamr drug classes. Outputs a file with rows populated with per-gene information.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_scan.add_argument(
        "--assembly",
        "-asm",
        nargs="*",
        # default="",
        help="Assembly file to use as input (*.fa*, *.gbk *.fa*.gz, *.gbk.gz). Can be multiple or wildcard.",
    )
    parser_sub_scan.add_argument(
        "--amrfinderplus",
        "-afp",
        nargs="*",
        # default="",
        help="EXPERIMENTAL - USE WITH CAUTION. AMR finder output file. Please note unexpected behaviour may result were versions and databases differ from abritamr. Can be mutliple",
    )
    parser_sub_scan.add_argument(
        '--sample-id',
        '-s',
        help = "sample identifier, this will be used to name output files and in line list reports. If not supplied, path to input file will be used as sample id",
        # default = ""
    )
    parser_sub_scan.add_argument(
        '--format',
        '-f',
        help = "Output format",
        choices=['csv', 'tab'],
        default = "csv"
    )
    parser_sub_scan.add_argument(
        '--output',
        '-o',
        help = "Filename to save output - default stdout.",

    )
    parser_sub_scan.add_argument(
        '--species',
        '-sp',
        help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence.",
       
    )
    parser_sub_scan.add_argument(
        '--threads',
        '--cpus',
        help="Number of max CPU cores to run.",
        default=1
    )
    parser_sub_scan.add_argument(
        "--min-identity",
        default = 0.9,
        help ="Minimum identity to reference gene for reporting a match."
    )
    parser_sub_scan.add_argument(
        "--min-coverage",
        default = 0.5,
        help ="Minimum coverage of the reference gene for reporting a complete match."
    )
    
    return parser_sub_scan