import pathlib, argparse, sys, os, logging,json


def references(parser):

    parser.add_argument(
        "--reference-folder",
        default = f"{pathlib.Path(__file__).parent.parent / 'configs' }",
        help = "Path to catalog for applying classes and amr types"
    )
    parser.add_argument(
        "--reference-catalog",
        default = f"{pathlib.Path(__file__).parent.parent / 'configs' / 'abritamr_reference_catalog.csv'}",
        help = "Path to catalog for applying classes and amr types"
    )

    return parser


def inputs(parser):

    parser.add_argument(
        "--contigs",
        "-c",
        nargs="*",
        # default="",
        help="Assembly file to use as input (*.fa*, *.gbk *.fa*.gz, *.gbk.gz). Can be multiple or wildcard.",
    )
    parser.add_argument(
        "--amrfinderplus",
        "-afp",
        nargs="*",
        # default="",
        help="EXPERIMENTAL - USE WITH CAUTION. AMR finder output file. Please note unexpected behaviour may result were versions and databases differ from abritamr. Can be mutliple",
    )
    parser.add_argument(
        '--sample-id',
        '-s',
        help = "sample identifier, this will be used to name output files and in line list reports. If not supplied, path to input file will be used as sample id",
        # default = ""
    )
    parser.add_argument(
        '--species',
        '-sp',
        help = "Species from which assemblies were derived. If not supplied, will be guessed using sourmash and used for SNP detection and inference.",
       
    )
    parser.add_argument(
        '--threads',
        '--cpus',
        help="Number of max CPU cores to run.",
        default=1
    )
    return parser


def detection_args(parser):
    parser.add_argument(
        "--min-identity",
        default = 0.9,
        help ="Minimum identity to reference gene for reporting a match."
    )
    parser.add_argument(
        "--min-coverage",
        default = 0.5,
        help ="Minimum coverage of the reference gene for reporting a complete match."
    )

    return parser

def basic_output(parser):

    parser.add_argument(
        '--output',
        '-o',
        help = "Filename to save output - default stdout.",

    )
    parser.add_argument(
        '--format',
        '-f',
        help = "Output format",
        choices=['csv', 'tab'],
        default = "csv"
    )

    return parser