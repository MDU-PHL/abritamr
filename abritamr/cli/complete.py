import pathlib, argparse, sys, os, logging,json




def complete_args(subparsers):
    parser_sub_complete = subparsers.add_parser('complete', help='Run the complete suite of abritamr functions. This will create a folder for each sample and generate a linelist report and inferred antibiogram (if supported for your species).', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_complete.add_argument(
        "--assembly",
        "-asm",
        default="",
        help="Assembly file to use as input (*.fa*, *.gbk *.fa*.gz, *.gbk.gz). If you would like to run with more than one sample - please use the --multi option.",
    )
    parser_sub_complete.add_argument(
        "--amrfinderplus",
        "-afp",
        default="",
        help="EXPERIMENTAL - USE WITH CAUTION. AMR finder output file. Please note unexpected behaviour may result were versions and databases differ from abritamr. If you would like to run with more than one sample - please use the --multi option.",
    )
    parser_sub_complete.add_argument(
        '--species',
        '-sp',
        help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence. If you would like to run with more than one sample - please use the --multi option.",
       
    )
    parser_sub_complete.add_argument(
        '--sample-id',
        '-s',
        help = "sample identifier, this will be used to name output folders and in line list reports. If you would like to run with more than one sample - please use the --multi option.",
    )
    parser_sub_complete.add_argument(
        '--multi',
        help = "Tab-delimited file with columns sample_id,assembly (or amrfinder) and species (optional).",
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
        '--prefix',
        '-op',
        help = "Prefix for collated output files.",
        default = "abritamr"

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
    parser_sub_complete.add_argument(
        '--threads',
        '--cpus',
        help="Number of max CPU cores to run.",
        default=1
    )
    parser_sub_complete.add_argument(
        "--min-identity",
        default = 0.9,
        help ="Minimum identity to reference gene for reporting a match."
    )
    parser_sub_complete.add_argument(
        "--min-coverage",
        default = 0.5,
        help ="Minimum coverage of the reference gene for reporting a complete match."
    )
    parser_sub_complete.add_argument(
        "--reportable-config",
        "-rc",
        default = f"{pathlib.Path(__file__).parent / 'configs' / 'abritamr_reporting.csv'}",
        help = "Path to config that defines the criteria for highlighting relevant genes and mechanisms for reporting. Supplying a new config will override the default abritamr criteria."
    )
    
    return parser_sub_complete