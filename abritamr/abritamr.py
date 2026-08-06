
import pathlib, argparse, sys, os, logging,json

from abritamr.commands import scan,linelist,amr_status
from abritamr.version import __version__, db
from abritamr.utils import output_results
from abritamr.logger import log
"""
abritamr is designed to implement AMRFinder and parse the results compatible for MDU use. It may also be used for other purposes where the format of output is compatible

"""
def run_scan(args):
    log.info(f"Running scan to identify genes and and drugclasses")
    amr = scan.scan(args)
    output_results(df = amr, output = args.output, _format = args.format)

def run_linelist(args):
    log.info(f"Running linelist to collate single sample results into a linelist for reporting.")
    result = linelist.linelist(args)
    # print(linelist)
    output_results(df = result, output = args.output, _format = args.format)

def run_status(args):
    
    log.info(f"Determining the amr status of the sequence based on genes and drug classes.")
    amr = amr_status.amr_status(args)
    output_results(df = amr, output = args.output, _format = args.format)

def cli():

    parser = argparse.ArgumentParser(
        description=f"****abritamr - AMR gene detection and reporting pipeline - version {__version__}****", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    
    subparsers = parser.add_subparsers(help="What would you like to do?")
    parser_sub_scan = subparsers.add_parser('scan', help='Run amrfinder and apply reporting logic. Outputs a file with rows populated with per-gene information.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    
    
    parser_sub_status = subparsers.add_parser('amr_status', help='Determine status of genes recovered from abritamr scan', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser_sub_llist = subparsers.add_parser('linelist', help='Generate a linelist report summarising genes observed by abritamr durgclass', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
# @click.option(
#     '--format',
#     '-f',
#     help = "Output format",
#     type = click.Choice(["csv", "tab"]),
#     default = "csv",
#     show_default = True
# )
# @click.option(
#     '--output',
#     '-o',
#     help = "Filename to save output - default stdout.",
#     default = "",
#     show_default = True
# )
# @click.option(
#     '--species',
#     '-sp',
#     help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence.",
#     default = "",
#     show_default = True
# )
# @click.option(
#     '--viewtype',
#     '-vt',
#     default = "full",
#     show_default=True,
#     type = click.Choice(["full", "compact"]),
#     help = "Format of abritamr report. Only applicable if --summary is set. Default - full will output all results for the sequence. Compact will only output summarised results."
# )
# @click.option(
#     '--genesonly',
#     is_flag=True,
#     default=False,
#     show_default=True,
#     help = "If only whole genes, not SNPs should be reported. Only applicable if --summary is set."
# )
# @click.option(
#     "--min-identity",
#     default = 0.9,
#     help ="Minimum identity to reference gene for reporting a match.",
#     show_default = True
# )
# @click.option(
#     "--min-coverage",
#     default = 0.5,
#     help ="Minimum coverage of the reference gene for reporting a complete match.",
#     show_default = True
# )
# @click.option(
#     "--reportable-config",
#     "-rc",
#     default = f"{pathlib.Path(__file__).parent / 'configs' / 'abritamr_reporting.csv'}",
#     help = "Path to config that defines the criteria for highlighting relevant genes and mechanisms for reporting. Supplying a new config will override the default abritamr criteria.",
#     show_default = True
# )




    parser_sub_scan.set_defaults(func=run_scan)
    parser_sub_status.set_defaults(func = run_status)
    parser_sub_llist.set_defaults(func = run_linelist)
    # parser_update.set_defaults(func = update_db)
    args = parser.parse_args()
    
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
    elif len(sys.argv) <= 2 and sys.argv[1] == 'scan':
        parser_sub_scan.print_help(sys.stderr)
    # elif len(sys.argv) <= 2 and sys.argv[1] == 'report':
    #     parser_mdu.print_help(sys.stderr)
    else:
        args.func(args)
	
if __name__ == '__main__':
    cli()
