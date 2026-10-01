"""Command-line entry point for running abriTAMR workflows."""

import argparse, sys

from abritamr.commands import (
    scan,
    linelist,
    amr_status,
    matrix,
    run,
    update_database,
    utils_catalog,
    utils_rules,
    infer,
)
from abritamr.cli.run import run_args
from abritamr.cli.scan import scan_args
from abritamr.cli.amr_status import status_args
from abritamr.cli.linelist import llist_args
from abritamr.cli.matrix import matrix_args
from abritamr.cli.infer import infer_args
from abritamr.cli.utils_database import database_args
from abritamr.cli.utils_rules import rules_args
from abritamr.cli.utils_catalog import catalog_args
from abritamr.version import __version__, db
from abritamr.utils import output_results
from abritamr.logger import log


def run_scan(args):
    """Run gene detection and write the scan results."""
    log.info(f"Running scan to identify genes and and drugclasses")
    amr = scan.scan(args)
    output_results(
        df=amr, output=args.output, _format=args.format, workdir=args.workdir
    )


def run_status(args):
    """Determine AMR status from scan results and write the report."""
    log.info(
        f"Determining the amr status of the sequence based on genes and drug classes."
    )
    amr = amr_status.amr_status(args)
    output_results(df=amr, output=args.output, _format=args.format)


def run_matrix(args):
    """Build and write a matrix from sample results."""
    log.info(f"Running matrix to collate single sample results into a matrix.")
    result = matrix.matrix(args)
    output_results(
        df=result, output=args.output, _format=args.format, workdir=args.workdir
    )


def run_linelist(args):
    """Build and write a reportable linelist from sample results."""
    log.info(
        f"Running linelist to collate single sample results into a linelist for reporting."
    )
    result = linelist.linelist(args)
    output_results(
        df=result, output=args.output, _format=args.format, workdir=args.workdir
    )


def run_gdst(args):
    """Infer and write genotypic drug susceptibility results."""
    log.info(f"Running gDST to collate single sample results into a gDST report.")
    result = infer.abritamr_gdst(args)
    output_results(
        df=result, output=args.output, _format=args.format, workdir=args.workdir
    )


def run_complete(args):
    """Run the complete abriTAMR workflow."""
    log.info(f"Running all amr modules. Please be patient this may take some time.")

    result = run.run(args)


def _update_databases(args):
    """Generate or update the local reference database."""
    log.info(f"Will generate an abritamr compatible reference gene catalog.")

    catalog = update_database.generate_database(args)


def _update_catalog(args):
    """Generate or update the reference gene and rules catalog."""
    log.info(f"Will generate an abritamr compatible reference gene and rules catalog.")

    catalog = utils_catalog.catalog(args)


def _update_rules(args):
    """Generate or update the AMR rules catalog."""
    log.info(f"Will generate an abritamr compatible rules catalog.")

    catalog = utils_rules.rules(args)


def cli():
    """Parse command-line arguments and dispatch the selected command."""
    parser = argparse.ArgumentParser(
        description=f"****abritamr - AMR gene detection and reporting pipeline - version {__version__}****",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )

    subparsers = parser.add_subparsers(help="What would you like to do?")
    # get parsers
    parser_sub_run = run_args(subparsers=subparsers)
    parser_sub_scan = scan_args(subparsers=subparsers)
    parser_sub_status = status_args(subparsers=subparsers)
    parser_sub_infer = infer_args(subparsers=subparsers)
    parser_sub_llist = llist_args(subparsers=subparsers)
    parser_sub_matrix = matrix_args(subparsers=subparsers)

    # utils modules sub sub parser

    parser_sub_utils = subparsers.add_parser(
        "utils",
        help="Utility functions for using AMR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_sub_utils_level2 = parser_sub_utils.add_subparsers(dest="level2")
    parser_sub_utils_database = database_args(subparsers=parser_sub_utils_level2)
    parser_sub_utils_rules = rules_args(subparsers=parser_sub_utils_level2)
    parser_sub_utils_catalog = catalog_args(subparsers=parser_sub_utils_level2)

    # tie parsers to functions
    parser_sub_run.set_defaults(func=run_complete)
    parser_sub_scan.set_defaults(func=run_scan)
    parser_sub_status.set_defaults(func=run_status)
    parser_sub_llist.set_defaults(func=run_linelist)
    parser_sub_infer.set_defaults(func=run_gdst)
    parser_sub_matrix.set_defaults(func=run_matrix)
    parser_sub_utils_database.set_defaults(func=_update_databases)
    parser_sub_utils_rules.set_defaults(func=_update_rules)
    parser_sub_utils_catalog.set_defaults(func=_update_catalog)
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
    # elif len(sys.argv) <= 2:
    #     parser_sub_scan.# print_help(sys.stderr)
    elif len(sys.argv) <= 2 and sys.argv[1] == "utils":
        parser_sub_utils.print_help(sys.stderr)
    else:
        args.func(args)


if __name__ == "__main__":
    cli()
