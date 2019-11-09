import pathlib, argparse, logging, sys, os

from abritamr.AmrSetup import Setupamr


"""
mdu_amr is designed to implement AMRFinder and parse the results compatible for MDU use. It may also be used for other purposes where the format of output is compatible

Input types:
1). Assemblies - for MDU purposes those generated by the MDU QC pipeline (shovill with spades)

"""


def run_pipeline(args):

    P = Setupamr(args)
    return P.run_amr()


# abstract out toml functions into a separate class


def set_parsers():
    parser = argparse.ArgumentParser(
        description="MDU AMR", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--mduqc",
        "-m",
        action="store_true",
        help="Set if running on MDU QC data. If set please provide the MDU QC .csv as further input and an additional output suitable for lims input will be provided.",
    )
    parser.add_argument(
        "--Singularity",
        "-S",
        action="store_true",
        help="If using singularity container for AMRfinderplus"
    )
    parser.add_argument(
        "--singularity_path",
        "-s",
        default="",
        help="Path to the singularity container for AMRfinderplus, if empty will defualt to shub://phgenomics-singularity/amrfinderplus"
    )
    parser.add_argument(
        "--contigs",
        "-c",
        default="",
        help="Tab-delimited file with sample ID as column 1 and path to assemblies as column 2",
    )
    parser.add_argument(
        "--amrfinder_output",
        "-o",
        default="",
        help="Tab-delimited file with sample ID as column 1 and path to amr_finder output files as column 2",
    )
    parser.add_argument(
        "--workdir",
        "-w",
        default=f"{pathlib.Path.cwd().absolute()}",
        help="Working directory, default is current directory",
    )
    parser.add_argument(
        "--resources",
        "-r",
        default=f"{pathlib.Path(__file__).parent }",
        help="Directory where templates are stored",
    )
    parser.add_argument(
        "--drug_classes",
        "-d",
        default=f"{pathlib.Path(__file__).parent / 'db' / 'refgenes.csv'}",
        help="Path to file ",
    )
    parser.add_argument(
        "--jobs", "-j", default=16, help="Number of AMR finder jobs to run in parallel."
    )
    parser.add_argument(
        "--keep",
        "-k",
        action="store_true",
        help="If you would like to keep intermediate files and snakemake log. Default is to remove",
    )
    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    return args


def main():
    """
    run pipeline
    """

    args = set_parsers()
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        # TODO see this page https://realpython.com/python-logging/ for logging
        # create logger with 'spam_application'
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        # create file handler which logs even debug messages
        fh = logging.FileHandler("qc.log")
        fh.setLevel(logging.WARNING)
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create formatter and add it to the handlers
        formatter = logging.Formatter(
            "[%(levelname)s:%(asctime)s] %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p"
        )
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(ch)
        logger.addHandler(fh)
        args.func(args)


if __name__ == "__main__":
    main()
