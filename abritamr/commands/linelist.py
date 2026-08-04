import click
import pathlib
import json
import logging
import pandas as pd

from abritamr.utils import check_assembly, check_amrfinder, check_any2fasta, wrangle_species, output_results
from abritamr.run_finder import run_amrf
from abritamr.parse_finder import amrf2dict
from abritamr.parse_reportable import add_abritamr_results
from abritamr.parse_amrtype import get_amr_type
from abritamr.amr_report import summary


logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

try:
    with open(f"{pathlib.Path(__file__).parent / 'species_config.json'}") as j:
        SPCFG = json.load(j)

except Exception as e:
    log.critical(f"Something has gone very wrong : {e}.")
    raise SystemExit

@click.command()
@click.option(
    '--assembly',
    '-asm',
    help = "Assembly file to use as input (*.fa*, *.gbk *.fa*.gz, *.gbk.gz)",
    default = "",
    show_default = True
)
@click.option(
    '--amrfinder',
    '-amf',
    help = "EXPERIMENTAL - USE WITH CAUTION. AMR finder output file. Please note unexpected behaviour may result were versions and databases differ from abritamr",
    default = "",
    show_default = True
)
@click.option(
    '--sample-id',
    '-s',
    help = "sample identifier, this will be used to name output files and in line list reports",
    default = "abritamr",
    show_default = True
)
@click.option(
    '--format',
    '-f',
    help = "Output format",
    type = click.Choice(["csv", "tab"]),
    default = "csv",
    show_default = True
)
@click.option(
    '--output',
    '-o',
    help = "Filename to save output - default stdout.",
    default = "",
    show_default = True
)
@click.option(
    '--summary/--no-summary',
    is_flag=True,
    default=False,
    show_default=True,
    help = "Set --summary if you require summarised report output. Report will be saved as sample-id_report.csv"
)
@click.option(
    '--species',
    '-sp',
    help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence.",
    default = "",
    show_default = True,
    type = click.Choice(SPCFG)
)
@click.option(
    '--threads',
    '--cpus',
    help="Number of max CPU cores to run.",
    default=1,
    show_default=True
)
@click.option(
    '--viewtype',
    '-vt',
    default = "full",
    show_default=True,
    type = click.Choice(["full", "compact"]),
    help = "Format of abritamr report. Only applicable if --summary is set. Default - full will output all results for the sequence. Compact will only output summarised results."
)
@click.option(
    '--genesonly',
    is_flag=True,
    default=False,
    show_default=True,
    help = "If only whole genes, not SNPs should be reported. Only applicable if --summary is set."
)
@click.option(
    "--min-identity",
    default = 0.9,
    help ="Minimum identity to reference gene for reporting a match.",
    show_default = True
)
@click.option(
    "--min-coverage",
    default = 0.5,
    help ="Minimum coverage of the reference gene for reporting a complete match.",
    show_default = True
)
@click.option(
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help = "Set to override existing outputs of the same name."
)
@click.option(
    "--reportable-config",
    "-rc",
    default = f"{pathlib.Path(__file__).parent / 'configs' / 'reportable_criteria.json'}",
    help = "Path to config that defines the criteria for highlighting relevant genes and mechanisms for reporting. Supplying a new config will override the default abritamr criteria.",
    show_default = True
)
@click.option(
    "--amr_type-config",
    "-ac",
    default = f"{pathlib.Path(__file__).parent / 'configs' / 'amr_type_criteria.json'}",
    help = "Path to config that defines the criteria for amr type based on genes and mechanisms detected. Supplying a new config will override the default abritamr criteria.",
    show_default = True
)
def scan(
    **kwargs
        ) -> dict:

    
    # print(kwargs)
   
    if kwargs['assembly'] == "" and kwargs['amrfinder'] == "":
        log.critical("You must supply an input file (assembly or amrfinder plus output). Exiting.")
        raise SystemExit(1)
    if kwargs['assembly'] != "" and pathlib.Path(f"{kwargs['assembly']}").exists() and check_assembly(f"{kwargs['assembly']}"):
        log.info(f"Running amrfinder plus")
        amr = run_amrf(
            min_identity = kwargs['min_identity'], 
            min_coverage = kwargs['min_identity'], 
            asm=kwargs['assembly'],
            threads=kwargs['threads'], 
            organism=kwargs['species']
            )
    elif kwargs['amrfinder'] != "":
        log.info(f"Opening existing amrfinder plus output")
        amr = amrf2dict(amrfinder = amrfinder)
    

    species,genus = wrangle_species(organism = kwargs['species'])
    
    amr = add_abritamr_results(amr = amr, cfgpath = kwargs['reportable_config'], species=species, genus = genus, sid = kwargs['sample_id'])
    amr= get_amr_type(amr=amr, species=species, genus=genus, cfgpath=kwargs['amr_type_config'])

    output_results(amr = amr, output = kwargs['output'])

    return amr