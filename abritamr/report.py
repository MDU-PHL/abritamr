import click
import pathlib


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
    '--sample_id',
    '-sid',
    help = "sample identifier, this will be used to name output files and in line list reports",
    default = "",
    show_default = True
)
@click.option(
    '--species',
    '-s',
    help = "Species from which assemblies were derived. Must be supplied for SNP detection and inferrence.",
    default = "",
    show_default = True,
    type = click.Choice(["Acinetobacter_baumannii", "Bordetella_pertussis", "Burkholderia_cepacia", "Burkholderia_mallei", "Burkholderia_pseudomallei", "Campylobacter", "Citrobacter_freundii", "Clostridioides_difficile", "Corynebacterium_diphtheriae", "Enterobacter_asburiae", "Enterobacter_cloacae", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Haemophilus_influenzae", "Helicobacter_pylori", "Klebsiella_oxytoca", "Klebsiella_pneumoniae", "Neisseria_gonorrhoeae", "Neisseria_meningitidis", "Pseudomonas_aeruginosa", "Salmonella", "Serratia_marcescens", "Staphylococcus_aureus", "Staphylococcus_epidermidis", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae", "Vibrio_parahaemolyticus", "Vibrio_vulnificus"])
)
@click.option(
    '--infer',
    is_flag = True,
    show_default=True,
    help = "Set for generating an inferred antibiogram. Note by default only 'Salmonella' is supported at the moment."
)
@click.option(
    '--threads',
    '--cpus',
    help="Number of max CPU cores to run.",
    default=1
)
@click.option(
    '--summary/--no-summary',
    is_flag=True,
    default=True,
    show_default=True,
    help = "Set --no-summary if you don't require summarised report output."
)
@click.option(
    '--viewtype',
    '-vt',
    default = "full",
    show_default=True,
    type = click.Choice(["full", "compact"]),
    help = "Format of output of abritamr. Default - full will output all results for the sequence. Compact will only output summarised results."
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
    "--reportable-config",
    "-rc",
    default = f"{pathlib.Path(__file__).parent / 'configs' / 'reportable_criteria.json'}",
    help = "Path to config that defines the criteria for highlighting relevant genes and mechanisms for reporting. Supplying a new config will override the default abritamr criteria.",
    show_default = True
)
@click.option(
    "--inferrence-config",
    "-ic",
    default = f"{pathlib.Path(__file__).parent / 'configs' / 'infer_criteria.json'}",
    help = "Path to config that defines the criteria for inferring genomic AST from detected genes and mechanisms. Supplying a new config will override the default abritamr criteria.",
    show_default = True
)
@click.option(
    "--amr_type-config",
    "-ac",
    default = f"{pathlib.Path(__file__).parent / 'configs' / 'amr_type_criteria.json'}",
    help = "Path to config that defines the criteria for amr type based on genes and mechanisms detected. Supplying a new config will override the default abritamr criteria.",
    show_default = True
)
def run(
    assembly:str,
    amrfinder:str,
    species:str,
    infer:bool,
    threads:int,
    raw:str,
    viewtype:str,
    output:str,
    min_identity:int,
    min_coverage:int,
    reportable_config:str,
    inferrence_config:str,
    amr_type_config:str
        ) -> bool:

    print('Will run abritamr')

