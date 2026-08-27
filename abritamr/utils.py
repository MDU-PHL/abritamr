from datetime import datetime
import pathlib
import subprocess
from abritamr.version import db
import pandas as pd
import sys
import logging 
import json
from abritamr.logger import log
from abritamr.run_sourmash import run_sourmash_search
# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
# handler = logging.StreamHandler(sys.stderr)

# log = logging.getLogger(__name__)
# log.addHandler(handler)
# log.setLevel(logging.DEBUG)

def _get_date():

    return datetime.today().strftime('%Y-%m-%d')


def check_path(pth:str) -> bool:
    log.info(f"Checking path: {pth}")
    if not pathlib.Path(pth).exists():
        return False
    return True

# this needs to be user defined!!
def get_refgenes(pth)-> pd.DataFrame:

    if check_path(pth = pth):
        rf = pathlib.Path(pth)
        return pd.read_csv(rf)
    else:
        log.critical(f"Something is very wrong - the reference catalog is not present.")
        raise SystemExit(1)


def check_any2fasta()-> bool:
    log.info(f'Checking that any2fasta is installed.')
    proc = subprocess.run("any2fasta -v", shell = True, capture_output = True, encoding = "utf-8")
    if proc.returncode == 0:
        log.info("any2fasta is installed.")
        return True
    else:
        log.critical(f"any2fasta is a dependency of abritamr. Please check your installation.")
        raise SystemExit(1)

def abritamr_scan_columns() -> list:
    return [
        'sample_id',
        'Element symbol',
        'abritamr_class', 
        'abritamr_subclass', 
        'abritamr_mechanism',
        # 'abritamr_priority_status',
        # 'abritamr_AMR_type', 
        # 'criteria_id', 
        # 'criteria_version',
        'species',
        'Element name', 
        'Scope', 
        'Type', 
        'Subtype',
        'Method',
        'pmid', 
        'abritamr_class_version',
        'amrfinderplus_db_version',
        'Target length', 
        'Reference sequence length',
        '% Coverage of reference', 
        '% Identity to reference',
        'Alignment length', 
        'Closest reference accession',
        'Closest reference name', 
        'HMM accession', 
        'HMM description',
        'abritamr_accession_key',
        'amrrules_mutation'
        ]

def abritamr_status_columns() -> list:
    return [
        'sample_id',
        'Element symbol',
        'abritamr_class', 
        'abritamr_subclass', 
        'abritamr_mechanism',
        'abritamr_priority_status',
        'abritamr_AMR_type', 
        'species',
        'criteria_id', 
        'criteria_version',
        'Element name', 
        'Scope', 
        'Type', 
        'Subtype',
        'Method',
        'pmid', 
        'abritamr_class_version',
        'amrfinderplus_db_version',
        'Target length', 
        'Reference sequence length',
        '% Coverage of reference', 
        '% Identity to reference',
        'Alignment length', 
        'Closest reference accession',
        'Closest reference name', 
        'HMM accession', 
        'HMM description',
        'abritamr_accession_key',
        'amrrules_mutation'
        ]

def check_sourmash() -> bool:
    log.info(f"Checking that sourmash is installed.")
    proc = subprocess.run("sourmash -v", shell = True, capture_output = True, encoding = "utf-8")
    if proc.returncode == 0:
        log.info("sourmash is installed.")
        return True
    else:
        log.critical(f"sourmash is a dependency of abritamr. Please check your installation.")
        raise SystemExit(1)

def check_amrfinder()-> bool:
    log.info(f"Checking that amrfinder plus is installed and database versions are compatible.")
    vrsn = subprocess.run("amrfinder -V", shell = True, capture_output = True, encoding = "utf-8")
    if vrsn.returncode != 0:
        log.critical(f"Something is wrong with your environment. amrfinderplus is required for abritamr to run. Please follow installation instructions and try again.")
    else:
        dbv = [i for i in vrsn.stdout.split('\n') if 'Database version' in i]
        # print(dbv)
        dbv = dbv[0].split(':')[-1].strip() if dbv != [] else ""
        if dbv == "":
            log.critical(f"It looks like the database version cannont be determined. There may be something wrong with your installation. Please check documentation and try again.")
            raise SystemExit(1)
        if db not in dbv:
            log.warning(f"Your amrfinder database is at {dbv}. This is different to what is expected ({db}). Please not unexpected behaviour may occur.")
            return dbv
        else:
            log.info(f"Your amrfinder installation is compatible with current abritamr.")
            return dbv

def check_assembly(pth) -> bool:
    log.info(f"Checking assembly is in a correct format")
    if check_any2fasta():
        # try:
            # print(f"any2fasta {pth}")
            proc = subprocess.run(f"any2fasta {pth}", shell = True, capture_output = True, encoding = "utf-8")
            if proc.returncode == 0:
                log.info(f"{pth} is a valid assembly file. Will no proceed with running amrfinder.")
                return True
            else:
                log.critical(f"Your assembly file appears to not be in the correct format. Please try again.")
                raise SystemExit(1)
        # except Exception as e:
        #     log.critical(f"Something has gone wrong with any2fasta checking your assembly. The following error was reported: {e}")
        #     raise SystemExit(1)

def guess_species(asm:str, sid:str="abritamr") -> str:
    """
    Guess the species of a given assembly using sourmash.

    Parameters
    ----------
    asm : str
        Path to the assembly file (e.g., FASTA or FASTQ).
    sid : str, optional
        Sample ID for the query signature. Default is "abritamr".

    Returns
    -------
    str
        The guessed species name.
    """
    sp = run_sourmash_search(query_filename = asm, SBT_filename = f"{pathlib.Path(__file__).parent / 'species_db' / 'abritamrdb.sbt.zip'}", sid = sid)
    return sp

def wrangle_species(organism:str, asm:str="", sid:str="abritamr", check_species:bool = True) -> tuple:

    try:
        with open(f"{pathlib.Path(__file__).parent / 'configs' / 'amrfinder_species.json'}") as j:
            SPCFG = json.load(j)

    except Exception as e:
        log.critical(f"Something has gone very wrong : {e}.")
        raise SystemExit

    if not organism and check_species and asm != "":
        log.info("No species supplied - will use sourmash to try to guess the best match for AMR classification.")

        organism = guess_species(asm, sid=sid)


    if organism != "":
        og = '_'.join(organism.split())
        if og in SPCFG:
            return og
        elif og[0] in SPCFG:
            return og[0]
        else:
            return ""

def output_results(df:pd.DataFrame, output:str= "",_format:str="csv") -> bool:


    dlm = "," if _format == 'csv' else '\t'
    suf = 'csv' if _format == 'csv' else 'txt'

    if not output:
        df.to_csv(sys.stdout, sep = dlm,index = False)
    else:
        df.to_csv(f"{output}.{suf}", sep = dlm, index = False)
        