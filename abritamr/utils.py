import pathlib
import subprocess
from abritamr.version import db
import pandas as pd
import sys
import logging 
import json
from abritamr.abritamr_logging import log
# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
# handler = logging.StreamHandler(sys.stderr)

# log = logging.getLogger(__name__)
# log.addHandler(handler)
# log.setLevel(logging.DEBUG)

def check_path(pth:str) -> bool:
    log.info(f"Checking path: {pth}")
    if not pathlib.Path(pth).exists():
        return True
    return True

def get_refgenes()-> pd.DataFrame:

    rf = pathlib.Path(__file__).parent / 'db' /'refgenes_latest.csv'
    if rf.exists():
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

def wrangle_species(organism:str) -> tuple:

    try:
        with open(f"{pathlib.Path(__file__).parent / 'configs' / 'species_config.json'}") as j:
            SPCFG = json.load(j)

    except Exception as e:
        log.critical(f"Something has gone very wrong : {e}.")
        raise SystemExit

    if not organism:
        log.info("No species specific criteria will be applied")
        return ""

    else:
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
        