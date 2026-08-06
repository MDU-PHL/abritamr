import subprocess
import pandas as pd
from abritamr.utils import check_amrfinder, check_any2fasta, check_assembly

import logging
from abritamr.logger import log

# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
# log = logging.getLogger(__name__)
# log.setLevel(logging.DEBUG)

# pd.DataFrame(rows[1:], columns = rows[0]).to_dict(orient= "records")


def generate_cmd(min_identity:float, min_coverage:float, asm:str,threads:int, organism:str) -> str:
    spc = ""
    if organism != "":
        spc= f"-O {organism}"
    
    cmd = f"amrfinder -n {asm} --plus --ident_min {min_identity} --coverage_min {min_coverage} --threads {threads} {spc}"

    return cmd

def run_cmd(cmd:str) -> str:
    log.info(f"Running amrfinder: {cmd}")
    proc = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
    # print(proc)
    if proc.returncode != 0:
        err = proc.stderr.split("\n")
        log.critical(f"The following error was reported:")
        log.critical("\n".join(err))

    else:
        log.info("AMRfinder plus was successful.")
        return proc.stdout

def parse_output(results:str) -> dict:

    rows = results.split('\n')
    rows = [row.split('\t') for row in rows if row != ""]
    return pd.DataFrame(rows[1:], columns = rows[0]).to_dict(orient= "records")

def run_amrf(min_identity:float, min_coverage:float, asm:str,threads:int, organism:str) -> dict:

    cmd = generate_cmd(min_identity = min_identity, min_coverage = min_coverage, asm = asm, threads = threads, organism = organism)
    stdout = run_cmd(cmd = cmd) 
    amr = parse_output(results = stdout)

    return amr