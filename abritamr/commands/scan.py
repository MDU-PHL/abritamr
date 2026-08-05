import click
import pathlib
import json
import logging
import pandas as pd
import sys

from abritamr.utils import output_results, check_assembly, check_amrfinder, check_any2fasta, wrangle_species, output_results, check_path
from abritamr.run_finder import run_amrf
from abritamr.parse_finder import amrf2dict
# from abritamr.parse_reportable import add_abritamr_results
from abritamr.abritamr_classes import apply_classes
# from abritamr.amr_report import summary

logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
handler = logging.StreamHandler(sys.stderr)

log = logging.getLogger(__name__)
log.addHandler(handler)

def generate_output(species:str,sample_id:str,amr:list)-> list:
    amr = apply_classes(amr = amr, species=species, sid = sample_id)
    
    return amr

def scan(
    args
        ) -> dict:

    if not args.assembly and not args.amrfinderplus :
        log.critical("You must supply an input file (assembly or amrfinder plus output). Exiting.")
        raise SystemExit(1)
    amr = []
    
    if args.assembly:
        log.info("Assembly(ies) have been supplied.")
        organism = wrangle_species(organism = args.species)
        for asm in args.assembly:  
            
            log.info(f"Will now try to run amrfinderplus on supplied assembly {asm}")
            if check_path(pth = f"{asm}") and check_assembly(f"{asm}"):
                full_path = f"{pathlib.Path(f'{asm}').absolute()}"
                sample_id = args.sample_id if args.sample_id else full_path
                log.info(f"Running amrfinder plus")
                res = run_amrf(
                    min_identity = args.min_identity, 
                    min_coverage = args.min_identity, 
                    asm=asm,
                    threads=args.threads, 
                    organism=organism
                    )
                
            # amr.extend(res)
                res = generate_output(species = args.species, sample_id = sample_id, amr = res )
                
                # print(res)
                amr.extend(res)
    if args.amrfinderplus:
        for afp in args.amrfinderplus:
            sample_id = args.sample_id if args.sample_id else full_path
            full_path = f"{pathlib.Path(f'{args.afp}').absolute()}"
            log.info(f"Opening existing amrfinder plus output")

            res = amrf2dict(amrfinder = args.amrfinder)
            res = generate_output(species = args.species, sample_id = sample_id, amr = res  )
            amr.append(res)
        
    abritamr_cols = [
        'sample_id',
        'Element symbol',
        'abritamr_class', 
        'abritamr_subclass', 
        # 'abritamr_AMR_reporting_status',
        # 'abritamr_AMR_type', 
        # 'criteria_id', 
        # 'criteria_version',
        'Element name', 
        'Scope', 
        'Type', 
        'Subtype',
        'Method',
        'pmid', 
        'abritamr_class_version',
        'Target length', 
        'Reference sequence length',
        '% Coverage of reference', 
        '% Identity to reference',
        'Alignment length', 
        'Closest reference accession',
        'Closest reference name', 
        'HMM accession', 
        'HMM description',
        ]
    amr = pd.DataFrame(amr)
    amr = amr[abritamr_cols]
    

    # print(pd.DataFrame(amr))
    
    # output_results(amr = amr, output = kwargs['output'])

    return amr