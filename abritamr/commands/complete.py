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
from abritamr.commands.scan import scan
from abritamr.drugclasses import apply_classes
from abritamr.logger import log

def generate_output(species:str,sample_id:str,amr:list)-> list:
    amr = apply_classes(amr = amr, species=species, sid = sample_id)
    
    return amr

def check_input(pth: str) -> bool:

    if check_path(pth):
        df = pd.read_csv(pth, sep = "\t")
        if 'species' not in df.columns:
            df['species'] = ""
        if 'sample_id' not in df.columns:
            log.critical(f"You must supply a sample_id column. Please check your inputs and try again.")
            raise SystemExit(1)
        if 'assembly' not in df.columns and 'amrfinder' not in df.columns:
            log.critical(f"You must supply a path to either assemblies or amrfinder output. Please try again.")
            raise SystemExit(1)
        if 'assembly' in df.columns and 'amrfinder' in df.columns:
            log.critical(f"You must supply a path to either assemblies OR amrfinder output. Please try again.")
            raise SystemExit(1)
        df = df.rename(columns = {'assembly': 'input', 'amrfinder':'input'})
        lines_sp = []
        lines_nosp = []
        sids = []
        for row in df.iterrows():
            if not check_path(row[1]['input']):
                log.warning(f"{row[1]['input']} is not a valid file path. abritamr will not run on this input.")
            else:
                if row[1]['sample_id'] in sids:
                    log.critical(f"It looks like you have duplicate IDs. Sample IDs must be unique when running complete workflow. abritamr will not run so you can remove these. Please try again.")
                    raise SystemExit(1)
                sids.append(row[1]['sample_id'])
                orgs = wrangle_species(row[1]['species'])
                line = [row[1]['sample_id'],row[1]['input'], orgs]
    else:
        log.critical(f"There is something wrong with your input file. Please check your inputs and try again.")
        raise SystemExit(1)

def complete(
    args
        ) -> dict:
    kargs = args
    if not kargs.assembly and not kargs.amrfinderplus and not kargs.multi:
        log.critical("You must supply an input file (input file with multiple samples OR as single assembly or single amrfinder plus output). Exiting.")
        raise SystemExit(1)
    if not kargs.multi:
        if args.assembly:
            kargs.assembly = [args.assembly]
        if args.amrfinderplus:
            kargs.amfinderplus = [args.amrfinderplus]
        amr = scan(kargs)
    else:

        pass
        # print(amr)
    # amr = []
    # dbv = "unknown"
    # if args.assembly:
    #     log.info("Assembly(ies) have been supplied.")
    #     organism = wrangle_species(organism = args.species)
    #     for asm in args.assembly:  
            
    #         log.info(f"Will now try to run amrfinderplus on supplied assembly {asm}")
    #         if check_path(pth = f"{asm}") and check_assembly(f"{asm}"):
    #             dbv = check_amrfinder()

    #             full_path = f"{pathlib.Path(f'{asm}').absolute()}"
    #             sample_id = args.sample_id if args.sample_id else full_path
    #             log.info(f"Running amrfinder plus")
    #             res = run_amrf(
    #                 min_identity = args.min_identity, 
    #                 min_coverage = args.min_identity, 
    #                 asm=asm,
    #                 threads=args.threads, 
    #                 organism=organism
    #                 )
                
    #         # amr.extend(res)
    #             res = generate_output(species = args.species, sample_id = sample_id, amr = res )
                
    #             # print(res)
    #             amr.extend(res)
    # if args.amrfinderplus:
    #     for afp in args.amrfinderplus:
    #         sample_id = args.sample_id if args.sample_id else full_path
    #         full_path = f"{pathlib.Path(f'{args.afp}').absolute()}"
    #         log.info(f"Opening existing amrfinder plus output")

    #         res = amrf2dict(amrfinder = args.amrfinder)
    #         res = generate_output(species = args.species, sample_id = sample_id, amr = res  )
    #         amr.append(res)
        
    # abritamr_cols = [
    #     'sample_id',
    #     'Element symbol',
    #     'abritamr_class', 
    #     'abritamr_subclass', 
    #     # 'abritamr_AMR_reporting_status',
    #     # 'abritamr_AMR_type', 
    #     # 'criteria_id', 
    #     # 'criteria_version',
    #     'species',
    #     'Element name', 
    #     'Scope', 
    #     'Type', 
    #     'Subtype',
    #     'Method',
    #     'pmid', 
    #     'abritamr_class_version',
    #     'amrfinderplus_db_version',
    #     'Target length', 
    #     'Reference sequence length',
    #     '% Coverage of reference', 
    #     '% Identity to reference',
    #     'Alignment length', 
    #     'Closest reference accession',
    #     'Closest reference name', 
    #     'HMM accession', 
    #     'HMM description',
    #     ]
    # amr = pd.DataFrame(amr)
    # amr['amrfinderplus_db_version'] = dbv
    # amr = amr[abritamr_cols]
    

    # # print(pd.DataFrame(amr))
    
    # # output_results(amr = amr, output = kwargs['output'])

    # return amr