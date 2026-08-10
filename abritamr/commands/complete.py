import click
import pathlib
import json
import logging
import pandas as pd
import sys

from abritamr.utils import output_results, check_assembly, check_amrfinder, check_any2fasta, wrangle_species, output_results, check_path,abritamr_scan_columns,abritamr_status_columns
from abritamr.run_finder import run_amrf
from abritamr.parse_finder import amrf2dict
# from abritamr.parse_reportable import add_abritamr_results
from abritamr.commands.scan import scan
from abritamr.drugclasses import apply_classes
from abritamr.logger import log
from abritamr.parse_reportable import add_abritamr_results
from abritamr.amr_report import summary
from abritamr.amr_matrix import summary as matrix_summary

def save_output(workdir:str,sample_id:str,result:pd.DataFrame, outname:str, _format:str="csv", no_keep:bool = False)-> bool:
    
    if not no_keep:
        dlm= ","
        if _format == "tab":
            dlm = "\t"
            _format = "txt"

        log.info(f"Creating output directory for {sample_id}")
        fldr = pathlib.Path(workdir, sample_id)
        try:
            fldr.mkdir(exist_ok=True)

            log.info(f"Creating output file {outname}")

            result.to_csv(f"{fldr}/{outname}.{_format}", sep = dlm, index = False)

            return True

        except Exception as e:
            log.critical(f"Something has gone wrong saving {outname}. The following error was encountered : {e}")
            raise SystemExit(1)

    return True


def check_multi(pth: str, workdir:pathlib.Path) -> bool:

    if check_path(f"{workdir}"):

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
            # df = df.rename(columns = {'assembly': 'input', 'amrfinder':'input'})
            df = df.fillna("")
            lines= []
            sids = []
            for row in df.iterrows():
                d = {}
                if row[1]['sample_id'] in sids:
                    log.critical(f"It looks like you have duplicate IDs. Sample IDs must be unique when running complete workflow. abritamr will not run. Please try again.")
                    raise SystemExit(1)
                sids.append(row[1]['sample_id'])
                d['sample_id'] = row[1]['sample_id']
                for ip in ['amrfinder','assembly']:
                    if ip in df.columns.tolist():
                        if not check_path(row[1][f"{ip}"]):
                            log.warning(f"{row[1][f'{ip}']} is not a valid file path. abritamr will not run on this input.")
                        else:
                            d[ip] = row[1][f"{ip}"]
                
                
                # print(sids)
                orgs = wrangle_species(row[1]['species'])
                d['species'] = species
                
                lines.append(d)
                
        
            return lines
            
        else:
            log.critical(f"{pth} does not exist. Please provide a valid path to the multi input file.")
            raise SystemExit(1)
    else:
        log.critical(f"{workdir} is not a valid path. Please provide a valid working directory.")
        raise SystemExit(1)
                # print(line)
    

def complete(
    args
        ) -> dict:
    
    simple = True if args.viewtype == 'compact' else False
    # dbv = "unknown"
    if not args.assembly and not args.amrfinderplus and not args.multi:
        log.critical("You must supply an input file (input file with multiple samples OR as single assembly or single amrfinder plus output). Exiting.")
        raise SystemExit(1)
    if not args.multi:
        if not args.sample_id: 
            log.critical(f"You must supply a sample_id column. Please check your inputs and try again.")
            raise SystemExit(1)
        inputs = [
            {
            'sample_id' : args.sample_id,
            'assembly' : args.assembly if args.assembly else "",
            'amrfinder': args.amrfinderplus if args.amrfinderplus else "",
            'species': args.species

        }
        ]
    else:
        inputs = check_multi(pth=args.multi, workdir = pathlib.Path(args.workdir))
    linelists = []
    matrices = []
    for i in inputs:
        amr = []
        dbv = "unknown"
    
        if i['assembly'] != "":
            if check_path(pth = f"{asm}") and check_assembly(f"{asm}"):
                dbv = check_amrfinder()
                res = run_amrf(
                    min_identity = args.min_identity, 
                    min_coverage = args.min_identity, 
                    asm=i['assembly'],
                    threads=args.threads, 
                    organism=i['species']
                    )
            
        elif i['amrfinder'] != "":
            res = amrf2dict(amrfinder = args.amrfinder)           
        

        # res = generate_output(species =  amr = res )
        amr = apply_classes(amr = res, species =  i['species'], sid = i['sample_id'],)
        # amr.extend(res)
        
        scanned = pd.DataFrame(amr)
        scanned['amrfinderplus_db_version'] = dbv
        scanned_cols = abritamr_scan_columns()
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = scanned[scanned_cols], outname = "abritamr_scan", _format=args.format, no_keep = args.no_keep)
        amr = add_abritamr_results(amr = amr,cfgpath=args.cfgpath)
        typed = pd.DataFrame(amr)
        typed_cols = abritamr_status_columns()
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = typed[typed_cols], outname = "abritamr_typed", _format=args.format, no_keep = args.no_keep)
        linelist = summary(
            results =typed_cols, 
            _format = args.format, 
            sid = sample_id,
            simple = simple, 
            genesonly = args.genesonly,
            minidentity = args.min_identity,
            mincoverage = args.min_coverage,
            )
        linelists.append(linelist)
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = linelist, outname = "abritamr_linelist", _format=args.format, no_keep = False)
        matrix = matrix_summary(
            results =typed, 
            facet =  args.facet,
            sid = i['sample_id'],
            minidentity = args.min_identity,
            mincoverage = args.min_coverage,
            )
            # print(res)
        matrices.append(matrix)
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = matrix, outname = "abritamr_matrix", _format=args.format, no_keep = False)
        # print(amr)
    dlm= ","
    if _format == "tab":
        dlm = "\t"
        _format = "txt"
    log.info(f"Saving a single file for output of linelist.")
    pd.concat(linelists).to_csv(f"{args.workdir}/{args.prefix}_linelist_report.{_format}", sep = {dlm}, index = False)
    log.info(f"Saving a single for output of matrix.")
    pd.concat(matrices).to_csv(f"{args.workdir}/{args.prefix}_matrix.{_format}", sep = {dlm}, index = False)

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