
import click
import pathlib
import json
import logging
import pandas as pd
import sys

from abritamr.utils import output_results, check_assembly, check_amrfinder, check_any2fasta, guess_species, wrangle_species, output_results, check_path,abritamr_scan_columns,abritamr_status_columns
from abritamr.run_finder import run_amrf
from abritamr.parse_finder import amrf2dict
# from abritamr.parse_reportable import add_abritamr_results
from abritamr.commands.scan import scan
from abritamr.drugclasses import apply_classes
from abritamr.logger import log
from abritamr.parse_reportable import add_abritamr_results
from abritamr.amr_report import summary
from abritamr.amr_matrix import summary as matrix_summary
from abritamr.amr_infer import gdst_results_to_df_long, gdst_results_to_df_wide,gdst

def save_output(workdir:str,sample_id:str,result:pd.DataFrame, outname:str, _format:str="csv", no_keep:bool = False)-> bool:
    
    if not no_keep:
        dlm= ","
        if _format == "tab":
            dlm = "\t"
            _format = "txt"

        if sample_id != "":

            log.info(f"Creating output directory for {sample_id}")
            fldr = pathlib.Path(workdir, sample_id)
        else:
            fldr = pathlib.Path(workdir)
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
                
                
                if row[1]['species'] != "":
                    guess = guess_species(asm = row[1]['assembly'], sid = row[1]['sample_id']) if row[1]['assembly'] != "" else ""
                else:
                    guess = row[1]['species']
                
                d['species'] = guess
                
                lines.append(d)
                
        
            return lines
            
        else:
            log.critical(f"{pth} does not exist. Please provide a valid path to the multi input file.")
            raise SystemExit(1)
    else:
        log.critical(f"{workdir} is not a valid path. Please provide a valid working directory.")
        raise SystemExit(1)
                # # print(line)



def run(
    args
        ) -> dict:
    
    simple = True if args.viewtype == 'compact' else False
    # dbv = "unknown"
    if not args.contigs and not args.amrfinderplus and not args.multi:
        log.critical("You must supply an input file (input file with multiple samples OR as single assembly or single amrfinder plus output). Exiting.")
        raise SystemExit(1)
    if not args.multi:
        if not args.sample_id: 
            log.critical(f"You must supply a sample_id column. Please check your inputs and try again.")
            raise SystemExit(1)
        species = guess_species(asm = args.contigs[0], sid = args.sample_id) if args.contigs else ""
        inputs = [
            {
            'sample_id' : args.sample_id,
            'assembly' : args.contigs if args.contigs else "",
            'amrfinder': args.amrfinderplus if args.amrfinderplus else "",
            'species': species

        }
        ]
    else:
        inputs = check_multi(pth=args.multi, workdir = pathlib.Path(args.workdir))
    # # print(inputs)
    linelists = []
    matrices = []
    infers = []
    for i in inputs:
        amr = []
        dbv = "unknown"
        # # print(i)
        res = []
        
        if 'assembly' in i and i['assembly'] != "":
            if check_path(pth = f"{i['assembly'][0]}") and check_assembly(f"{i['assembly'][0]}"):
                dbv = check_amrfinder()
                log.info(f"Running amrfinder plus on supplied assembly {i['assembly'][0]} with {args.threads} threads.")
                res = run_amrf(
                    min_identity = args.min_identity, 
                    min_coverage = args.min_coverage,
                    asm=i['assembly'][0],
                    threads=args.threads, 
                    organism=i['species']
                    )

        elif 'amrfinder' in i and i['amrfinder'] != "":
            res = amrf2dict(amrfinder = i['amrfinder'])
        for r in res:
            r['amrfinderplus_db_version'] = dbv
        # # print(i)
        # res = generate_output(species = i['species'], amr = res )
        amr = apply_classes(amr = res, species =  i['species'], sid = i['sample_id'],catalog = args.reference_catalog)
    #     # amr.extend(res)
        # # print(amr)
        scanned = pd.DataFrame(amr)
        # # print(scanned.columns.tolist())
    #     # scanned['amrfinderplus_db_version'] = dbv
        scanned_cols = abritamr_scan_columns()
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = scanned[scanned_cols], outname = "abritamr_scan", _format=args.format, no_keep = args.no_keep)
        amr = add_abritamr_results(amr = amr,catalog = args.reference_catalog)
        # # print(amr)
        typed = pd.DataFrame(amr)
        typed_cols = abritamr_status_columns()
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = typed[typed_cols], outname = "abritamr_typed", _format=args.format, no_keep = args.no_keep)
        linelist = summary(
            results =typed, 
            _format = args.format, 
            sid = i['sample_id'],
            simple = simple, 
            genesonly = args.genesonly,
            minidentity = args.min_identity,
            mincoverage = args.min_coverage,
            )
        linelists.append(linelist)
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = linelist, outname = "abritamr_linelist", _format=args.format, no_keep = False)
        infer = gdst(
            results =scanned, 
            species = species,
            reference_folder = args.reference_folder,
            dflt_result = args.dflt_result
            )
        if infer != []:
            if args.reporttype == 'long':
                gdstlinelist = gdst_results_to_df_long(infer)
            elif args.reporttype == 'wide':
                gdstlinelist = gdst_results_to_df_wide(infer)
            infers.extend(gdstlinelist)
        # # print(linelist.columns.tolist())
            save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = gdstlinelist, outname = "abritamr_gdst", _format=args.format, no_keep = False)
        matrix = matrix_summary(
            results =typed, 
            facet =  args.facet,
            sid = i['sample_id'],
            minidentity = args.min_identity,
            mincoverage = args.min_coverage,
            refgenes = args.reference_catalog
            )
            # # print(res)
        matrices.append(matrix)
        save_output(workdir=f"{args.workdir}",sample_id=i['sample_id'],result = matrix, outname = "abritamr_matrix", _format=args.format, no_keep = False)
        # # print(amr)
    # # print(pd.concat(linelists, axis = 0, ignore_index=True))
    log.info(f"Saving a single file for output of linelist.")
    save_output(workdir=f"{args.workdir}",sample_id="",result = pd.concat(linelists).reset_index(drop=True), outname = f"{args.prefix}_linelist", _format=args.format, no_keep = False)
    
    log.info(f"Saving a single for output of matrix.")
    save_output(workdir=f"{args.workdir}",sample_id="",result = pd.concat(matrices).reset_index(drop=True), outname = f"{args.prefix}_matrix", _format=args.format, no_keep = False)

    
    return amr