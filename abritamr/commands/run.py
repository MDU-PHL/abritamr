"""Orchestrate complete AMR workflows and save per-sample outputs."""

import pathlib

import pandas as pd

from collections import namedtuple

from abritamr.amr_infer import gdst, gdst_results_to_df_long, gdst_results_to_df_wide
from abritamr.amr_matrix import summary as matrix_summary
from abritamr.amr_report import summary

from abritamr.drugclasses import apply_classes
from abritamr.logger import log
from abritamr.parse_finder import amrf2dict
from abritamr.parse_reportable import add_abritamr_results
from abritamr.run_finder import run_amrf
from abritamr.utils import (
    check_path,
    guess_species,
)
from abritamr.commands.scan import run_scan, generate_inputs
from abritamr.commands.amr_status import do_typing
from abritamr.commands.linelist import generate_linelist
from abritamr.commands.matrix import make_matrix
from abritamr.commands.infer import do_gdst


def save_output(
    workdir: str,
    sample_id: str,
    result: pd.DataFrame,
    outname: str,
    _format: str = "csv",
) -> bool:
    """Save a result dataframe under the selected sample and output name."""
    dlm = ","
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

        result.to_csv(f"{fldr}/{outname}.{_format}", sep=dlm, index=False)

        return True

    except Exception as e:
        log.critical(
            f"Something has gone wrong saving {outname}. The following error was encountered : {e}"
        )
        raise SystemExit(1)

    return True


def check_multi(pth: str, workdir: pathlib.Path) -> bool:
    """Validate a multi-sample input table and build per-sample settings."""
    if check_path(f"{workdir}"):
        if check_path(pth):
            df = pd.read_csv(pth, sep="\t")
            if "species" not in df.columns:
                df["species"] = ""
            if "sample_id" not in df.columns:
                log.critical(
                    f"You must supply a sample_id column when using multi mode. Please check your inputs and try again."
                )
                raise SystemExit(1)
            if "contigs" not in df.columns and "amrfinder" not in df.columns:
                log.critical(
                    f"You must supply a path to either assemblies or amrfinder output. Please try again."
                )
                raise SystemExit(1)
            # df = df.rename(columns = {'assembly': 'input', 'amrfinder':'input'})
            df = df.fillna("")
            lines = []
            sids = []
            for row in df.iterrows():
                d = {}

                if row[1]["sample_id"] in sids:
                    log.critical(
                        f"It looks like you have duplicate IDs. Sample IDs must be unique when running complete workflow. abritamr will not run. Please try again."
                    )
                    raise SystemExit(1)
                sids.append(row[1]["sample_id"])
                d["sample_id"] = row[1]["sample_id"]
                d["outdir"] = row[1]["sample_id"]
                for ip in ["amrfinder", "contigs"]:
                    if ip in df.columns.tolist():
                        if not check_path(row[1][f"{ip}"]):
                            log.warning(
                                f"{row[1][f'{ip}']} is not a valid file path. abritamr will not run on this input."
                            )
                        else:
                            d[ip] = row[1][f"{ip}"]

                if row[1]["species"] != "":
                    guess = (
                        guess_species(asm=row[1]["contigs"], sid=row[1]["sample_id"])
                        if row[1]["contigs"] != ""
                        else ""
                    )
                else:
                    guess = row[1]["species"]

                d["species"] = guess

                lines.append(d)

            return lines

        else:
            log.critical(
                f"{pth} does not exist. Please provide a valid path to the multi input file."
            )
            raise SystemExit(1)
    else:
        log.critical(
            f"{workdir} is not a valid path. Please provide a valid working directory."
        )
        raise SystemExit(1)
        # # print(line)


def wrangle_outputs(amr: list, cols: list = []) -> pd.DataFrame:
    """Convert result records into a dataframe with the requested columns."""
    if amr != []:
        df = pd.DataFrame(amr)
        cols = cols if cols != [] else df.columns.tolist()
        df = df[cols]
    else:
        df = pd.DataFrame(columns=cols, data=amr)

    return df


def run(args) -> dict:
    """Run scan, typing, reporting, matrix, and inference workflows."""
    simple = True if args.viewtype == "compact" else False
    # dbv = "unknown"
    if not args.contigs and not args.amrfinderplus and not args.multi:
        log.critical(
            "You must supply an input file (input file with multiple samples OR as single assembly or single amrfinder plus output). Exiting."
        )
        raise SystemExit(1)

    if not args.multi:
        if not args.sample_id:
            log.warning(
                f"You have not supplied a sample id - please note that path to input will be used as sample id"
            )
        if args.outdir == "":
            log.warning(
                f"You have not supplied an output directory. Output files will be generated in your working directory: {args.workdir}"
            )
        # species = (
        #     guess_species(asm=args.contigs[0], sid=args.sample_id)
        #     if args.contigs
        #     else ""
        # )
        inputs = []
        if args.contigs:
            for contig in args.contigs:
                d = generate_inputs(pth=contig, sid=args.sample_id, key="contigs")
                inputs.append(d)
        if args.amrfinderplus:
            for afp in args.amrfinderplus:
                d = generate_inputs(pth=afp, sid=args.sample_id, key="amrfinderplus")
                inputs.append(d)

        for i in inputs:
            i["species"] = args.species
            i["outdir"] = args.outdir
    else:
        inputs = check_multi(pth=args.multi, workdir=pathlib.Path(args.workdir))

    for i in inputs:
        i["threads"] = args.threads
        i["min_identity"] = args.min_identity
        i["min_coverage"] = args.min_coverage
        i["reference_catalog"] = args.reference_catalog

    scanned = run_scan(inputs=inputs)
    typed = do_typing(amr=scanned, reference_catalog=args.reference_catalog)
    linelist = generate_linelist(
        amr=typed,
        _format=args.format,
        viewtype=args.viewtype,
        genesonly=args.genesonly,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
    )
    matrix = make_matrix(
        amr=typed,
        facet=args.facet,
        reference_catalog=args.reference_catalog,
        min_coverage=args.min_coverage,
        min_identity=args.min_identity,
    )
    infer_res = do_gdst(
        amr=typed,
        dflt_result=args.dflt_result,
        reference_folder=args.reference_folder,
        reporttype=args.reporttype,
    )
    outputs = {
        "scan": scanned,
        "typed": typed,
        "linelist": linelist,
        "matrix": matrix,
        "infer": infer_res,
    }

    if args.multi or (len(args.contigs) == 1 and args.sample_id != ""):
        for o in outputs:
            if not outputs[o].empty:
                for i in inputs:
                    idcol = outputs[o].columns.tolist()[0]
                    ts = outputs[o][outputs[o][f"{idcol}"] == i["sample_id"]]
                    save_output(
                        workdir=f"{args.workdir}",
                        sample_id=i["outdir"],
                        result=ts,
                        outname=f"abritamr_{o}",
                        _format=args.format,
                    )
    else:
        for o in outputs:
            if not outputs[o].empty:
                save_output(
                    workdir=f"{args.workdir}",
                    sample_id=args.outdir,
                    result=outputs[o],
                    outname=f"abritamr_{o}",
                    _format=args.format,
                )


# linelists = []
# matrices = []
# infers = []
# r i in inputs:
#     amr = []
#     t
#     dbv = "unknown"
#     # # print(i)
#     res = []
#
#     if "assembly" in i and i["assembly"] != "":
#         if check_path(pth=f"{i['assembly'][0]}") and check_assembly(
#             f"{i['assembly'][0]}"
#         ):
#             dbv = check_amrfinder()
#             log.info(
#                 f"Running amrfinder plus on supplied assembly {i['assembly'][0]} with {args.threads} threads."
#             )
#             res = run_amrf(
#                 min_identity=args.min_identity,
#                 min_coverage=args.min_coverage,
#                 asm=i["assembly"][0],
#                 threads=args.threads,
#                 organism=i["species"],
#             )
#
#     elif "amrfinder" in i and i["amrfinder"] != "":
#         res = amrf2dict(amrfinder=i["amrfinder"])
#     for r in res:
#         r["amrfinderplus_db_version"] = dbv
#     # res = generate_output(species = i['species'], amr = res )
#     amr = apply_classes(
#         amr=res,
#         species=i["species"],
#         sid=i["sample_id"],
#         catalog=args.reference_catalog,
#     )
#     if amr != []:
#         scanned_cols = abritamr_scan_columns()
#         scanned = wrangle_outputs(amr=amr, cols=scanned_cols)
#         save_output(
#             workdir=f"{args.workdir}",
#             sample_id=i["sample_id"],
#             result=scanned,
#             outname="abritamr_scan",
#             _format=args.format,
#             no_keep=args.no_keep,
#         )
#         amr = add_abritamr_results(amr=amr, catalog=args.reference_catalog)
#         typed_cols = abritamr_status_columns()
#         typed = wrangle_outputs(amr=amr, cols=typed_cols, sid=i["sample_id"])
#         save_output(
#             workdir=f"{args.workdir}",
#             sample_id=i["sample_id"],
#             result=typed[typed_cols],
#             outname="abritamr_typed",
#             _format=args.format,
#             no_keep=args.no_keep,
#         )
#         linelist = summary(
#             results=typed,
#             _format=args.format,
#             sid=i["sample_id"],
#             simple=simple,
#             genesonly=args.genesonly,
#             minidentity=args.min_identity,
#             mincoverage=args.min_coverage,
#         )
#         linelists.append(linelist)
#         save_output(
#             workdir=f"{args.workdir}",
#             sample_id=i["sample_id"],
#             result=linelist,
#             outname="abritamr_linelist",
#             _format=args.format,
#             no_keep=False,
#         )
#         infer = gdst(
#             results=scanned,
#             species=species,
#             reference_folder=args.reference_folder,
#             dflt_result=args.dflt_result,
#         )
#         if infer != []:
#             if args.reporttype == "long":
#                 gdstlinelist = gdst_results_to_df_long(infer)
#             elif args.reporttype == "wide":
#                 gdstlinelist = gdst_results_to_df_wide(infer)
#             infers.extend(gdstlinelist)
#             # # print(linelist.columns.tolist())
#             save_output(
#                 workdir=f"{args.workdir}",
#                 sample_id=i["sample_id"],
#                 result=gdstlinelist,
#                 outname="abritamr_gdst",
#                 _format=args.format,
#                 no_keep=False,
#             )
#         matrix = matrix_summary(
#             results=typed,
#             facet=args.facet,
#             sid=i["sample_id"],
#             minidentity=args.min_identity,
#             mincoverage=args.min_coverage,
#             refgenes=args.reference_catalog,
#         )
#         # # print(res)
#         matrices.append(matrix)
#         save_output(
#             workdir=f"{args.workdir}",
#             sample_id=i["sample_id"],
#             result=matrix,
#             outname="abritamr_matrix",
#             _format=args.format,
#             no_keep=False,
#         )
#     else:
#         log.warning(f"{i['sample_id']} did not return any AMR mechanisms")
#         # # print(amr)
# # # print(pd.concat(linelists, axis = 0, ignore_index=True))
# if linelists != []:
#     log.info(f"Saving a single file for output of linelist.")
#     save_output(
#         workdir=f"{args.workdir}",
#         sample_id="",
#         result=pd.concat(linelists).reset_index(drop=True),
#         outname=f"{args.prefix}_linelist",
#         _format=args.format,
#         no_keep=False,
#     )
# if matrices != []:
#     log.info(f"Saving a single for output of matrix.")
#     save_output(
#         workdir=f"{args.workdir}",
#         sample_id="",
#         result=pd.concat(matrices).reset_index(drop=True),
#         outname=f"{args.prefix}_matrix",
#         _format=args.format,
#         no_keep=False,
#     )
# return amr
