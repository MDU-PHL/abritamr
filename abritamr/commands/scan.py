"""Run AMRFinderPlus scans or annotate existing AMRFinder results."""

import pathlib
import pandas as pd

from collections import namedtuple
from abritamr.utils import (
    check_assembly,
    check_amrfinder,
    guess_species,
    check_path,
    abritamr_scan_columns,
)
from abritamr.run_finder import run_amrf
from abritamr.parse_finder import amrf2dict
from abritamr.drugclasses import apply_classes
from abritamr.logger import log


def generate_output(species: str, sample_id: str, amr: list, catalog: str) -> list:
    """Annotate AMRFinder hits with catalog classes and sample metadata."""
    if len(amr) == 1 and set(amr[0].values()) == set("-"):
        log.warning(f"There are no genes detected for {sample_id}.")
    amr = apply_classes(amr=amr, species=species, sid=sample_id, catalog=catalog)

    return amr


def runfindr(
    asm: str,
    sid: str,
    min_coverage: float,
    min_identity: float,
    threads: int,
    species: str,
    reference_catalog: str,
) -> list:
    """Validate an assembly, run AMRFinderPlus, and annotate its results."""
    res = []
    log.info("Assembly(ies) have been supplied.")
    if check_path(pth=f"{asm}") and check_assembly(f"{asm}"):
        species = (
            species
            if species
            else guess_species(asm=asm, sid=sid if sid else "abritamr")
        )
        dbv = check_amrfinder()

        full_path = f"{pathlib.Path(f'{asm}').absolute()}"
        log.info(f"Running amrfinder plus on {sid}")
        res = run_amrf(
            min_identity=min_identity,
            min_coverage=min_coverage,
            asm=asm,
            threads=threads,
            organism=species,
        )
        res = generate_output(
            species=species,
            sample_id=sid,
            amr=res,
            catalog=reference_catalog,
        )
    else:
        log.critical(f"{asm} is not a valid file. Nothing further will be run")
    return res


def prsfindr(afp: str, species: str, sid: str, reference_catalog: str) -> list:
    """Read existing AMRFinderPlus output and annotate its results."""
    species = species if species else ""
    full_path = f"{pathlib.Path(f'{afp}').absolute()}"
    log.info(f"Opening existing amrfinder plus output")

    res = amrf2dict(amrfinder=afp)
    res = generate_output(
        species=species,
        sample_id=sid,
        amr=res,
        catalog=reference_catalog,
    )

    return res


def run_scan(
    inputs: list,
) -> pd.DataFrame:
    """Scan all supplied inputs and return standardized AMR hit records."""
    amr = []
    dbv = "unknown"

    for input in inputs:
        Data = namedtuple("Data", input.keys())
        data = Data(**input)
        if not data.contigs and not data.amrfinderplus:
            log.critical(
                "You must supply an input file (assembly or amrfinder plus output)."
            )

        if data.contigs:
            res = runfindr(
                asm=data.contigs,
                sid=data.sample_id,
                min_identity=data.min_identity,
                min_coverage=data.min_coverage,
                threads=data.threads,
                reference_catalog=data.reference_catalog,
                species=data.species,
            )

            amr.extend(res)
        elif data.amrfinderplus:
            res = prsfindr(
                sid=data.sample_id,
                reference_catalog=data.reference_catalog,
                afp=data.amrfinderplus,
                species=data.species,
            )
            amr.extend(res)
    abritamr_columns = abritamr_scan_columns()
    if amr != []:
        amr = pd.DataFrame(amr)
    else:
        amr = pd.DataFrame(columns=abritamr_columns)
    amr["amrfinderplus_db_version"] = dbv
    amr = amr[abritamr_columns]

    return amr


def generate_inputs(pth: str, sid: str, key: str) -> dict:
    """Create an input record using a supplied or path-derived sample ID."""
    filepath = f"{pathlib.Path(pth).resolve()}"
    s_id = sid if sid else filepath
    d = {key: pth, "sample_id": s_id}
    return d


def scan(args) -> list:
    """Build scan inputs from CLI arguments and return their results."""
    inputs = []
    if not args.contigs and not args.amrfinderplus:
        log.critical(
            f"There is no input supplied. You must supply contigs or amrfinder outputs as input files. Please try again."
        )
        raise SystemExit(1)
    if args.contigs:
        for contig in args.contigs:
            d = generate_inputs(pth=contig, sid=args.sample_id, key="contigs")
            inputs.append(d)
    if args.amrfinderplus:
        for afp in args.amrfinderplus:
            d = generate_inputs(pth=afp, sid=args.sample_id, key="amrfinderplus")
            inputs.append(d)
    for i in inputs:
        i["threads"] = args.threads
        i["species"] = args.species
        i["min_coverage"] = args.min_coverage
        i["min_identity"] = args.min_identity
        i["reference_catalog"] = args.reference_catalog

    amr = run_scan(inputs=inputs)

    return amr
