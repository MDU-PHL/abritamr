"""Shared validation, catalog, species, and output utilities."""

from datetime import datetime
import pathlib
import subprocess
from abritamr.version import db
import pandas as pd
import sys
from abritamr.logger import log
from abritamr.run_sourmash import run_sourmash_search
from abritamr.cel_functions import evaluate_rule, create_cel_context


def _get_date():
    """Return today's date in ISO format."""
    return datetime.today().strftime("%Y-%m-%d")


def check_path(pth: str) -> bool:
    """Return whether a filesystem path exists."""
    log.info(f"Checking path: {pth}")
    if not pathlib.Path(pth).exists():
        return False
    return True


# this needs to be user defined!!
def get_refgenes(pth) -> pd.DataFrame:
    """Read the reference gene catalog at the given path."""
    if check_path(pth=pth):
        rf = pathlib.Path(pth)
        return pd.read_csv(rf)
    else:
        log.critical(f"Something is very wrong - the reference catalog is not present.")
        raise SystemExit(1)


def check_any2fasta() -> bool:
    """Check that any2fasta is installed, exiting if it is unavailable."""
    log.info(f"Checking that any2fasta is installed.")
    proc = subprocess.run(
        "any2fasta -v", shell=True, capture_output=True, encoding="utf-8"
    )
    if proc.returncode == 0:
        log.info("any2fasta is installed.")
        return True
    else:
        log.critical(
            f"any2fasta is a dependency of abritamr. Please check your installation."
        )
        raise SystemExit(1)


def abritamr_matrix_columns(refgenes: str, group: str = "abritamr_subclass") -> list:
    """Return sorted matrix columns from the selected catalog field."""
    refs = get_refgenes(pth=refgenes)
    refs["gene"] = refs["allele"].fillna(refs["gene_family"])
    cols_final = sorted(refs[group].unique().tolist())

    return cols_final


def abritamr_amrtype_columns() -> list:
    """Return the ordered columns used by AMR linelist reports."""
    return [
        "Sample_id",
        "Priority AMR mechansims",
        "AMR type",
        "Species provided",
        "Other acquired AMR mechansims",
        "Priority AMR mechanisms (low coverage/identity)",
        "Other core AMR mechansims",
        "Other AMR mechanisms",
        "Other",
    ]


def abritamr_scan_columns() -> list:
    """Return the ordered columns used by scan output."""
    return [
        "sample_id",
        "Element symbol",
        "abritamr_class",
        "abritamr_subclass",
        "abritamr_mechanism",
        # 'abritamr_priority_status',
        # 'abritamr_AMR_type',
        # 'criteria_id',
        # 'criteria_version',
        "species",
        "Element name",
        "Scope",
        "Type",
        "Subtype",
        "Method",
        "pmid",
        "abritamr_class_version",
        "amrfinderplus_db_version",
        "Target length",
        "Reference sequence length",
        "% Coverage of reference",
        "% Identity to reference",
        "Alignment length",
        "Closest reference accession",
        "Closest reference name",
        "HMM accession",
        "HMM description",
        "abritamr_accession_key",
        "amrrules_mutation",
    ]


def abritamr_status_columns() -> list:
    """Return the ordered columns used by AMR status output."""
    return [
        "sample_id",
        "Element symbol",
        "abritamr_class",
        "abritamr_subclass",
        "abritamr_mechanism",
        "abritamr_priority_status",
        "abritamr_AMR_type",
        "species",
        "criteria_id",
        "criteria_version",
        "Element name",
        "Scope",
        "Type",
        "Subtype",
        "Method",
        "pmid",
        "abritamr_class_version",
        "amrfinderplus_db_version",
        "Target length",
        "Reference sequence length",
        "% Coverage of reference",
        "% Identity to reference",
        "Alignment length",
        "Closest reference accession",
        "Closest reference name",
        "HMM accession",
        "HMM description",
        "abritamr_accession_key",
        "amrrules_mutation",
    ]


def check_sourmash() -> bool:
    """Check that sourmash is installed, exiting if it is unavailable."""
    log.info(f"Checking that sourmash is installed.")
    proc = subprocess.run(
        "sourmash -v", shell=True, capture_output=True, encoding="utf-8"
    )
    if proc.returncode == 0:
        log.info("sourmash is installed.")
        return True
    else:
        log.critical(
            f"sourmash is a dependency of abritamr. Please check your installation."
        )
        raise SystemExit(1)


def check_amrfinder() -> bool:
    """Check AMRFinderPlus installation and return its database version."""
    log.info(
        f"Checking that amrfinder plus is installed and database versions are compatible."
    )
    vrsn = subprocess.run(
        "amrfinder -V", shell=True, capture_output=True, encoding="utf-8"
    )
    if vrsn.returncode != 0:
        log.critical(
            f"Something is wrong with your environment. amrfinderplus is required for abritamr to run. Please follow installation instructions and try again."
        )
    else:
        dbv = [i for i in vrsn.stdout.split("\n") if "Database version" in i]
        # # print(dbv)
        dbv = dbv[0].split(":")[-1].strip() if dbv != [] else ""
        if dbv == "":
            log.critical(
                f"It looks like the database version cannont be determined. There may be something wrong with your installation. Please check documentation and try again."
            )
            raise SystemExit(1)
        if db not in dbv:
            log.warning(
                f"Your amrfinder database is at {dbv}. This is different to what is expected ({db}). Please not unexpected behaviour may occur."
            )
            return dbv
        else:
            log.info(
                f"Your amrfinder installation is compatible with current abritamr."
            )
            return dbv


def check_assembly(pth) -> bool:
    """Validate an assembly file with any2fasta."""
    log.info(f"Checking assembly is in a correct format")
    if check_any2fasta():
        # try:
        # # print(f"any2fasta {pth}")
        proc = subprocess.run(
            f"any2fasta {pth}", shell=True, capture_output=True, encoding="utf-8"
        )
        if proc.returncode == 0:
            log.info(
                f"{pth} is a valid assembly file. Will no proceed with running amrfinder."
            )
            return True
        else:
            log.critical(
                f"Your assembly file appears to not be in the correct format. Please try again."
            )
            raise SystemExit(1)
    # except Exception as e:
    #     log.critical(f"Something has gone wrong with any2fasta checking your assembly. The following error was reported: {e}")
    #     raise SystemExit(1)


def guess_species(asm: str, sid: str = "abritamr") -> str:
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
    sp = run_sourmash_search(
        query_filename=asm,
        SBT_path=pathlib.Path(__file__).parent / "species_db",
        sid=sid,
    )
    return sp


def clean_gtdb_species(species: str) -> str:
    """Remove the gtdb prefix and suffices from species"""

    species = species.replace("s__", "").split("_")[0]

    return species


def wrangle_species(
    organism: str, asm: str = "", species_rules: str = "", check_species: bool = True
) -> tuple:
    """Map a species name to the AMRFinderPlus organism option."""
    try:
        species_rules = pd.read_csv(f"{species_rules}").to_dict(orient="records")
    except Exception as e:
        log.critical(f"Something has gone very wrong : {e}.")
        raise SystemExit

    if not organism and check_species and asm != "":
        log.info(
            "No species supplied - will use sourmash to try to guess the best match for AMR classification."
        )

        organism = guess_species(asm)
    organism = clean_gtdb_species(species=organism)
    ctx = create_cel_context(data={"species": [organism]}, name="row")
    for rule in species_rules:
        # log.info(rule)
        species_criteria = rule["criteria"]
        # log.info(f"Checking rule: {species_criteria}")
        if species_criteria != "":
            rpt = evaluate_rule(species_criteria, ctx)
            if rpt:
                log.info(f"{rule['criteria_id']} has returned true")

                return f"-O {rule['afpspecies']}"

    return ""


def output_results(
    df: pd.DataFrame,
    output: str = "",
    _format: str = "csv",
    workdir: str = f"{pathlib.Path.cwd()}",
) -> bool:
    """Write a dataframe to stdout or a file in CSV or tab-delimited format."""
    dlm = "," if _format == "csv" else "\t"
    suf = "csv" if _format == "csv" else "txt"

    if not output:
        df.to_csv(sys.stdout, sep=dlm, index=False)
    else:
        p = pathlib.Path(workdir / f"{output}.{suf}")
        df.to_csv(p, sep=dlm, index=False)
