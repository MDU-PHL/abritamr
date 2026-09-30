"""Look up AMR classes and mechanisms in the reference gene catalog."""

import pandas as pd
import numpy as np
from abritamr.logger import log
from abritamr.utils import get_refgenes
import warnings

warnings.simplefilter(action="ignore", category="UserWarning")


def get_class(key: str, val: str, refgenes: pd.DataFrame, _type: str) -> str:
    """Look up a catalog value by key and return the requested field."""
    try:
        if refgenes[refgenes[key] == val].empty:
            return "NA"
        return refgenes[refgenes[key] == val][f"{_type}"].values[0]
    except:
        log.critical("The entry is not present in the reference catalog")
        return "unknown"


def get_amrrules_mutation(val: str, refgenes: pd.DataFrame) -> str:
    """Look up the AMR-rule mutation associated with an accession key."""
    try:
        if refgenes[refgenes["abritamr_accession_key"] == val].empty:
            return "-"
        return refgenes[refgenes["abritamr_accession_key"] == val][
            f"amrrules_mutation"
        ].values[0]
    except:
        return "-"


def get_mechanism(val: str, refgenes: pd.DataFrame, symbol: str) -> str:
    """Find the mechanism matching an accession key and element symbol."""
    try:
        if refgenes[
            refgenes["abritamr_accession_key"].str.contains(val, regex=False)
        ].empty:
            return "-"

        elif not refgenes[
            (refgenes["abritamr_accession_key"].str.contains(val, regex=False))
            & (refgenes["abritamr_mechanism"].str.contains(symbol, regex=False))
        ].empty:
            return refgenes[
                (refgenes["abritamr_accession_key"].str.contains(val, regex=False))
                & (refgenes["abritamr_mechanism"].str.contains(symbol, regex=False))
            ]["abritamr_mechanism"].values[0]

        return refgenes[
            refgenes["abritamr_accession_key"].str.contains(val, regex=False)
        ][f"abritamr_mechanism"].values[0]
    except Exception as e:
        log.critical(
            f"The following error has occured {e} and mechanism {val}({symbol}) cannot be identified not present in the reference catalog"
        )
        return "unknown"


def get_accession_key(key: str, val: str, symbol: str, refgenes: pd.DataFrame) -> str:
    """Resolve a hit to its catalog accession key."""
    try:
        if refgenes[
            refgenes["abritamr_accession_key"].str.contains(val, regex=False)
        ].empty:
            return "-"
        elif not refgenes[
            (refgenes["abritamr_accession_key"].str.contains(val, regex=False))
            & (refgenes["abritamr_mechanism"].str.contains(symbol, regex=False))
        ].empty:
            return refgenes[
                (refgenes["abritamr_accession_key"].str.contains(val, regex=False))
                & (refgenes["abritamr_mechanism"].str.contains(symbol, regex=False))
            ]["abritamr_accession_key"].values[0]

        return refgenes[refgenes[key] == val][f"abritamr_accession_key"].values[0]
    except Exception as e:
        log.critical(
            f"The following error has occured {e} and accession for {val} ({symbol} cannot not identified in the reference catalog"
        )
        return "unknown"


def get_classes(key: str, val: str, refgenes: pd.DataFrame, symbol: str) -> tuple:
    """Return class, subclass, provenance, accession, mutation, and mechanism."""
    abacc = get_accession_key(key=key, val=val, refgenes=refgenes, symbol=symbol)
    # log.info(f"Found {abacc}")
    _class = get_class(key=key, val=val, refgenes=refgenes, _type="abritamr_class")
    _subclass = get_class(
        key=key, val=val, refgenes=refgenes, _type="abritamr_subclass"
    )
    pmid = get_class(key=key, val=val, refgenes=refgenes, _type="pubmed_reference")
    db_version = get_class(key=key, val=val, refgenes=refgenes, _type="db_version")
    key = get_class(key=key, val=val, refgenes=refgenes, _type="abritamr_accession_key")
    amrrules_mut = get_amrrules_mutation(val=abacc, refgenes=refgenes)
    mech = get_mechanism(val=abacc, refgenes=refgenes, symbol=symbol)

    return _class, _subclass, pmid, db_version, key, amrrules_mut, mech


def find_classes(refgenes: pd.DataFrame, accession: str, symbol: str) -> str:
    """Look up catalog annotations by any supported reference accession."""
    _class = _subclass = pmid = db_version = akey = amrrules_mut = mech = "-"

    for key in [
        "refseq_protein_accession",
        "refseq_nucleotide_accession",
        "genbank_protein_accession",
        "genbank_nucleotide_accession",
    ]:
        if accession in refgenes[key].unique().tolist():
            _class, _subclass, pmid, db_version, akey, amrrules_mut, mech = get_classes(
                key=key, val=accession, refgenes=refgenes, symbol=symbol
            )
            break

    return _class, _subclass, pmid, db_version, akey, amrrules_mut, mech


def apply_classes(amr: dict, species: str, sid: str, catalog: str) -> dict:
    """Annotate AMRFinder hits with catalog classes and sample metadata."""
    refgenes = get_refgenes(pth=catalog)
    for row in amr:
        _class, _subclass, pmid, db_version, key, amrrules_mut, mech = find_classes(
            refgenes=refgenes,
            accession=row["Closest reference accession"],
            symbol=row["Element symbol"],
        )
        row["abritamr_class"] = _class
        row["abritamr_subclass"] = _subclass
        row["species"] = species
        row["sample_id"] = sid
        row["pmid"] = pmid
        row["abritamr_class_version"] = db_version
        row["abritamr_accession_key"] = key
        row["amrrules_mutation"] = amrrules_mut
        row["abritamr_mechanism"] = mech
    return amr
