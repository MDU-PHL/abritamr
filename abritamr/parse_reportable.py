import pandas as pd

from abritamr.utils import get_refgenes

from abritamr.filter_reportable import construct_filter
from abritamr.logger import log


def get_class(key: str, val: str, refgenes: pd.DataFrame, _type: str) -> str:

    try:
        return refgenes[refgenes[key] == val][f"{_type}"].values[0]
    except:
        log.critcal("The entry is not present in the reference catalog")
        return "unknown"


def get_classes(key: str, val: str, refgenes: pd.DataFrame) -> tuple:

    _class = get_class(key=key, val=val, refgenes=refgenes, _type="abritamr_class")
    _subclass = get_class(
        key=key, val=val, refgenes=refgenes, _type="abritamr_subclass"
    )

    return _class, _subclass


def find_classes(refgenes: pd.DataFrame, accession: str) -> str:

    _class = _subclass = "unknown"

    for key in [
        "refseq_protein_accession",
        "refseq_nucleotide_accession",
        "genbank_protein_accession",
        "genbank_nucleotide_accession",
    ]:
        if accession in refgenes[key].unique().tolist():
            _class, _subclass = get_classes(key=key, val=accession, refgenes=refgenes)
            break

    return _class, _subclass


def add_abritamr_results(amr: dict, sid: str = "", catalog: str = "") -> pd.DataFrame:
    log.info(f"Adding abritamr classes and determining gene status.")
    refgenes = get_refgenes(pth=catalog)
    for row in amr:
        if row["sample_id"] == "":
            log.critcal(f"You must have a sample_id value.")
            raise SystemExit(1)
        row["gene"] = row["abritamr_mechanism"]
        row = construct_filter(result=row, refgenes=refgenes)
    return amr
