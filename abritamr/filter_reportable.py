"""Apply reportability and AMR typing criteria to AMRFinder results."""

import pandas as pd
from abritamr.logger import log
from abritamr.cel_functions import create_cel_context, evaluate_rule


def construct_filter(result: dict, refgenes: pd.DataFrame):
    """Annotate one AMRFinder result with its catalog reporting criteria."""
    result["abritamr_priority_status"] = "not-reportable"
    result["criteria_id"] = ""
    result["criteria_version"] = ""
    result["abritamr_AMR_type"] = ""
    refgenes = refgenes.fillna("")
    row = refgenes[
        (
            refgenes["abritamr_accession_key"].str.contains(
                result["Closest reference accession"], na=False
            )
        )
        & (result["Closest reference accession"] != "-")
    ]
    if not row.empty:
        rw = row.to_dict(orient="records")[0]
        sts = rw["priority_status"]
        tpe = rw["amrtype"]
        cid = rw["criteria_id"]
        civ = rw["criteria_version"]
        status_criteria = rw["additional_status_criteria"]
        tpe_criteria = rw["additional_type_criteria"]
        ctx = create_cel_context(data=result, name="row")
        if status_criteria != "":
            rpt = evaluate_rule(status_criteria, ctx)
            if not rpt:
                sts = f"not-{sts}"
        if tpe_criteria != "":
            trpt = evaluate_rule(tpe_criteria, ctx)
            if not trpt:
                tpe = ""
        result["abritamr_priority_status"] = sts
        result["criteria_id"] = cid
        result["criteria_version"] = civ
        result["abritamr_AMR_type"] = tpe

    else:
        log.warning(
            f"{result['Element symbol']} with accession {result['Closest reference accession']} could not be found."
        )

    return result
