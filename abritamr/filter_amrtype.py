"""Build and evaluate filters for AMR typing criteria."""

from abritamr.criteria import get_abritamr_configs
from dataclasses import dataclass, asdict
from cel import evaluate


def filter_string(crt: dict, dlm: str = " && ") -> str:
    """Convert a criteria mapping into a CEL expression."""
    fl = []

    for k in crt:
        if f"{crt[k]}" != "None" and k != "exception" and k != "amrtype":
            if k not in ["abritamr_class", "abritamr_subclass", "gene"]:
                f = f"{k} in {crt[k]}"
            else:
                lst = eval(crt[k]) if not isinstance(crt[k], list) else crt[k]
                tmp = []
                for l in lst:
                    tmp.append(f"'{l}' in {k}")
                f = " || ".join(tmp)

            fl.append(f)
    return f"{dlm}".join(fl)


def extract_amrtype(criteria: dict):
    """Return a criteria object's fields except its AMR type label."""
    crt = {k: str(v) for k, v in asdict(criteria).items() if k not in ["amrtype"]}

    return crt


def construct_filter(result: dict, cfgpath: str = ""):
    """Evaluate configured typing rules and return the matching AMR type."""
    rpt = "No known type"
    listofamrtypes = get_abritamr_configs(cfgtype="amr_type", cfgpath=cfgpath)
    for criteria in listofamrtypes:
        # # print(criteria)
        crt = extract_amrtype(criteria=criteria)
        amrtp = {k: str(v) for k, v in asdict(criteria).items() if k == "amrtype"}
        primary_filter = filter_string(crt)

        res = evaluate(primary_filter, result)
        if res:
            rpt = amrtp["amrtype"]
        if "exception" in crt and crt["exception"] != "None":
            sfilterres = set()
            for xcptn in eval(crt["exception"]):
                flt = filter_string(xcptn, dlm=" && ")

                resxptn = evaluate(flt, result)
                if resxptn:
                    sfilterres.add(xcptn["amrtype"])
            if sfilterres:
                rpt = "|".join(list(sfilterres))

    return rpt
