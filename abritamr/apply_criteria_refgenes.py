import pandas as pd
import pathlib
from dataclasses import asdict
from logger import log
from utils import get_refgenes
from criteria import get_abritamr_configs
from cel import evaluate
import sys

# reporting_criteria = pd.read_csv(f"{sys.argv[1]}")
rcdict = get_abritamr_configs(cfgpath = sys.argv[1])
refgenes = get_refgenes()
refgenes['gene'] = refgenes['allele'].fillna(refgenes['gene_family'])

dflt_rs = 'not-reportable'

refgenes = refgenes.to_dict(orient = "records")
rows = []
for row in refgenes:
    row['reportable_status'] = dflt_rs
    for cr in rcdict:
        crt = asdict(cr)
        rs = evaluate(crt['criteria'], row)
        if rs:
            print(row)
            print(crt['status'])
            row['reportable_status'] = crt['status']
            row['amrtype'] = crt['amrtype']
            row['criteria_id'] = crt['criteria_id']
            row['criteria_version'] = crt['criteria_version']
            row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
            row['additional_status_criteria'] = crt['additional_status_criteria'] if crt['additional_status_criteria'] else ""
            row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
            row['additional_type_criteria'] = crt['additional_type_criteria'] if crt['additional_type_criteria'] else ""
            break
            
        # else:
        #     row['reportable_status'] = dflt_rs
    rows.append(row)
refgenes2 = pd.DataFrame(rows)

refgenes2.to_csv(f"refgenes_v2.csv", index = False)
