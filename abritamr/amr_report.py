import pandas as pd
import numpy as np

from abritamr.inferrence import infer
from abritamr.logger import log
from abritamr.utils import get_refgenes
import json
import pathlib
import logging

# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
# log = logging.getLogger(__name__)
# log.setLevel(logging.DEBUG)


def save_report(report: pd.DataFrame, outname:str= "", _format:str="csv") -> bool:
    
    dlm= ","
    if _format == "tab":
        dlm = "\t"
        _format = "txt"
    
    report.to_csv(f'{outname}.{_format}', sep = dlm, index = False)

def wrangle_cols(repdf:pd.DataFrame, repmechs:dict, cols:list, simple:bool=True) -> tuple:
    cols_final = sorted(repdf['abritamr_subclass'].unique().tolist())
    # print(cols_final)
    # cols_final = []
    if not simple:
        refs = get_refgenes()
        cols_final = sorted(refs['abritamr_subclass'].unique().tolist())
    for dc in cols_final:
        
        tmp = repdf[repdf['abritamr_subclass'] == dc]
        if not tmp.empty:
            repmechs[dc] = ';'.join(sorted(tmp['Element symbol'].unique().tolist()))
        else:
            repmechs[dc] = ""
        cols.append(dc)
    
    return repmechs,cols

def summary(
    results:pd.DataFrame,
    _format:str="csv", 
    simple:bool=False, 
    sid:str="", 
    genesonly:bool = False, 
    minidentity:float = 90, 
    mincoverage:float = 90,
    outname:str = "abritamr_report"
    ) -> bool:
    
    mincoverage = float(mincoverage)
    minidentity = float(minidentity)
    results = results.fillna("")
    log.info(f"Generating report for {sid}")
    results['% Coverage of reference'] =  pd.to_numeric(results['% Coverage of reference'] , errors='coerce')
    results['% Identity to reference'] = pd.to_numeric(results['% Identity to reference'] , errors='coerce')
    
    species = results['species'].unique().tolist()
    if len(species) >1:
        log.critical(f"It seems that there is something wrong with your input files. There is more than one species value ({', '.join(species)})")
        raise SystemExit(1)
    repmechs = {'Sample_id':sid, 'Species provided': species[0] if species != [] else ""}
    reportable = results[
        (results['abritamr_AMR_reporting_status'] == 'reportable') & 
        (results['% Identity to reference']>=minidentity) &
        (results['% Coverage of reference']>=mincoverage)]
    if genesonly:
        reportable = reportable[~reportable['Subtype'].str.contains('POINT')]
    reportable_amr_low=results[
        (results['abritamr_AMR_reporting_status'] == 'reportable') & 
        (
            (results['% Identity to reference']<float(minidentity)) |
            (results['% Coverage of reference']<float(mincoverage))
        ) &
        (results['Type'] == 'AMR')]
    nonreportable_amr = results[
        (~results['Element symbol'].isin(reportable['Element symbol'].unique().tolist())) & 
        (results['% Identity to reference']>=float(minidentity)) &
        (results['% Coverage of reference']>=float(mincoverage)) &
        (results['Type'] == 'AMR')]
    nonreportable_other=results[(results['Type'] != 'AMR')]
    amrtype = results["abritamr_AMR_type"].unique().tolist()
    repmechs['Reportable AMR mechansims']=  ','.join(sorted(reportable['Element symbol'].unique().tolist()))
    if len(amrtype) == 1 and amrtype[0] == "":
        amrtype= ""
    else:
        amrtype = ";".join([a for a in amrtype if a])
    cols = ['Sample_id','Reportable AMR mechansims', 'AMR type','Species provided', 'Non-reportable AMR mechanisms','Reportable AMR mechanisms (low coverage/identity)','Non-reportable other']
    # print(cols)
    # for df in [reportable, nonreportable_amr, nonreportable_other]:
    df = results[
        (results['% Identity to reference']>=minidentity) &
        (results['% Coverage of reference']>=mincoverage)]
    repmechs,cols = wrangle_cols(df, repmechs, cols, simple = simple)
        
    repmechs['Non-reportable AMR mechanisms']=','.join(sorted(nonreportable_amr['Element symbol'].unique().tolist()))
    repmechs['Reportable AMR mechanisms (low coverage/identity)'] = ','.join(sorted(reportable_amr_low['Element symbol'].unique().tolist()))
    repmechs['Non-reportable other']=','.join(sorted(nonreportable_other['Element symbol'].unique().tolist()))
    repmechs["Species provided"] = species
    repmechs["AMR type"] = amrtype
    
    report = pd.DataFrame(repmechs, index = [0])
    report=report[cols]

    return report


def infer_phenotype(amr:dict, species: str= "", genus:str="", default:str="Susceptible", sid:str="abritamr", output:str="abritamr_inferred_ast.csv") -> bool:
    
    colorder = ["seq_id"]
    inferred = infer(amr = amr, species=species, default = default)
    indf = pd.DataFrame(inferred, index = [0])
    cols = sorted(indf.columns.tolist())
    indf["ID"] = sid
    colorder.extend(cols)
    
    indf[colorder].to_csv(f"{output}", index = False)

    return True
