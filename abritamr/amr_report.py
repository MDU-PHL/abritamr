import pandas as pd
import numpy as np

from inferrence import infer

from parse_amrtype import get_amr_type
from parse_reportable import add_abritamr_results

import json
import pathlib



def save_report(report: pd.DataFrame, outname:str= "abritamr_reportable", long:bool = True) -> bool:
    if long:
        report.T.to_csv(f'{outname}.csv', header = False)
    else:
        report.to_csv(f'{outname}.csv', index = False)

def wrangle_cols(repdf:pd.DataFrame, repmechs:dict, cols:list) -> tuple:

    for dc in repdf['abritamr_subclass'].unique().tolist():
            tmp = repdf[repdf['abritamr_subclass'] == dc]
            repmechs[dc] = ','.join(sorted(tmp['Element symbol'].unique().tolist()))
            cols.append(dc)
    
    return repmechs,cols

def summary(results:pd.DataFrame,species:str="", genus:str="", simple:bool=False, sid:str="abritamr", genesonly:bool = False, minidentity:float = 90, mincoverage:float = 90) -> pd.DataFrame:

    repmechs = {'Sample_id':sid}
    reportable = results[
        (results['abritamr_AMR_reporting'] == 'reportable') & 
        (results['% Identity to reference']>=minidentity) &
        (results['% Coverage of reference']>=mincoverage)]
    if genesonly:
        reportable = reportable[~reportable['Element subtype'].str.contains('POINT')]
    reportable_amr_low=results[
        (results['abritamr_AMR_reporting'] == 'reportable') & 
        (
            (results['% Identity to reference']<minidentity) |
            (results['% Coverage of reference']<mincoverage)
        ) &
        (results['Type'] == 'AMR')]
    nonreportable_amr = results[
        (~results['Element symbol'].isin(reportable['Element symbol'].unique().tolist())) & 
        (results['% Identity to reference']>=minidentity) &
        (results['% Coverage of reference']>=mincoverage) &
        (results['Type'] == 'AMR')]
    nonreportable_other=results[(results['Type'] != 'AMR')]
    amrtype = results["abritamr AMR type"].unique().tolist()
    repmechs['Reportable AMR mechansims']=  ','.join(sorted(reportable['Element symbol'].unique().tolist()))
    if len(amrtype) == 1 and amrtype[0] == "No known type":
        amrtype= "No known type"
    else:
        amrtype = ", ".join([a for a in amrtype if a != "No known type"])
    cols = ['Sample_id','Reportable AMR mechansims', 'AMR type', 'Non-reportable AMR mechanisms','Reportable AMR mechanisms (low coverage/identity)','Non-reportable other','Species provided']
    if not simple:
        for df in [reportable, nonreportable_amr, nonreportable_other]:
            repmechs,cols = wrangle_cols(df, repmechs, cols)
        
    repmechs['Non-reportable AMR mechanisms']=','.join(sorted(nonreportable_amr['Element symbol'].unique().tolist()))
    repmechs['Reportable AMR mechanisms (low coverage/identity)'] = ','.join(sorted(reportable_amr_low['Element symbol'].unique().tolist()))
    repmechs['Non-reportable other']=','.join(sorted(nonreportable_other['Element symbol'].unique().tolist()))
    repmechs["Species provided"] = species
    repmechs["AMR type"] = amrtype
    
    report = pd.DataFrame(repmechs, index = [0])
    report=report[cols]
    # print(report.T)
    return report
    # report.T.to_csv('abritamr_report.csv', header = False)

def reportable(amrfinder:str='amrfinder.out', species:str="", genus:str="", sid:str="abritamr") -> pd.DataFrame:
    try:
        amrdf = pd.read_csv(amrfinder, sep = "\t")
    except Exception as e:
        print(f"An error has occured opening the amrfinder output : {e}.\n Please try again.")
        raise SystemExit

    amr = amrdf.to_dict(orient = "records")
    
    reportable_amr = add_abritamr_results( amr = amr )
    reportable_amr = get_amr_type( amr = reportable_amr)
    # print(reportable_amr)

    return reportable_amr


def infer_phenotype(amr:dict, species: str= "", genus:str="", default:str="Susceptible", sid:str="abritamr", output:str="abritamr_inferred_ast.csv") -> bool:
    
    colorder = ["seq_id"]
    inferred = infer(amr = amr, species=species, default = default)
    indf = pd.DataFrame(inferred, index = [0])
    cols = sorted(indf.columns.tolist())
    indf["ID"] = sid
    colorder.extend(cols)
    
    indf[colorder].to_csv(f"{output}", index = False)

    return True

rep = reportable(amrfinder = "salmo.out", species = "Salmonella enterica", genus = "Salmonella")

summary(results = pd.DataFrame(rep),species = "Salmonella enterica", genus = "Salmonella")
infer_phenotype(amr = rep, species = "Salmonella enterica")
# for r in rep:
#     print(r)
# report = report(result, simple= True)

# save_report(report= report,long = False)