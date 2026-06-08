import pandas as pd
import numpy as np

import json
import pathlib



def save_report(report: pd.DataFrame, outname:str= "abritamr_reportable", long:bool = True) -> bool:
    if long:
        report.T.to_csv(f'{outname}.csv', header = False)
    else:
        report.to_csv(f'{outname}.csv', index = False)

    

def report(results:pd.DataFrame,species:str="", simple:bool=False, sid:str="abritamr", genesonly:bool = False, minidentity:float = 90, mincoverage:float = 90) -> pd.DataFrame:

    repmechs = {'seq_id':sid}
    reportable = results[
        (results['abritamr AMR reporting'] == 'reportable') & 
        (results['% Identity to reference sequence']>=minidentity) &
        (results['% Coverage of reference sequence']>=mincoverage)]
    if genesonly:
        reportable = reportable[~reportable['Element subtype'].str.contains('POINT')]
    reportable_amr_low=results[
        (results['abritamr AMR reporting'] == 'reportable') & 
        (
            (results['% Identity to reference sequence']<minidentity) |
            (results['% Coverage of reference sequence']<mincoverage)
        ) &
        (results['Element type'] == 'AMR')]
    nonreportable_amr = results[
        (~results['Gene symbol'].isin(reportable['Gene symbol'].unique().tolist())) & 
        (results['% Identity to reference sequence']>=minidentity) &
        (results['% Coverage of reference sequence']>=mincoverage) &
        (results['Element type'] == 'AMR')]
    nonreportable_other=results[(results['Element type'] != 'AMR')]
    cols = ['seq_id']
    if simple:
        repmechs['Reportable AMR mechansims']=  ','.join(sorted(reportable['Gene symbol'].unique().tolist()))
        cols.append('Reportable AMR mechansims')
    else:
        for dc in reportable['abritamr subclass'].unique().tolist():
            tmp = reportable[reportable['abritamr subclass'] == dc]
            repmechs[dc] = ','.join(sorted(tmp['Gene symbol'].unique().tolist()))
            cols.append(dc)
    repmechs['Non-reportable AMR mechanisms']=','.join(sorted(nonreportable_amr['Gene symbol'].unique().tolist()))
    repmechs['Reportable AMR mechanisms (low coverage/identity)'] = ','.join(sorted(reportable_amr_low['Gene symbol'].unique().tolist()))
    repmechs['Non-reportable other']=','.join(sorted(nonreportable_other['Gene symbol'].unique().tolist()))
    repmechs["Species provided"] = species
    cols.extend(['Non-reportable AMR mechanisms','Reportable AMR mechanisms (low coverage/identity)','Non-reportable other','Species provided'])
    
    report = pd.DataFrame(repmechs, index = [0])
    report=report[cols]
    return report
    # report.T.to_csv('abritamr_report.csv', header = False)

result = pd.read_csv('abritamr.csv')

report = report(result, simple= True)

save_report(report= report,long = False)