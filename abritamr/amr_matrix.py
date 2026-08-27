import pandas as pd
import numpy as np

# from abritamr.inferrence import infer
from abritamr.logger import log
from abritamr.utils import get_refgenes
import json
import pathlib
import logging



def wrangle_cols(repdf:pd.DataFrame, repmechs:dict, cols:list, group:str= 'abritamr_subclass',refgenes:str="") -> tuple:
    cols_final = sorted(repdf['abritamr_subclass'].unique().tolist())

    refs = get_refgenes(pth = refgenes)
    refs['gene'] = refs['allele'].fillna(refs['gene_family'])
    cols_final = sorted(refs[group].unique().tolist())
    
    for dc in cols_final:
        g = group if group != "gene" else 'Element symbol'
        tmp = repdf[repdf[g] == dc]
        if not tmp.empty and group != "gene":
            repmechs[dc] = ';'.join(sorted(tmp['Element symbol'].unique().tolist()))
        elif not tmp.empty and group == "gene":
            repmechs[dc] = len(sorted(tmp['Element symbol'].unique().tolist()))
        elif tmp.empty and group == "gene":
            repmechs[dc] = 0
        else:
            repmechs[dc] = '-'
        cols.append(dc)
    
    return repmechs,cols

def summary(
    results:pd.DataFrame,
    facet:str,
    sid:str="", 
    minidentity:float = 90, 
    mincoverage:float = 90,
    refgenes:str = ""
    ) -> bool:
    
    mincoverage = float(mincoverage)
    minidentity = float(minidentity)
    results = results.fillna("")
    log.info(f"Generating matrix for {sid}")
    results['% Coverage of reference'] =  pd.to_numeric(results['% Coverage of reference'] , errors='coerce')
    results['% Identity to reference'] = pd.to_numeric(results['% Identity to reference'] , errors='coerce')
    
    repmechs = {'Sample_id':sid}
    reportable = results[
        (results['% Identity to reference']>=minidentity) &
        (results['% Coverage of reference']>=mincoverage)]
    cols = ['Sample_id']
    
    
    repmechs,cols = wrangle_cols(repdf = reportable, repmechs=repmechs, cols= cols,group = facet,refgenes = refgenes)
    report = pd.DataFrame(repmechs, index = [0])
    report=report[cols]

    return report

