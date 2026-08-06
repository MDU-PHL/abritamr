import pandas as pd
import numpy as np
import logging
import json
import pathlib
import sys
from abritamr.logger import log
from abritamr.utils import get_refgenes

from abritamr.filter_reportable import construct_filter


# logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', filename=sys.stderr, datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
# handler = logging.StreamHandler(sys.stderr)

# log = logging.getLogger(__name__)
# log.addHandler(handler)
# log.setLevel(logging.DEBUG)

def add_abritamr_results(amr:dict,  cfgpath:str="") -> pd.DataFrame:
    log.info(f"Adding abritamr classes and determining gene status.")
    refgenes = get_refgenes()
    for row in amr:
        row['gene'] = row['Element symbol']
        row = construct_filter( result = row , cfgpath = cfgpath)
        
        
        
        # print(row)
    return amr



    