import click
import pathlib
import json
import logging
import pandas as pd
import sys
from io import StringIO

from abritamr.utils import check_assembly, check_amrfinder, check_any2fasta, wrangle_species, output_results
# from abritamr.run_finder import run_amrf
# from abritamr.parse_finder import amrf2dict
# from abritamr.parse_reportable import add_abritamr_results
# from abritamr.parse_amrtype import get_amr_type
from abritamr.amr_report import summary


logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
handler = logging.StreamHandler(sys.stderr)

log = logging.getLogger(__name__)
log.addHandler(handler)



def linelist(
    args
        ) -> dict:

    # try:
    # amr = pd.DataFrame(StringIO(args.amr.read()), sep = ",")
    print(args.amr.read())
    # except:
    #     print('no input')
    