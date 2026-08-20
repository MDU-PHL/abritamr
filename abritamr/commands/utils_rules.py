import click
import pathlib
import json
import logging
import pandas as pd
import sys

from abritamr.logger import log
from abritamr.parse_amrrules import generate_rules
from abritamr.utils import check_path


def utils_rules(args) -> bool:

    with open(f"{pathlib.Path(__file__).parent.parent / 'configs' / 'amr_rules_config.json'}", "r") as j:
        cfg = json.load(j)

    
    eg = cdg["grade"][args.evidence_grade]
    rules = generate_rules(species=args.species, evidence_grade=ev, cfg=cfg)


