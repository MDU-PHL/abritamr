import click
import pathlib
import json
import logging
import pandas as pd
import sys

from abritamr.logger import log
from abritamr.catalog import wrangle_catalog


def update_catalog(args) -> bool:

    # try:
        created = wrangle_catalog(
            catalog = args.catalog,
            previous_catalog = args.previous_catalog,
            amrtyping_definitions = args.amrtyping_definitions,
            class_definitions = args.class_definitions, 
            output= args.output)

        log.info(f"The basic reference gene catalogue has been created: {args.output}")
        return created
    # except Exception as e:
    #     log.critical(f"Looks like something has gone wrong with generating the reference catalogue. The following error was reported : {e}. Please try again.")
    # catalog:str, previous_catalog:str, amrtyping_definitions:str, class_defintitions:str output:str, src:str = "abritamr"