
import pathlib
import json
import logging

import sys

from abritamr.logger import log
from abritamr.parse_gdstrules import generate_rules
from abritamr.utils import check_path
from abritamr.commands.update_database import create_db_folder,update_catalog


def catalog(args) -> bool:
    if create_db_folder(args.output_dir):
        log.info(f"Database folder created at {args.output_dir}")

        try:
            update_catalog(args) 
        except Exception as e:
            log.critical(f"Looks like something has gone wrong with generating the catalog. The following error was reported : {e}. Please try again.")
            raise SystemExit(1)