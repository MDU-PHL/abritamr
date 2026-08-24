
import pathlib
import json
import logging
import sys

from abritamr.logger import log
from abritamr.parse_gdstrules import generate_rules
from abritamr.utils import check_path
from abritamr.commands.update_database import create_db_folder,update_rules


def rules(args) -> bool:
    if create_db_folder(args.output_dir):
        log.info(f"Database folder created at {args.output_dir}")
        if check_path(args.catalog):
            log.info(f"Will generate abritamr compatible rules from {args.catalog}")
            try:
                update_rules(args) 
            except Exception as e:
                log.critical(f"Looks like something has gone wrong with generating the rules. The following error was reported : {e}. Please try again.")
                raise SystemExit(1)
        else:
            log.critical(f"Cannot find the reference gene catalog at {args.catalog}. Please check your input and try again.")
            raise SystemExit(1)
        