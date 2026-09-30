"""Read AMRFinderPlus output into record dictionaries."""

import pandas as pd
import pathlib
from abritamr.logger import log


def amrf2dict(amrfinder: str) -> dict:
    """Read a tab-delimited AMRFinderPlus result file into records."""
    if pathlib.Path(amrfinder).exists():
        try:
            df = pd.read_csv(amrfinder, sep="\t")

            amr = df.to_dict(orient="records")
            return amr

        except Exception as e:
            log.critical(
                f"Something has gone wrong opening your amrfinder output file. The following error was reported."
            )
            raise SystemExit(1)
