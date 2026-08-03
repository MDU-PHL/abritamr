import pandas as pd
import pathlib


def amrf2dict(amrfinder:str)->dict:

    if pathlib.Path(amrfinder).exists():
        try:
            df = pd.read_csv(amrfinder, sep = "\t")

            amr = df.to_dict(orient = "records")
            return amr

        except Exception as e:
            log.critical(f"Something has gone wrong opening your amrfinder output file. The following error was reported.")
            raise SystemExit(1)