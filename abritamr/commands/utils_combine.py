import pathlib
import pandas as pd

from abritamr.logger import log


def save(
    df: pd.DataFrame, suffix: str, output_directory: str, prefix: str = "abritamr"
):

    if not df.empty:
        oname = f"{prefix}_{suffix}.csv"

        try:
            df.to_csv(f"{oname}", index=False)
            log.info(f"Saved {oname}")
        except Exception as e:
            log.warning(
                f"Could not save {oname}. The foloowing error was reported: {e}"
            )


def combine_files(args):

    output_suffixes = ["scan", "linelist", "matrix", "typed", "infer"]

    for s in output_suffixes:
        files = sorted(pathlib.Path(f"{args.input_directory}").rglob(f"*_{s}.csv"))
        if files != []:
            dfs = []
            log.info(f"Found files for {s}. Will now combine and save.")
            for f in files:
                df = pd.read_csv(f)
                dfs.append(df)
            save(
                df=pd.concat(dfs),
                suffix=s,
                output_directory=args.output_directory,
                prefix=args.prefix,
            )
        else:
            log.info(f"No files found for {s}")
