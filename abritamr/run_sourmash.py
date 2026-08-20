import subprocess
from abritamr.logger import log
from abritamr.utils import check_sourmash



def make_sketch(assembly_path:str, sketch_path:str, ksize:int = 31, scaled:int = 1000) -> None:
    """
    Make a sourmash sketch of the assembly.

    Parameters
    ----------
    assembly_path : str
        Path to the assembly file.
    sketch_path : str
        Path to the output sketch file.
    ksize : int, optional
        K-mer size for the sketch. Default is 31.
    scaled : int, optional
        Scaled value for the sketch. Default is 1000.

    Returns
    -------
    sketch_path : str
        Path to the output sketch file.
    """
    log.info(f"Making sourmash sketch of {assembly_path}.")
    cmd = f"sourmash sketch -k {ksize} --scaled {scaled} -o {sketch_path} {assembly_path}"
    proc = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
    
    if proc.returncode != 0:
        log.error(f"Error making sourmash sketch: {proc.stderr}")
        raise RuntimeError(f"Error making sourmash sketch: {proc.stderr}")
    
    log.info(f"Sourmash sketch saved to {sketch_path}.")


    return sketch_path


def compare_sketches(query_sketch:str, subject_sketch:str) -> None:
    """
    Compare two sourmash sketches.

    Parameters
    ----------
    query_sketch : str
        Path to the query sketch file.
    subject_sketch : str
        Path to the subject sketch file.
    
    Returns
    -------
    None
    """
    log.info(f"Comparing sourmash sketches: {query_sketch} vs {subject_sketch}.")
    cmd = f"sourmash search {query_sketch} {subject_sketch}"
    proc = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
    
    if proc.returncode != 0:
        log.error(f"Error comparing sourmash sketches: {proc.stderr}")
        raise RuntimeError(f"Error comparing sourmash sketches: {proc.stderr}")
    else:
        log.info(f"Sourmash comparison completed.")