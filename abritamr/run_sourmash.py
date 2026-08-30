

import sourmash
import screed
import tempfile
import pathlib  
import subprocess
from abritamr.logger import log


def load_sourmash_index(SBT_filename:str):
    """
    Load a sourmash index from a file.

    Parameters
    ----------
    SBT_filename : str
        Path to the sourmash index file.

    Returns
    -------
    sourmash.Index
        The loaded sourmash index.
    """
    # check_sourmash()    
    tree = sourmash.load_file_as_index(SBT_filename)
    return tree

def sourmash_sig(query_filename:str, sid:str):
    """
    Run a sourmash search using a query file and a sourmash index.

    Parameters
    ----------
    query_filename : str
        Path to the query file (e.g., FASTA or FASTQ).
    SBT_filename : str
        Path to the sourmash index file.

    Returns
    -------
    None
    """
    # check_sourmash()
    
    # Load the sourmash index
    minhash = sourmash.MinHash(ksize=31, n=0, scaled=10000) 
    query_seq = next(iter(screed.open(query_filename))).sequence
    minhash.add_sequence(query_seq)
    query_sig = sourmash.SourmashSignature(minhash, name=sid)

    return query_sig

def create_sourmash_index(tdir:str,SBT_path:str):
    """
    Create a sourmash index from a list of signatures.

    Parameters
    ----------
    SBT_filename : str
        Path to the output sourmash index file.
    sigs : list
        List of sourmash signatures.

    Returns
    -------
    None
    """
    # check_sourmash()
    cmd =  f"sourmash index --ksize 31 {tdir}/abritamrdb2 {SBT_path}/*.sig"
    proc = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
    if proc.returncode != 0:
        err = proc.stderr.split("\n")
        log.critical(f"The following error was reported:")
        log.critical("\n".join(err))
        raise SystemExit(1)
    return f"{tdir}/abritamrdb2.sbt.zip"


def run_sourmash_search(query_filename:str, SBT_path:str, sid:str):
    """
    Run a sourmash search using a query signature and a sourmash index.

    Parameters
    ----------
    query_filename : str
        Path to the query file (e.g., FASTA or FASTQ).
    SBT_filename : str
        Path to the sourmash index file.
    sid : str
        Sample ID for the query signature.

    Returns
    -------
    list of tuples
        A list of tuples containing the matching signatures and their similarity scores.
    """
    # check_sourmash()
    sp = ""
    # mx = 0
    # Load the sourmash index
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = pathlib.Path(temp_dir)
        
        tree = sourmash.load_file_as_index(create_sourmash_index(tdir = temp_dir_path, SBT_path = SBT_path))
        query_sig = sourmash_sig(query_filename, sid)
        # # Perform the search
        results = tree.search(query_sig, threshold=0.1)
        spc = []
        for similarity, found_sig, filename in tree.search(query_sig, threshold=0.0001):
            # qname = query_sig
            sp = f"{found_sig}"
            sim = similarity
            spc.append({sim: sp})
            # # print(f"Query: {qname}, Found: {' '.join(sp.split('_'))}, Similarity: {sim}")   
        try:
            mx = max([list(x.keys())[0] for x in spc])
            
            if mx >= 0.01:
                sp = [list(x.values())[0] for x in spc if list(x.keys())[0] == mx][0]
                sp = ' '.join(sp.split('_'))
                
            
        except ValueError:
            log.warning(f"No matches found for {sid} using sourmash. ")
            
    return sp
