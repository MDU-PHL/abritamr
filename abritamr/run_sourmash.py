import subprocess
# from abritamr.logger import log
# from abritamr.utils import check_sourmash


import sourmash
import screed


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


def run_sourmash_search(query_filename:str, SBT_filename:str, sid:str):
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
    # Load the sourmash index
    tree = load_sourmash_index(SBT_filename)
    query_sig = sourmash_sig(query_filename, sid)
    # # Perform the search
    results = tree.search(query_sig, threshold=0.1)

    for similarity, found_sig, filename in tree.search(query_sig, threshold=0.01):
        qname = query_sig
        sp = f"{found_sig}"
        sim = similarity*1000
        print(f"Query: {qname}, Found: {' '.join(sp.split('_'))}, Similarity: {sim}")

    return sp
    # return results
