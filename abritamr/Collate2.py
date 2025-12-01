#!/usr/bin/env python3
import pathlib, math, sys,  re, logging
import pandas as pd
import numpy as np
import warnings
pd.options.mode.chained_assignment = None
from CustomLog import CustomFormatter


class Collate:

    def __init__(self,
                prefix: str,
                # species:str,
                input:str,
                refgenes:str
                ):

        self.logger =logging.getLogger(__name__) 
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(CustomFormatter())
        fh = logging.FileHandler('abritamr.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
        fh.setFormatter(formatter)
        self.logger.addHandler(ch) 
        self.logger.addHandler(fh)
        self.prefix = prefix
        # self.species = species
        self.input = input
        self.refgenes = self.get_reference(pth = refgenes)
        self.cfg = {
            "gene_name": "Element symbol",
            "default_name": "allele",
            "backup_name":"gene_family",
            "default_class": "Unknown",
            "closest_ref_accession": "Closest reference accession",
            "colnames":[
                'refseq_protein_accession',
                'genbank_protein_accession',
                'refseq_nucleotide_accession',
                'genbank_nucleotide_accession'],
            "output_delimiter": ",",
            "input_delimiter": "\t",
        }
        
    
    def get_reference(self, pth: str) -> pd.DataFrame:
        """
        Get reference genes as a dataframe
        """
        if pathlib.Path(pth).exists():
            refgenes = pd.read_csv(pth, sep = None, engine = "python")
            refgenes = refgenes.fillna("")
            return refgenes
        else:
            self.logger.critical(f"{pth} does not exist. You must have a reference gene collection. Please try again.")
            sys.exit(1)

      

    def get_drugclass(self, row_data:dict, dflt_subclass:str) -> str:
        """
        Identifies the enhanced drug class from the reference genes based on the input row data.

        Parameters:
        row_data (dict): A dictionary representing a row of input data.
        dflt_subclass (str): The default subclass to return if no match is found.

        Returns:
        str: The identified enhanced drug class, or the default subclass if no match is found.

        Raises:
        SystemExit: If the reference gene collection does not exist or if multiple drug classes are found for a single accession.
        """
        enhanced_drugclass = self.cfg["default_subclass"]
        accession = row_data[self.cfg["closest_ref_accession"]]
        for col in self.cfg["colnames"]:
            print(col)
            if accession in list(self.refgenes[col]):
                enhanced_drugclass = self.refgenes[self.refgenes[col] == accession]["enhanced_subclass"].unique()
                if len(enhanced_drugclass) > 1:
                    self.logger.warning(f"{accession} has multiple drug classes. This is not expected. Please check your reference gene collection.")
                    sys.exit(1)
                
                return enhanced_drugclass[0]
            
    def generate_dict(self, tab:pd.DataFrame, smp:str) -> dict:
        """
        Generates a dictionary from the given DataFrame, categorizing genes by their enhanced drug class.

        Parameters:
        tab (pd.DataFrame): A DataFrame containing the input data.
        smp (str): A sample identifier to be included in the dictionary.

        Returns:
        dict: A dictionary where keys are drug classes and values are lists of gene symbols associated with those classes.
        """
        dct = {"ID":smp}
        
        for _, row in tab.iterrows():
            row_data = row.to_dict()
            
            drugclass = self.get_drugclass(row_data = row_data, dflt_subclass = "Unknown")
            if drugclass not in dct:
                dct[drugclass] = [row_data[self.cfg["gene_name"]]]
            else:
                dct[drugclass].append(row_data[self.cfg["gene_name"]])
        
        
        return dct

    def collate(self):

        tab = pd.read_csv(self.input, sep = "\t")
        raw_dict = self.generate_dict(tab = tab, smp = "test")

        

fndr = sys.argv[1]

C = Collate(
    prefix= "test",
                # species:str,
                input = fndr,
                refgenes = "/home/khhor/dev/abritamr/abritamr/db/refgenes_latest.csv"
)

C.collate()