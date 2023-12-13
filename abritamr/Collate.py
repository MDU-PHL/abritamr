#!/usr/bin/env python3
import pathlib, pandas, math, sys,  re, logging, numpy
import warnings
pandas.options.mode.chained_assignment = None
# from pandas.core.algorithms import isin
from abritamr.CustomLog import CustomFormatter

class Collate:

    """
    a base class for collation of amrfinder results - to be used when not doing MDU QC
    """
    
    ANNOTATIONS = {'blast':'*','partial':'^','exact':''}
    REFGENES = pathlib.Path(__file__).parent / "db" / "refgenes_latest.csv"
    MATCH = ["ALLELEX", "BLASTX", "EXACTX", "POINTX"]

    def __init__(self, args):
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
        self.prefix = args.prefix
        self.run_type = args.run_type
        self.input = args.input

    def joins(self, dict_for_joining):
        """
        make them a comma separated list
        """
        for i in dict_for_joining:
            if i != "Isolate":
                dict_for_joining[i] = sorted(list(set(dict_for_joining[i])))
                dict_for_joining[i] = ",".join(dict_for_joining[i])

        return dict_for_joining

    def get_drugclass(self, reftab, row, colname):

        """
        if the enhanced subclass is in either NONRTM or MACROLIDES then then use the groups specified by Norelle. If it is empty (-) then fall back on the AMRFinder subclass, else report the extended subclass
        """
        gene_id_col = "Gene symbol" if colname != "refseq_protein_accession" else "Accession of closest sequence" # to get the name of drug and the drugclass
        
        if reftab[reftab[colname] == row[1][gene_id_col]].empty:
            d = 'Unknown'
            for i in ['genbank_protein_accession','refseq_nucleotide_accession','genbank_nucleotide_accession']:
                # print(reftab[reftab[i] == row[1][gene_id_col]])
                if not reftab[reftab[i] == row[1][gene_id_col]].empty:
                    d= reftab[reftab[i] == row[1][gene_id_col]]['enhanced_subclass'].values[0]
                    break
        else:
            
            d = reftab[reftab[colname] == row[1][gene_id_col]]['enhanced_subclass'].values[0]
        
        return d

    def extract_bifunctional_name(self, protein, reftab):
        """
        extract the joint name of bifunctional genes
        """
        return reftab[reftab["refseq_protein_accession"] == protein]["gene_family"].values[0]

    def extract_gene_name(self, protein, reftab, pointn = False):
        
        # add in to address the calling of point mutations that are not an exact match to reference
        # amrfinderplus returns the nucleotide accession with the range for what was detected.. these
        # so extract just the name of accession and set col to work with to nucleotide rather than protein
        # Add POINTN to end so users know this is a different type of match 
        suff = f"_POINTN" if pointn else ""
        protein = protein.split(':')[0] if pointn else protein
        col ="refseq_nucleotide_accession" if pointn else "refseq_protein_accession"
        
        try:
            if reftab[reftab[col] == protein]["allele"].values[0] != '-':
                nme = f"{reftab[reftab[col] == protein]['allele'].values[0]}{suff}"
                return nme
            else:
                return reftab[reftab[col] == protein]["gene_family"].values[0]
        except IndexError:
            for i in ['genbank_protein_accession','refseq_nucleotide_accession','genbank_nucleotide_accession']:
                if len(reftab[reftab[i] == protein]["gene_family"].unique()) >= 1:
                    return reftab[reftab[i] == protein]["gene_family"].unique()[0]
            
    def setup_dict(self, drugclass_dict, reftab, row, _type = 'exact', pointn = False):
        """
        return the dictionary for collation
        """
        if row[1]["Gene symbol"] in list(reftab["allele"]) and 'POINT' not in row[1]['Method']:
            drugclass = self.get_drugclass(
                    reftab=reftab, row=row, colname="allele"
                    )
            drugname = self.extract_gene_name(protein = row[1]["Accession of closest sequence"], reftab = reftab, pointn = pointn)
        elif row[1]["Gene symbol"] in list(reftab["allele"]) and 'POINT' in row[1]['Method']:
            drugclass = self.get_drugclass(
                    reftab=reftab, row=row, colname="allele"
                    )
            drugname = row[1]["Gene symbol"]
            
        elif row[1]["Gene symbol"] in list(reftab["gene_family"]):
            
            drugclass = self.get_drugclass(
                reftab=reftab, row=row, colname="refseq_protein_accession"
            )
            drugname = f"{self.extract_gene_name(protein = row[1]['Accession of closest sequence'], reftab = reftab)}*" if not row[1]["Method"] in ["EXACTX", "ALLELEX"] else f"{self.extract_gene_name(protein = row[1]['Accession of closest sequence'], reftab = reftab)}"
            
        elif row[1]["Accession of closest sequence"] in list(reftab["refseq_protein_accession"]):
            
            drugclass = self.get_drugclass(
                reftab = reftab, row = row, colname = "refseq_protein_accession"
            )
            drugname = self.extract_bifunctional_name(protein = row[1]['Accession of closest sequence'], reftab = reftab)
        else:
            drugname = row[1]["Gene symbol"]
            drugclass = "Unknown"

        if drugclass in drugclass_dict:
            drugclass_dict[drugclass].append(drugname)
        elif drugclass not in drugclass_dict:
            drugclass_dict[drugclass] = [drugname]

        return drugclass_dict

    def _other_dict(self, other_dict, row):
        """
        report virulence and stress genes
        """
        if row[1]['Element subtype'].capitalize() in other_dict:
            other_dict[row[1]['Element subtype'].capitalize()].append(row[1]['Gene symbol'])
        else:
            other_dict[row[1]['Element subtype'].capitalize()] = [row[1]['Gene symbol']]
        return other_dict

    def get_per_isolate(self, reftab, df, isolate):
        """
        make three dictionaries for each isolate that contain the drug class assignments for each match that is one of ALLELEX,POINTX, EXACTX or BLASTX, another dictionary which lists all partial mathces and a dictionary of virulence factors
        """
        drugclass_dict = {"Isolate": isolate}
        partials = {"Isolate": isolate}
        other = {"Isolate": isolate}
        for row in df.iterrows():
            # if the match is good then generate a drugclass dict
            if row[1]["Gene symbol"] == "aac(6')-Ib-cr" and row[1]["Method"] in ["EXACTX", "ALLELEX"]: # This is always a partial - unclear
                partials = self.setup_dict(drugclass_dict = partials, reftab = reftab, row = row)
            elif row[1]["Method"] in self.MATCH and row[1]["Element type"] == "AMR" and row[1]['Element subtype'] != "AMR-SUSCEPTIBLE":
                drugclass_dict = self.setup_dict(drugclass_dict = drugclass_dict, reftab = reftab, row = row)
            elif "POINTN" in row[1]["Method"] and row[1]["Element type"] == "AMR" and row[1]['Element subtype'] != "AMR-SUSCEPTIBLE":
                drugclass_dict = self.setup_dict(drugclass_dict = drugclass_dict, reftab = reftab, row = row, pointn = True)
            elif row[1]["Method"] not in self.MATCH and row[1]["Element type"] == "AMR" and row[1]['Element subtype'] != "AMR-SUSCEPTIBLE":
                partials = self.setup_dict(drugclass_dict = partials, reftab = reftab, row = row)
            else:
                other = self._other_dict(other_dict = other, row = row)
        drugclass_dict = self.joins(dict_for_joining=drugclass_dict)
        partials = self.joins(dict_for_joining=partials)
        other = self.joins(dict_for_joining = other)
        return drugclass_dict, partials, other

    def _get_cols(self, df, _not = 'Isolate'):
        
        cols = [c for c in df.columns if c != _not]

        return cols

    def _add_caret(self,df,cols):

        for c in cols:
            df[c] = numpy.where(df[c]!="", df[c].str.strip("*") + "^", df[c])
        return df

    def _merge(self, match, partial):

        df = pandas.DataFrame()
        tmp = match.merge(partial, on = 'Isolate', how = 'outer')
        tmp_cols = self._get_cols(df = tmp)
        
        for t in tmp_cols:
            if t.endswith('_x') or t.endswith('_y'): 
                # if this is a a col that has appeared in more than one df
                if df.empty:
                    df = tmp[['Isolate']]
                c = t.split('_')[0] # strip the _x or _y for col name
                if c not in list(df.columns): # if this has not already been added - should never happen but best safe than sorry
                    df[c] = tmp[[f"{c}_x", f"{c}_y"]].apply(lambda x: ','.join([i for i in x if i != '']), axis = 1) # combine the values and join with a comma
            else:
                # if this is a column that only appears in one df
                if df.empty:
                    df = tmp[['Isolate', t]]
                else:
                    df = df.merge(tmp[['Isolate', t]], on = 'Isolate',how = 'outer') 
        
        return df if not df.empty else tmp

    def _combine_dfs(self, match, partial, virulence):
        
        match = match.fillna('')
        partial = partial.fillna('')
        pcols = self._get_cols(df = partial)
        mcols = self._get_cols(df = match)
        vcols = self._get_cols(df = virulence)
        
        if vcols == [] and pcols == [] and mcols == []:
            return pandas.DataFrame()

        if pcols != []:
            partial = self._add_caret(df = partial, cols = pcols)      
        
        df = self._merge(match = match, partial = partial)
        
        if not isinstance(virulence, pandas.DataFrame) and not df.empty:
            df = df.merge(virulence, on = 'Isolate', how = 'outer') 
        elif isinstance(virulence, pandas.DataFrame) and df.empty:
            if not virulence.empty:
                df = virulence
            else:
                df = pandas.DataFrame()

        return df

    def save_files(self, path, match, partial, virulence):
        
        files = {'summary_matches.txt': match, 'summary_partials.txt': partial, 'summary_virulence.txt':virulence}
        
        for f in files:
            out = f"{path}/{f}" if path != '' else f"{f}"
            self.logger.info(f"Saving {out}")
            files[f].set_index('Isolate').to_csv(f"{out}", sep = '\t')
        combd = self._combine_dfs(match = match, partial = partial, virulence = virulence)
        combd_out = f"{path}/abritamr.txt" if path != '' else f"abritamr.txt"
        if not combd.empty:
            self.logger.info(f"Saving combined file : {combd_out}")
            combd.set_index('Isolate').to_csv(f"{combd_out}", sep = '\t')
        
        return True
        

    def _get_reftab(self):
        """
        get reftab
        """

        reftab = pandas.read_csv(self.REFGENES)
        reftab = reftab.fillna("-")

        return reftab

    def collate(self, prefix = ''):
        """
        if the refgenes.csv is present then proceed to collate data and save the csv files.
        """

        
        reftab = self._get_reftab()
        
        df = pandas.read_csv(f"{prefix}/amrfinder.out", sep="\t")
        self.logger.info(f"Opened amrfinder output for {prefix}")
        drug, partial, virulence = self.get_per_isolate(
            reftab=reftab, df=df, isolate=prefix
        )
        
        summary_drugs = pandas.DataFrame(drug, index = [0])
        summary_partial = pandas.DataFrame(partial, index = [0])
        summary_virulence = pandas.DataFrame(virulence, index = [0])
        return summary_drugs, summary_partial,summary_virulence
        
    
    def _combine_df(self, existing, temp):
        """
        combine result dataframes for batch
        """
        if existing.empty:
            existing = temp
        else:
            existing = pandas.concat([existing,temp])
        
        return existing

    def _batch_collate(self,input_file):


        summary_matches = pandas.DataFrame()
        summary_partial = pandas.DataFrame()
        summary_virulence = pandas.DataFrame()

        df = pandas.read_csv(input_file, sep = '\t', header = None)
        for row in df.iterrows():
            prefix = f"{row[1][0]}"
            self.logger.info(f"Collating results for {prefix}")
            temp_match, temp_partial, temp_virulence = self.collate(prefix = prefix)
            summary_matches = self._combine_df(existing = summary_matches, temp = temp_match)
            summary_partial = self._combine_df(existing = summary_partial, temp = temp_partial)
            summary_virulence = self._combine_df(existing = summary_virulence, temp = temp_virulence)
        
        return summary_matches, summary_partial, summary_virulence

    def run(self):


        if not pathlib.Path(self.REFGENES).exists():
            self.logger.critical(f"The refgenes DB ({self.REFGENES}) seems to be missing.")
            raise SystemExit

        if self.run_type != 'batch':
            self.logger.info(f"This is a single sample run.")
            summary_drugs, summary_partial, virulence = self.collate(prefix = self.prefix)
        else:
            self.logger.info(f"You are running abritamr in batch mode. Your collated results will be saved.")
            summary_drugs, summary_partial, virulence = self._batch_collate(input_file = self.input)
        self.logger.info(f"Saving files now.")
        self.save_files(path='' if self.run_type == 'batch' else f"{self.prefix}", match = summary_drugs,partial=summary_partial, virulence = virulence)
        
class MduCollate(Collate):
    
    def __init__(self, args):
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
        self.sop = args.sop
        self.sop_name = args.sop_name
        self.mduqc = args.qc
        self.db = args.db
        self.partials = args.partials
        self.match = args.matches
        self.runid = args.runid
        self.NONE_CODES = {
            "Salmonella":"CPase_ESBL_AmpC_16S_NEG",
            "Shigella":"CPase_ESBL_AmpC_16S_NEG",
            "Staphylococcus":"Mec_VanAB_Linez_NEG",
            "Enterococcus":"Van_Linez_NEG",
            "Other":"CPase_16S_mcr_NEG"
        }
        self.REPORTING = {"Salmonella enterica":self.mdu_reporting_salmonella}
    def mdu_qc_tab(self):
        self.logger.info(f"Checking the format of the QC file")
        cols = ["ISOLATE", 'SPECIES_EXP', 'SPECIES_OBS', 'TEST_QC']
        tab = pandas.read_csv(self.mduqc)
        tab = tab.rename(columns = {tab.columns[0]: 'ISOLATE'})
        for c in cols:
            if c not in list(tab.columns):
                self.logger.critical(f"There seems to be a problem with your QC file. This file must have {','.join(cols)}. Please check your input and try again.")
                raise SystemExit

        pos = pandas.DataFrame(data = {"ISOLATE": "9999-99888", "TEST_QC" : "PASS", "SPECIES_EXP":"Staphylococcus aureus", "SPECIES_OBS":"Staphylococcus aureus" }, index = [0])
        
        return pandas.concat([tab,pos])

    def strip_bla(self, gene):
        '''
        strip bla from front of genes except
        '''
        if gene.startswith("bla") and len(gene) >6 and gene.endswith("*"):
            gene = gene.replace("bla", "")
        elif gene.startswith("bla") and len(gene) >5 and not gene.endswith("*"):
            gene = gene.replace("bla", "")
        return gene


    # Carbapenemase to be reported for all cases
    # Carbapenemase (MBL) all in all HOWEVER if blaL1 should  not be reported in Stenotrophomonas maltophilia
    # Carbapenemase (OXA-51 family) REPORTED IN ALL except in Acinetobacter baumannii,Acinetobacter calcoaceticus,Acinetobacter nosocomialis,Acinetobacter pittii,Acinetobacter baumannii complex,
    # All ESBL in Salmonella and Shigella
    # ESBL (Amp C type) in Salmonella and Shigella
    # Aminoglycosides (Ribosomal methyltransferases) ALL
    # Colistin ALL
    # Oxazolidinone & phenicol resistance if Genus = Enterococcus or Staphylococcus aureus and Staphylococcus argenteus
    # report vanA*, B*, C*, vanD*, vanE*, vanG*, vanL*, vanM*, vanN* 
    # Methicillin ALL

    def get_all_genes(self, row):
        all_genes = []
        for r in row[1]:
            if isinstance(r, str):
                if len(r.split(",")) > 1:
                    for a in r.split(","):
                        all_genes.append(a)
                else:
                    all_genes.append(r)
        return all_genes

    def none_replacement_code(self, genus):

        if genus in self.NONE_CODES:
            return self.NONE_CODES[genus]
        else:
            return "CPase_16S_mcr_NEG"

    def assign_itemcode(self,mduid, reg):
        self.logger.info(f"Checking for item code")
        m = reg.match(mduid)
        try:
            itemcode = m.group('itemcode') if m.group('itemcode') else ''
        except AttributeError:
            itemcode = ''
        return itemcode

    def assign_mduid(self, mduid, reg):
        self.logger.info(f"Extracting MDU sample ID")
        m = reg.match(mduid)
        try:
            mduid = m.group('id')
        except AttributeError:
            mduid = mduid.split('/')[-1]
        return mduid

    def _ampicillin_res_sal(self, col, gene):
        
        if col in [ 'Beta-lactam','ESBL','AmpC'] or 'Ampicillin' in col:
            return gene
        return ''

    def _chloramphenicol_res_sal(self, col, gene):

        if 'phenicol' in col.lower():
            return gene
        return ''

    def _cefo_esbl_res_sal(self, col, gene):

        if col == 'ESBL':
            return gene
        return ''

    def _cefo_ampc_res_sal(self, col, gene):

        if 'AmpC' in col :
            return gene
        return ''
    
    def _carbapenem_res_salmo(self, col, gene):

        if "Carbapenemase" in col: #not KPC
            return gene
        return ''

    def _azi_res_salmo(self, col, gene):

        if 'Azithromycin' in  col or col in ['Macrolide', 'Lincosamide/Macrolide/Streptogramin']:
            return gene
        return ''

    def _gentamicin_res_salm(self, col, gene):
        if col in ['Gentamicin', 'Aminoglycosides (Ribosomal methyltransferase)']: 
            return gene
        return ''
    
    def _kanamycin_res_salm(self, col, gene):
    
        if col in ['Kanamycin','Aminoglycosides (Ribosomal methyltransferase)']:
            return gene
        return ''
    
    def _streptomycin_res_salm(self, col, gene):
        if 'Streptomycin' in col:
            return gene
        return ''
    
    def _spectinomycin_res_salm(self, col, gene):
        if 'Streptomycin' in col:
            return gene
        return ''
    
    def _tetra_res_salmo(self, col, gene):

        if col == 'Tetracycline':
            return gene
        return ''

    def _cipro_res_salmo(self, col, gene):

        if 'Quinolone' in col:
            return gene
        return ''

    def _sulf_res_salmo(self, col, gene):
        
        if col == "Sulfonamide":
            return gene
        return ''
    
    def _trimet_res_salmo(self, col, gene):

        if col == "Trimethoprim":
            return gene
        return ''
    
    def _trim_sulpha_salmo(self, trim_gene, sul_gene):

        return f"{trim_gene},{sul_gene}"

    def _rmt_res_salmo(self, col, gene):

        if col == "Aminoglycosides (Ribosomal methyltransferase)":
            return gene
        return ''

    def _colistin_res_salmo(self, col, gene):
        if col == 'Colistin':
            return gene
        return ''
    
    def reporting_logic_salmonella(self, row):

        # Ampicillin - ResMech	Ampicillin - Interpretation	
        # Cefotaxime (ESBL) - ResMech	Cefotaxime (ESBL) - Interpretation	
        # Cefotaxime (AmpC) - ResMech	Cefotaxime (AmpC) - Interpretation	
        # Tetracycline - ResMech	Tetracycline - Interpretation	
        # Gentamycin - ResMech	Gentamycin - Interpretation	
        # Sulfathiazole - ResMech	Sulfathiazole - Interpretation	
        # Trimethoprim - ResMech	Trimethoprim - Interpretation	
        # Ciprofloxacin - ResMech	Ciprofloxacin - Interpretation	
        # Azithromycin - ResMech	Azithromycin - Interpretation

        mduidreg = re.compile(r'(?P<id>[0-9]{4}-[0-9]{5,6})-?(?P<itemcode>.{1,2})?')
        
        all_genes = self.get_all_genes(row)
        all_genes = [a for a in all_genes if a != row[1]['Isolate'] and a != '']
                
        isodict = row[1].to_dict()
        item_code = self.assign_itemcode(row[1]['Isolate'], mduidreg)
        md = self.assign_mduid(row[1]['Isolate'], mduidreg)

        abx = {
            "Ampicillin" : self._ampicillin_res_sal,
            "Cefotaxime (ESBL)":self._cefo_esbl_res_sal,
            "Cefotaxime (AmpC)":self._cefo_ampc_res_sal,
            "Meropenem" :self._carbapenem_res_salmo,
            "Azithromycin":self._azi_res_salmo,
            "Gentamicin":self._gentamicin_res_salm,
            "Kanamycin":self._kanamycin_res_salm,
            "Streptomycin":self._streptomycin_res_salm,
            "Spectinomycin":self._spectinomycin_res_salm,
            "Tetracycline":self._tetra_res_salmo,
            "Ciprofloxacin":self._cipro_res_salmo,
            "Sulfathiazole":self._sulf_res_salmo,
            "Trimethoprim":self._trimet_res_salmo,
            "Trim-Sulpha":self._trim_sulpha_salmo,
            "Chloramphenicol":self._chloramphenicol_res_sal,
            "Aminoglycosides (RMT)": self._rmt_res_salmo,
            "Colistin":self._colistin_res_salmo
        }

        
        tmp_results = {
            "Ampicillin":[],	
            "Cefotaxime (ESBL)":[],
            "Cefotaxime (AmpC)":[],
            "Tetracycline":[],
            "Gentamicin":[],
            "Sulfathiazole":[],
            "Trimethoprim":[],
            "Trim-Sulpha":[],
            "Kanamycin":[],
            "Streptomycin":[],
            "Spectinomycin":[],
            "Ciprofloxacin":[],
            "Azithromycin":[],
            "Meropenem":[],
            "Chloramphenicol":[],
            "Aminoglycosides (RMT)":[],"Colistin":[], "Other":[]
        }

        
        for ab in abx:
            if ab != 'Trim-Sulpha':
                for i in isodict:
                    if i != "Isolate":
                        
                        g = abx[ab](col = i, gene = isodict[i])
                        
                        if g != '':
                            gene_list = g.split(',')
                            for gene in gene_list:
                                tmp_results[ab].append(gene)
                                if gene in all_genes:
                                    all_genes.remove(gene)
            
        # interprete results
        tmp_results['Other'] = all_genes
        results = {'Isolate': row[1]['Isolate'], 'MDU Sample ID':md, 'Item code': item_code}
        # for Trim-Sulpha
        
        if tmp_results['Trimethoprim'] != [] and tmp_results['Sulfathiazole'] != []:
            tmp_results['Trim-Sulpha'] = list(set(tmp_results['Trimethoprim']).union(tmp_results['Sulfathiazole']))
        
        for res in tmp_results:
            results[f"{res} - ResMech"] = ';'.join(tmp_results[res]) if tmp_results[res] != [] else "None detected"              
            if res in ["Aminoglycosides (RMT)","Colistin", "Other"]:
                results[f"{res} - Interpretation"] = ''
            elif res == 'Ciprofloxacin':
                if tmp_results[res] == []:
                    results[f"{res} - Interpretation"] = 'Susceptible'
                elif len(tmp_results[res]) == 1:
                    results[f"{res} - Interpretation"] = 'Intermediate'
                else:
                    results[f"{res} - Interpretation"] = 'Resistant'
            elif res not in ["Aminoglycosides (RMT)","Colistin", "Other"]:
                if tmp_results[res] == []:
                    results[f"{res} - Interpretation"] = 'Susceptible'
                else:
                    results[f"{res} - Interpretation"] = 'Resistant'
        
        return results

            
    def reporting_logic_general(self, row, species, neg_code = True):
        # get all genes found
        all_genes = self.get_all_genes(row)
        isodict = row[1].to_dict()
        # determine the genus EXPECTED
        genus = species.split()[0]
        reportable = [
            "Carbapenemase",
            "Carbapenemase (MBL)",
            "ESBL (KPC variant)"
            "Carbapenemase (OXA-51 family)",
            "ESBL",
            "AmpC", # will need to change to reflect new subclass AmpC
            "Aminoglycosides (Ribosomal methyltransferase)",
            "Colistin",
            "Chloramphenicol/Florfenicol/Linezolid", 
            "Florfenicol/Oxazolidinone",
            "Vancomycin",
            "Methicillin"
        ]
        non_caveat_reportable = [
            "Carbapenemase",
            "ESBL (KPC variant)",
            "Aminoglycosides (Ribosomal methyltransferase)",
            "Colistin"
        ]
        # Carbapenemase (KPC variant)
# Carbapenemase (MBL)
# Carbapenemase (OXA-51 family)

        abacter_excluded = [
            "Acinetobacter baumannii",
            "Acinetobacter calcoaceticus",
            "Acinetobacter nosocomialis",
            "Acinetobacter pittii",
            "Acinetobacter baumannii complex"
        ]

        
        van_match = re.compile(r"van[A,B,C,D,E,G,L,M,N][\S]*")
        mec_match = re.compile(r"mec[^IR]")
        
        genes_reported = []  # genes for reporting
        genes_not_reported = []  # genes found but not reportable
        for i in isodict:
            genes = []
            if i != 'Isolate':
                if isinstance(isodict[i], str):
                    genes = isodict[i].split(',')
                    genes = [gene for gene in genes if '_' not in gene] # this is to remove the point mutations for MMS118
                if genes != []: # for each bin we do things to genes
                    if i in reportable: #don't report point mutations
                        # print(i)
                        if i in non_caveat_reportable:
                            genes_reported.extend(genes)
                        elif i == "Carbapenemase (MBL)" and species != "Stenotrophomonas maltophilia":
                            genes_reported.extend(genes)
                        elif i == "Carbapenemase (MBL)" and species == "Stenotrophomonas maltophilia":
                            # if species is "Stenotrophomonas maltophilia" don't report blaL1
                            genes_reported.extend([g for g in genes if not g.startswith("blaL1")])
                            genes_not_reported.extend([g for g in genes if g.startswith("blaL1")])
                        elif i == "Carbapenemase (OXA-51 family)" and species not in abacter_excluded:
                            genes_reported.extend(genes)
                        elif i in ["ESBL","AmpC"] and genus in ["Salmonella"]: 
                            genes_reported.extend(genes)
                        elif i in ["ESBL","AmpC"] and genus in ["Shigella"]: 
                            genes_reported.extend([g for g in genes if "blaEC" not in g])
                            genes_not_reported.extend([g for g in genes if "blaEC" in g]) # don't report blaEC for shigella
                        elif i == "Vancomycin":
                            genes_reported.extend([g for g in genes if van_match.match(g)])
                            genes_not_reported.extend([g for g in genes if not van_match.match(g)])
                        elif i == "Methicillin":
                            genes_reported.extend([g for g in genes if mec_match.match(g)])
                            genes_not_reported.extend([g for g in genes if not mec_match.match(g)])
                        
                        elif "Oxazolidinone" in i or "Linezolid" in i:
                            if species in ["Staphylococcus aureus","Staphylococcus argenteus"] or genus == "Enterococcus":
                                genes_reported.extend(genes)
                            else:
                                genes_not_reported.extend(genes)
                        else:
                            genes_not_reported.extend(genes)

                    else:
                        genes_not_reported.extend(genes)
        if genes_reported == []:
            genes_reported = [self.none_replacement_code(genus= genus)] if neg_code else ''
        if genes_not_reported == []:
            genes_not_reported = ["No non-reportable genes found."] if neg_code else ''
            # break
        
        self.logger.info(f"{row[1]['Isolate']} has {len(genes_reported)} reportable genes.")
        return genes_reported, genes_not_reported

    def mdu_reporting_salmonella(self, match, isolates):

        self.logger.info(f"Applying MDU business logic for interpretation of  Salmonella AST")
        cols = ["MDU Sample ID", "Item code", 
        "Ampicillin - ResMech", "Ampicillin - Interpretation",
        "Cefotaxime (ESBL) - ResMech","Cefotaxime (ESBL) - Interpretation",
        "Cefotaxime (AmpC) - ResMech","Cefotaxime (AmpC) - Interpretation",
        "Tetracycline - ResMech","Tetracycline - Interpretation",
        "Gentamicin - ResMech","Gentamicin - Interpretation",
        "Kanamycin - ResMech",	"Kanamycin - Interpretation",
        "Streptomycin - ResMech",	"Streptomycin - Interpretation",
        "Sulfathiazole - ResMech","Sulfathiazole - Interpretation",
        "Trimethoprim - ResMech","Trimethoprim - Interpretation",
        "Trim-Sulpha - ResMech",	"Trim-Sulpha - Interpretation",
        "Chloramphenicol - ResMech",	"Chloramphenicol - Interpretation",
        "Ciprofloxacin - ResMech","Ciprofloxacin - Interpretation",
        "Meropenem - ResMech",	"Meropenem - Interpretation",
        "Azithromycin - ResMech","Azithromycin - Interpretation",
        "Aminoglycosides (RMT) - ResMech",
        "Aminoglycosides (RMT) - Interpretation",
        "Colistin - ResMech",
        "Colistin - Interpretation", 
        "Other - ResMech",
        "Other - Interpretation"]
        # select passed Salmonella
        
        df = pandas.read_csv(match, sep = '\t')
        df = df[df['Isolate'].isin(isolates)]
        result_df = pandas.DataFrame()
        df = df.fillna('')
        for row in df.iterrows():
            isolate_results = self.reporting_logic_salmonella(row = row)
            tmpdf = pandas.DataFrame(isolate_results, index = [0])
            if result_df.empty:
                result_df = tmpdf
            else:
                result_df = result_df.append(tmpdf)
        
        return result_df[cols]

    def _extract_plus_isolates(self,species):

        qc = self.mdu_qc_tab()
        qc = qc[(qc['SPECIES_OBS'] == species) & (qc['TEST_QC'] == 'PASS')]

        return list(qc["ISOLATE"])


    def mdu_reporting_general(self, match, neg_code = True):

        self.logger.info(f"Applying MDU business logic {'matches' if neg_code else 'partials'}.")
        mduidreg = re.compile(r'(?P<id>[0-9]{4}-[0-9]{5,6})-?(?P<itemcode>.{1,2})?')
        reporting_df = pandas.DataFrame()
        qc = self.mdu_qc_tab()
        match_df = pandas.read_csv(match, sep = '\t')
        for row in match_df.iterrows():
            isolate = row[1]['Isolate']
            item_code = self.assign_itemcode(isolate, mduidreg)
            md = self.assign_mduid(isolate, mduidreg)
            d = {"MDU sample ID": md, "Item code" : item_code}
            qcdf = qc[qc['ISOLATE'].str.contains(isolate)]
            exp_species = qcdf["SPECIES_EXP"].values[0]
            obs_species = qcdf["SPECIES_OBS"].values[0]
           
            species = obs_species if obs_species == exp_species else exp_species
            genes_reported, genes_not_reported = self.reporting_logic_general(
                row=row, species=species, neg_code=neg_code
            )
            
            # strip bla
            genes_not_reported = [self.strip_bla(g) for g in genes_not_reported]
            genes_not_reported = [g for g in genes_not_reported if g != isolate]
            genes_reported = [self.strip_bla(g) for g in genes_reported]
            genes_reported = [g for g in genes_reported if g != isolate]
            d["Resistance genes (alleles) detected"] = ",".join(genes_reported)
            d["Resistance genes (alleles) det (non-rpt)"] = ",".join(genes_not_reported)
            if qcdf["TEST_QC"].values[0] == 'FAIL': # species not needed for MDU LIMS upload
                d["Species_exp"] = exp_species
            d["Species_obs"] = obs_species
            d["Species_exp"] = exp_species
            d['db_version'] = self.db
            tempdf = pandas.DataFrame(d, index=[0])
            tempdf = tempdf.set_index("MDU sample ID")
            
            if reporting_df.empty:
                reporting_df = tempdf
            else:
                reporting_df = pandas.concat([reporting_df,tempdf])
        
        return reporting_df.reindex(labels = ['Item code','Resistance genes (alleles) detected','Resistance genes (alleles) det (non-rpt)','Species_obs', 'Species_exp', 'db_version'], axis = 'columns')

    
    def save_spreadsheet_general(
        self,
        passed_match,
        passed_partials,
        
    ):
        self.logger.info(f"Saving {self.sop_name}.")
        writer = pandas.ExcelWriter(f"{self.runid}_{self.sop_name}.xlsx", engine="xlsxwriter")
        passed_match.to_excel(writer, sheet_name=f"{self.sop_name}")
        passed_partials.to_excel(writer, sheet_name="Passed QC partial")
        writer.close()

    def save_spreadsheet_interpreted(self, results):
        sheets = {"Salmonella enterica":f"{self.sop_name}-01"}
        self.logger.info(f"Saving MMS184")
        writer = pandas.ExcelWriter(f"{self.runid}_{self.sop_name}.xlsx", engine = "xlsxwriter")
        for result in results:
            result[1].to_excel(writer, sheet_name = sheets[result[0]], index = False)
        writer.close()

    def run(self):
        if self.sop == 'general' and pathlib.Path(self.partials).exists():
            passed_match_df = self.mdu_reporting_general(match=self.match)
            passed_partials_df = self.mdu_reporting_general(match = self.partials)
            self.save_spreadsheet_general(
                passed_match_df,
                passed_partials_df
            )
        elif self.sop == 'plus':
            dfs = []
            for r in self.REPORTING:
                self.logger.info(f"Checking if there are any {r} in this run.")
                isolates = self._extract_plus_isolates(species = r)
                if isolates != []:
                    self.logger.info(f"There are {len(isolates)} {r} in this run.")
                    plus_df = self.mdu_reporting_salmonella(match = self.match, isolates=isolates)
                    dfs.append((r,plus_df))
                else:
                    self.logger.info(f"There are no {r} in this run. Collation will be skipped.")
            self.save_spreadsheet_interpreted(results = dfs)
