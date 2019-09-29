import pathlib, pandas, math, sys, toml, re


class Collate:

    """
    a base class for collation of amrfinder results - to be used when not doing MDU QC
    """

    REFGENES = pathlib.Path(__file__).parent.parent / "db" / "refgenes_latest.csv"
    MATCH = ["ALLELEX", "BLASTX", "EXACTX"]
    NONRTM = [
        "Amikacin/Gentamicin/Kanamycin/Tobramycin",
        "Amikacin/Kanamycin",
        "Amikacin/Kanamycin/Tobramycin",
        "Amikacin/Quinolone",
        "Amikacin/Tobramycin",
        "Aminoglycosides",
        "Gentamicin",
        "Gentamicin/Kanamycin/Tobramycin",
        "Gentamicin/Tobramcyin",
        "Kanamycin",
        "Kanamycin/Tobramycin",
        "Spectinomycin",
        "Streptogramin",
        "Streptomycin",
        "Streptomycin/Spectinomycin",
        "Tobramycin",
    ]
    MACROLIDES = [
        "Erythromycin",
        "Erythromycin/Telithromycin",
        "Lincosamides",
        "Lincosamides/Streptogramin",
        "Macrolides",
        "Streptogramin",
    ]
    RTM = "Other aminoglycoside resistance (non-RMT)"
    MAC = "Macrolide, lincosamide & streptogramin resistance"

    
    def joins(self, dict_for_joining):
        """
        make them a comma separated list
        """
        for i in dict_for_joining:
            if i != "Isolate":
                dict_for_joining[i] = list(set(dict_for_joining[i]))
                dict_for_joining[i] = ",".join(dict_for_joining[i])

        return dict_for_joining

    def get_drugclass(self, reftab, row, colname):

        """
        if the enhanced subclass is in either NONRTM or MACROLIDES then then use the groups specified by Norelle. If it is empty (-) then fall back on the AMRFinder subclass, else report the extended subclass
        """
        gene_id_col = "Gene symbol" if colname != "refseq protein" else "Accession of closest sequence"
        d = reftab[reftab[colname] == row[1][gene_id_col]][
            "enhanced_subclass"
        ].unique()[0]

        if d in self.NONRTM:
            return self.RTM
        elif d in self.MACROLIDES:
            return self.MAC
        elif d == "-":
            return reftab[reftab[colname] == row[1]["Gene symbol"]][
                "subclass"
            ].unique()[0]
        else:
            return d

    def extract_bifunctional_name(self, protein, reftab):
        """
        extract the joint name of bifunctional genes
        """
        return reftab[reftab["refseq protein"] == protein]["gene family"].values[0]
    
    def setup_dict(self, drugclass_dict, reftab, row):
        """
        return the dictionary for collation
        """

        if row[1]["Gene symbol"] in list(reftab["#allele"]):

            drugclass = self.get_drugclass(
                    reftab=reftab, row=row, colname="#allele"
                    )
            drugname = row[1]["Gene symbol"]
        elif row[1]["Gene symbol"] in list(reftab["gene family"]):
            drugclass = self.get_drugclass(
                reftab=reftab, row=row, colname="gene family"
            )
            drugname = row[1]["Gene symbol"]
        elif row[1]["Accession of closest sequence"] in list(reftab["refseq protein"]):
            drugclass = self.get_drugclass(
                reftab = reftab, row = row, colname = "refseq protein"
            )
            drugname = self.extract_bifunctional_name(protein = row[1]["Accession of closest sequence"], reftab = reftab)
        else:
            drugclass = "Unknown"
        if drugclass in drugclass_dict:
            drugclass_dict[drugclass].append(drugname)
        else:
            drugclass_dict[drugclass] = [drugname]

        return drugclass_dict

    def get_per_isolate(self, reftab, df, isolate):
        """
        make two dictionaries for each isolate that contain the drug class assignments for each match that is one of ALLELEX, EXACTX or BLASTX and then another dictionary which lists all partial mathces
        """
        drugclass_dict = {"Isolate": isolate}
        partials = {"Isolate": isolate}
        
        for row in df.iterrows():
            
            # if the match is good then generate a drugclass dict
            if row[1]["Method"] in self.MATCH and row[1]["Scope"] == "core" and row[1]["Element subtype"] == "AMR":
                drugclass_dict = self.setup_dict(drugclass_dict = drugclass_dict, reftab = reftab, row = row)
            elif row[1]["Method"] not in self.MATCH and row[1]["Scope"] == "core" and row[1]["Element subtype"] == "AMR":
                partials = self.setup_dict(drugclass_dict = partials, reftab = reftab, row = row)
        
        drugclass_dict = self.joins(dict_for_joining=drugclass_dict)
        partials = self.joins(dict_for_joining=partials)
        
        return drugclass_dict, partials

    def get_amr_output(self, path):
        """
        check that AMR finder output is present in the correct format.
        """
        amr_output = sorted(path.glob("*/*.out"))
        if amr_output != []:
            return amr_output
        else:
            print(
                "You do not have any AMR finder output. Please check your settings and try again."
            )
            raise SystemExit

    def save_files(self, path, tosave):

        tosave.to_csv(path)

    def collate(self):
        """
        if the refgenes.csv is present then proceed to collate data and save the csv files.
        """
        if self.REFGENES.exists():
            reftab = pandas.read_csv(self.REFGENES)
            reftab = reftab.fillna("-")
            # print(reftab.head())
            p = pathlib.Path.cwd()
            amr_output = self.get_amr_output(path=p)
            summary_drugs = pandas.DataFrame()
            summary_partial = pandas.DataFrame()
            for a in amr_output:
                
                df = pandas.read_csv(a, sep="\t")
                isolate = f"{a.parts[-2]}"
                
                drug, partial = self.get_per_isolate(
                    reftab=reftab, df=df, isolate=isolate
                )
                temp_drug = pandas.DataFrame(drug, index=[0])
                
                temp_partial = pandas.DataFrame(partial, index=[0])
                if summary_drugs.empty:
                    summary_drugs = temp_drug
                else:
                    summary_drugs = summary_drugs.append(temp_drug)

                if summary_partial.empty:
                    summary_partial = temp_partial
                else:
                    summary_partial = summary_partial.append(temp_partial)
            summary_drugs = summary_drugs.set_index("Isolate")
            summary_partial = summary_partial.set_index("Isolate")
            

            return summary_drugs, summary_partial
        else:
            print("The refgenes DB seems to be missing.")
            raise SystemExit

    def run(self):

        summary_drugs, summary_partial = self.collate()

        self.save_files(path="summary_matches.csv", tosave=summary_drugs)
        self.save_files(path="summary_partials.csv", tosave=summary_partial)


class MduCollate(Collate):
    def __init__(self):
        self.workdir = pathlib.Path.cwd()
        self.mduqc = self.workdir / "mdu_qc_checked.csv"
        self.check_for_mduqc()

    def check_for_mduqc(self):

        if self.mduqc.exists():
            return True
        else:
            print(
                f"It appears you are running mdu-amr in the context of MDU QC, the mdu_qc_checked.csv file is not present in {self.workdir}. Please check your settings and try again."
            )
            raise SystemExit

    def get_passed_isolates(self, qc_name):
        """
        generate lists of isolates that passed QC and need AMR, failed QC and should have AMR and all other isolates
        """
        

        qc = pandas.read_csv(qc_name, sep=None, engine="python")
        qc = qc.rename(columns={qc.columns[0]: "ISOLATE"})
        
        failed = list(
            qc[qc["TEST_QC"] == "FAIL"]["ISOLATE"]
        )  # isolates that failed qc and should have amr
        passed = list(
            qc[qc["TEST_QC"] == "PASS"]["ISOLATE"]
        )  # isolates that failed qc and should have amr
        
        return (passed, failed)

    def split_dfs(
        self, passed, failed, summary_drugs, summary_partial
    ):

        passed_match = summary_drugs[summary_drugs.index.isin(passed)]
        passed_partials = summary_partial[summary_partial.index.isin(passed)]
        failed_match = summary_drugs[summary_drugs.index.isin(failed)]
        failed_partials = summary_partial[summary_partial.index.isin(failed)]
    
        return (
            passed_match,
            passed_partials,
            failed_match,
            failed_partials
        )

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

    def strip_bla(self, gene):
        '''
        strip bla from front of genes except
        '''
        if gene.startswith("bla") and gene != "blaZ":
            gene = gene.replace("bla", "")
        return gene

    def reporting_logic(self, row, species):
        # get all genes found
        all_genes = self.get_all_genes(row)
        # determine the genus EXPECTED
        genus = species.split()[0]
        # A list of reportable genes - TODO move to a class variable
        reportable = [
            "Carbapenemase",
            "Carbapenemase (MBL)",
            "Carbapenemase (OXA-51 family)",
            "ESBL",
            "ESBL (Amp C type)",
            "Aminoglycosides (Ribosomal methyltransferases",
            "Colistin",
            "Oxazolidinone & phenicol resistance",
            "Vancomycin",
            "Methicillin",
        ]
        van_match = re.compile("van[A,B,C,D,E,G,L,M,N].[^\S]*")
        mec_match = re.compile("mec[^IR]")
        genes_reported = []  # genes for reporting
        genes_not_reported = []  # genes found but not reportable
        for r in reportable:
            try:
                genes = row[1][r].split(",")
                for gene in genes:
                    if r == "Carbapenemase (OXA-51 family)" and species not in [
                        "Acinetobacter baumannii",
                        "Acinetobacter calcoaceticus",
                        "Acinetobacter nosocomialis",
                        "Acinetobacter pittii",
                        "Acinetobacter baumannii complex",
                    ]:
                        genes_reported.append(gene)
                    elif (
                        r == "Carbapenemase (MBL)"
                        and species != "Stenotrophomonas maltophilia"
                    ):
                        genes_reported.append(gene)
                    elif (
                        r == "Carbapenemase (MBL)"
                        and species == "Stenotrophomonas maltophilia"
                        and gene != "blaL1"
                    ):
                        genes_reported.append(gene)
                    elif r in ["ESBL", "ESBL (Amp C type)"] and genus in [
                        "Salmonella",
                        "Shigella",
                    ]:
                        genes_reported.append(gene)
                    elif r == "Oxazolidinone & phenicol resistance" and (
                        genus == "Enterococcus"
                        or species
                        in ["Staphylococcus aureus", "Staphylococcus argenteus"]
                    ):
                        genes_reported.append(gene)
                    elif r == "Vancomycin" and gene.match(van_match):
                        genes_reported.append(gene)
                    elif r == "Methicilin" and gene.match(mec_match):
                        genes_reported.append(gene)
                    elif r in [
                        "Carbapenemase",
                        "Aminoglycosides (Ribosomal methyltransferases",
                        "Colistin"
                    ]:
                        genes_reported.append(gene)
                    else:
                        genes_not_reported.append(gene)

            except:
                print(f"{r} not found")

        genes_not_reported.extend([a for a in all_genes if a not in genes_reported])

        return genes_reported, genes_not_reported

    def mdu_reporting(self, match):

        reporting_df = pandas.DataFrame()
        qc = pandas.read_csv(self.mduqc, sep=None, engine="python")
        qc = qc.rename(columns={qc.columns[0]: "ISOLATE"})
        for row in match.iterrows():
            d = {"MDU sample ID": row[0]}
            qcdf = qc[qc["ISOLATE"] == row[0]]
            exp_species = qcdf["SPECIES_EXP"].values[0]
            obs_species = qcdf["SPECIES_OBS"].values[0]
            if qcdf["TEST_QC"].values[0] == 'FAIL':
                d["Species_exp"] = exp_species
            d["Species_obs"] = obs_species
            species = obs_species if obs_species == exp_species else exp_species
            genes_reported, genes_not_reported = self.reporting_logic(
                row=row, species=species
            )
            d["Item code"] = qcdf["ITEM_CODE"].values[0] 
            d["Resistance genes (alleles) detected"] = ",".join(genes_reported)
            d["Resistance genes (alleles) det (non-rpt)"] = ",".join(genes_not_reported)
            tempdf = pandas.DataFrame(d, index=[0])
            tempdf = tempdf.set_index("MDU sample ID")
            # print(tempdf)
            if reporting_df.empty:
                reporting_df = tempdf
            else:
                reporting_df = reporting_df.append(tempdf)
        return reporting_df

    def mdu_partials(self, partials):
        '''
        split the partial matches into reportable/nonreportable
        input is a dataframe with partial type as colnames 
        '''
        partial_names = ['PARTIALX', 'PARTIAL_CONTIG_ENDX']
        partial_df = pandas.DataFrame()
        for row in partials.iterrows:
            qc = pandas.read_csv(self.mduqc, sep=None, engine="python")
            qc = qc.rename(columns={qc.columns[0]: "ISOLATE"})
            exp_species = qcdf["SPECIES_EXP"].values[0]
            obs_species = qcdf["SPECIES_OBS"].values[0]
            species = obs_species if obs_species == exp_species else exp_species
            d = {"MDU sample ID": row[0]}
            
    def save_spreadsheet(
        self,
        passed_match,
        passed_partials,
        failed_match,
        failed_partials
    ):
        writer = pandas.ExcelWriter("MMS118.xlsx", engine="xlsxwriter")

        
        passed_match.to_excel(writer, sheet_name="MMS118")
        passed_partials.to_excel(writer, sheet_name="MMS118 - passed QC partial")
        failed_match.to_excel(writer, sheet_name="failed QC matches")
        failed_partials.to_excel(writer, sheet_name="failed QC partial")

        writer.close()

    def get_toml(self):
        '''
        get the toml file
        '''
        pass

    def write_to_toml(self):
        '''
        write the output of amrfinder to the toml
        '''
        pass

    def run(self):
        # get isolates binned into groups.
        for_amr, for_amr_failed = self.get_passed_isolates(self.mduqc)
        # get collated data
        summary_drugs, summary_partial = self.collate()
        # get dfs for specific groups
        passed_match, passed_partials, failed_match, failed_partials= self.split_dfs(
            for_amr, for_amr_failed, summary_drugs, summary_partial
        )
        # generate front tab for MDU reporting
        passed_match_df = self.mdu_reporting(match=passed_match)
        passed_partials_df = self.mdu_reporting(match = passed_partials)
        failed_match_df = self.mdu_reporting(match = failed_match)
        failed_partials_df = self.mdu_reporting(match = failed_partials)
        self.save_spreadsheet(
            passed_match_df,
            passed_partials_df,
            failed_match_df,
            failed_partials_df
        )


if __name__ == "__main__":

    if len(sys.argv) > 1:
        if sys.argv[1] == "mduqc":
            c = MduCollate()
    else:
        c = Collate()
    c.run()
