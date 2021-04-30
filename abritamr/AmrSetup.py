import pathlib, pandas, datetime, getpass, jinja2, re, subprocess, os, logging
from abritamr.version import db

"""
A class for setting up mdu-amr
"""


class Setupamr(object):
    
    def __init__(self, args):
        # some variables to be use
        # create file handler which logs even debug messages
        # print(logging.__file__)
        self.species_list = ['Acinetobacter_baumannii', "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"]
        self.db = db
        self.logger =logging.getLogger(__name__) 
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('abritamr.log')
        fh.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)
        self.logger.addHandler(ch) 
        self.logger.addHandler(fh)
        # self.logger = logger
        self.workdir = pathlib.Path(args.workdir)
        self.resources = pathlib.Path(args.resources)
        self.jobs = args.jobs
        self.mduqc = args.mduqc
        self.positive_control = True if self.mduqc else args.positive_control
        # self.drugs = pathlib.Path(args.drug_classes)
        self.contigs = args.contigs
        self.from_contigs = True if args.contigs != '' else False
        self.prefix = args.prefix
        self.keep = args.keep
        # self.run_singulairty = args.Singularity
        # self.singularity_path = f"shub://phgenomics-singularity/amrfinderplus" if args.singularity_path == '' else args.singularity_path # this needs to addressed before formal release.
        self.finaloutput = (
            ['MMS118.xlsx']
            if self.mduqc
            else ['summary_matches.csv', 'summary_partials.csv']
        )
        self.qc = args.qc
        self.species = args.species if args.species in self.species_list else ""



    def file_present(self, name):
        """
        check file is present
        """
        
        if name == "":
            return False
        elif pathlib.Path(name).exists():
            self.logger.info(f"Checking if file {name} exists")
            return True
        else:
            return False

    def check_input_exists(self):
        """
        check that all files required are present. 
        If both contigs and amrfinder are absent (or both present) provide warning and exit.
        If mduqc is true ensure that the mdu_qc_checked file is present
        """
        self.logger.info(f"Checking that all correct input files are present.")
        if self.mduqc and not self.file_present("mdu_qc_checked.csv"):
            self.logger.warning(
                "You appear to be trying to run mdu-amr in the context of mdu-qc, but the mdu_qc_checked.csv file is not present. Please check your settings and try again."
            )
            raise SystemExit

        if self.mduqc and self.positive_control == "":
            self.logger.warning(
                "You appear to be trying to run mdu-amr in the context of mdu-qc, but you have not provided a path to positive control. Please check your settings and try again."
            )
            raise SystemExit

        if (
            not self.file_present(self.contigs)
            and not self.file_present(self.amrfinder_output)
            and not self.mduqc
        ):
            self.logger.warning(
                "You have not provided a valid path to any input files. Please provide a file containing paths to assemblies or amrfinder outputs."
            )
            raise SystemExit
        elif self.file_present(self.contigs) and self.file_present(
            self.amrfinder_output
        ):
            self.logger.warning(
                "You seem to have provided both assemblies and amrfinder results. Only one is required. Please check your setting and try again."
            )
            raise SystemExit

        if self.file_present(self.contigs):
            self.from_contigs = True
        elif self.file_present(self.amrfinder_output):
            self.from_contigs = False
        # print(self.from_contigs)
        self.logger.info(f"All files seem to be present and accounted for. Well done.")
        return True

    def make_links(self, first_column, second_column):
        """
        link files
        """
        target_name = "contigs.fa" if self.from_contigs else f"{first_column}.out"
        
        isolate_dir = self.workdir / f"{first_column}"
        
        if not isolate_dir.exists():
            isolate_dir.mkdir()
        source = pathlib.Path(f"{second_column}").absolute()
        
        target = isolate_dir / target_name
        # cmd = f"cp {source} {target}"
        # subprocess.run(cmd, shell = True)
        if not target.exists():
            target.symlink_to(source)

    def check_input_tab(self, tab):
        
        self.logger.info(f"Checking the structure of your input file.")
        if tab.shape[1] == 2:
            self.logger.info(f"The input file seems to be in the correct format. Thank you.")
            return 'batch'
        elif tab.shape[1] == 22:
            self.logger.info(f"You are running abriTAMR on individual samples. Have you provided a prefix for your output structure?")
            return 'amrfinder_output'
        elif tab.shape[1] == 1:
            self.logger.info(f"It seems you might be trying to run abriTAMR directly from an assembly file")
            return 'assembly'
        else:
            logging.warning(
                "Your input file should be a tab delimited file with two columns. Please check your input and try again."
            )
            raise SystemExit

    def link_input_files(self):
        """
        Ensure that the files (either contigs or amrfinder output) exist and generate structure for links
        """
        input_file = pathlib.Path(self.contigs)
        tab = pandas.read_csv(input_file, engine="python", header=None, sep = '\t')
        running_type = self.check_input_tab(tab)
        # print(running_type)
        if running_type == 'batch':
            self.logger.info(f"Checking that the input data is present. If present will link to {self.workdir}")
            for row in tab.iterrows():
                if self.file_present(row[1][1]):
                    self.make_links(first_column=row[1][0], second_column=row[1][1])
            if self.positive_control:
                self.logger.info(f"You have provided a path the a positive control. Checking if path exists, if present will link to {self.workdir}")
                if self.file_present(self.positive_control):
                    self.make_links(first_column = "9999-99888", second_column = self.resources / 'control' / 'contigs.fa')
        else:
            self.make_links(first_column = self.prefix, second_column = input_file)     
        
        

    def generate_workflow_files(self):
        # varaiables for config.yaml
        self.logger.info(f"Setting up workflow files. ")
        print(self.species)
        mduqc = "mduqc" if self.mduqc else ""
        if self.mduqc:
            self.file_present(self.qc)
        # if running singleton put summary files in prefix dir
        config_source = self.resources / "templates" / "nextflow.config.j2"
        self.logger.info(f"Writing config file")
        config_template = jinja2.Template(config_source.read_text())
        config_target = self.workdir / "nextflow.config"
        config_target.write_text(
            config_template.render(
                mduqc=mduqc, # true or false if true MUST have provided a qc result table
                qc=self.qc, # name of qc result table
                db = self.db, # version of DB
                outdir = f"{self.workdir}", # the directory where all things should be saved to
                species = self.species if self.species != "" else "none" # if this is not an empty string then run amrfinder plus - for MDU this should not be needed.
            )
        )
        self.logger.info(f"Written nextflow.config to {self.workdir}")

    def run_wflw(self):

        
        cmd = f"nextflow {self.resources / 'main.nf'} -resume 2>&1"
        self.logger.info(f"Running pipeline using command {cmd}. This may take some time.")
        wkfl = subprocess.run(cmd, shell=True)
        while True:
                if wkfl.stdout != None:
                    line = wkfl.stdout.readline().strip()
                    if not line:
                        break
                line = ''
                break
                self.logger.info(f"{line}")
                
        if wkfl.returncode == 0:
            return True
        else:
            return False

    # def clean(self):

    #     logs = self.workdir / ".snakemake" / "log"
    #     if logs.exists():
    #         rmlog = subprocess.run(f"rm -rf {logs}", shell=True, capture_output=True)
    #         if rmlog.returncode == 0:
    #             self.logger.info("Removed old log files")
    #     conda = self.workdir / ".snakemake" / "conda"
    #     if conda.exists():
    #         cleanconda = subprocess.run(
    #             f"snakemake --cleanup-conda", shell=True, capture_output=True
    #         )
    #         if cleanconda.returncode == 0:
    #             self.logger.info("Cleaned unused conda environments")

    def run_amr(self):
        # setup the pipeline
        self.link_input_files()
        # write snakefile
        self.generate_workflow_files()
        # run snakefile
        wkflow = self.run_wflw()
        # if wkflow:
        #     self.logger.info(f"Pipeline completed")
        #     if self.mduqc:
        #         outfile = self.workdir / "MMS118.xlsx"
        #         if outfile.exists():
        #             self.logger.info(f"MMS118.xlsx found, pipeline successfully completed. Come again soon.")
        #         else:
        #             self.logger.warning(f"MMS118.xlsx is not present. Please check logs and try again.")
        #             raise SystemExit
        #     else:
        #         outfile_matches = self.workdir / "summary_matches.csv"
        #         outfile_partials = self.workdir / "summary_partials.csv"
        #         if outfile_matches.exists() and outfile_partials.exists():
        #             self.logger.info(f"'summary_matches.csv' and 'summary_partials.csv' found, pipeline successfully completed. Come again soon.")
        #         else:
        #             self.logger.warning(f"'summary_matches.csv' and 'summary_partials.csv' are not present. Please check logs and try again.")
        #             raise SystemExit
        #     if not self.keep:
        #         self.logger.info(f"Cleaning up the working directory.")
        #         self.clean()
        #     self.logger.info("Thank you and come again!")
        # else:
        #     self.logger.info(f"Pipeline did not complete successfully. Check logs and try again")
