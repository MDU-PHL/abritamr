import pathlib, pandas, datetime, getpass, jinja2, re, subprocess, os, logging,subprocess
from abritamr.version import db
from abritamr.CustomLog import CustomFormatter


class Setup(object):
    """
    A base class for setting up abritamr return a valid input object for subsequent steps
    """
    def __init__(self, args):
        

        self.logger =logging.getLogger(__name__) 
        self.logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(CustomFormatter())
        fh = logging.FileHandler('abritamr.log')
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
        fh.setFormatter(formatter)
        self.logger.addHandler(ch) 
        self.logger.addHandler(fh)

        
        # self.mduqc = args.mduqc # if this is a MDU collation
        
        
        # self.keep = args.keep
        # mdu
        self.qc = args.qc
        
    def file_present(self, name):
        """
        check file is present
        :name is the path of the file to check
        """
        
        if name == "":
            return False
        elif pathlib.Path(name).exists():
            self.logger.info(f"Checking if file {name} exists")
            return True
        else:
            return False

class SetupAMR(Setup):
    """
    setup amr inputs for amrfinder run
    """
    def __init__(self, args):
        

        # for amr
        self.species_list = ['Acinetobacter_baumannii', "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"]
        self.db = db

        self.logger =logging.getLogger(__name__) 
        self.logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(CustomFormatter())
        fh = logging.FileHandler('abritamr.log')
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
        fh.setFormatter(formatter)
        self.logger.addHandler(ch) 
        self.logger.addHandler(fh)

        
        self.jobs = args.jobs # number of amrfinderplus to run at a time
        
        self.positive_control = True if self.mduqc else args.positive_control 
        self.contigs = args.contigs
        self.from_contigs = True if args.contigs != '' else False
        self.prefix = args.prefix
        # amr
        self.species = args.species if args.species in self.species_list else ""


    def _check_run_type(self, tab):
        """
        Establish the input is a table of ids and contigs or if it is a single path etc - return batch or assembly;
        """
        self.logger.info(f"Checking the structure of your input file.")
        if tab.shape[1] == 2:
            self.logger.info(f"The input file seems to be in the correct format. Thank you.")
            return 'batch'
        elif tab.shape[1] == 1:
            self.logger.info(f"It seems you might be trying to run abriTAMR directly from an assembly file")
            return 'assembly'
        else:
            logging.critical(
                "Your input file should either be a tab delimited file with two columns or the path to contigs. Please check your input and try again."
            )
            raise SystemExit

    def _get_run_type(self):
        """
        Get the type of run return type and object for checking - either a path or a df
        """
        input_file = pandas.read_csv(pathlib.Path(self.contigs))
        running_type = self._check_run_type(input_file)
        return running_type, input_file

    def _input_files(self):
        """
        Ensure that the files (either contigs or amrfinder output) exist and return running type
        """
        
        running_type,tab = self._get_run_type()
        if running_type == 'batch':
            self.logger.info(f"Checking that the input data is present.")
            for row in tab.iterrows():
                if not self.file_present(row[1][1]):
                    self.logging.critical(f"{row[1][1]} is not a valid file path. Please check your input and try again.")
                    raise SystemExit
        elif running_type == 'assembly' and self.file_present(self.contigs):
            self.logging.info(f"{self.contigs} is present. abritamr can proceed.")
        else:
            self.logging.critical(f"Something has gone wrong with your inputs. Please try again.")
        
        return running_type
            
              

        
        

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
