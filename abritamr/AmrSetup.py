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
        self.snakefile = self.resources / "templates" / "Snakefile.smk"
        self.jobs = args.jobs
        self.mduqc = args.mduqc
        self.positive_control = args.positive_control
        # self.drugs = pathlib.Path(args.drug_classes)
        self.contigs = args.contigs
        self.amrfinder_output = args.amrfinder_output
        self.from_contigs = True
        self.prefix = args.prefix
        self.keep = args.keep
        self.run_singulairty = args.Singularity
        self.singularity_path = f"shub://phgenomics-singularity/amrfinderplus" if args.singularity_path == '' else args.singularity_path # this needs to addressed before formal release.
        self.finaloutput = (
            ['MMS118.xlsx']
            if self.mduqc
            else ['summary_matches.csv', 'summary_partials.csv']
        )
        self.qc = args.qc
        self.species_detect = args.species_detect



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
        input_file = (
            pathlib.Path(self.contigs)
            if self.from_contigs
            else pathlib.Path(self.amrfinder_output)
        )
        tab = pandas.read_csv(input_file, engine="python", header=None, sep = '\t')
        running_type = self.check_input_tab(tab)
        # print(running_type)
        if running_type == 'batch':
            isos = list(tab[0])
            self.logger.info(f"Checking that the input data is present. If present will link to {self.workdir}")
            for row in tab.iterrows():
                if self.file_present(row[1][1]):
                    self.make_links(first_column=row[1][0], second_column=row[1][1])
            if self.positive_control != "":
                self.logger.info(f"You have provided a path the a positive control. Checking if path exists, if present will link to {self.workdir}")
                if self.file_present(self.positive_control):
                    self.make_links(first_column = "9999-99888", second_column = self.positive_control)
                    isos.append("9999-99888")
        elif running_type == 'amrfinder_output':
            self.prefix = self.prefix if self.prefix != '' else f'abritamr'
            self.from_contigs = False
            self.make_links(first_column=self.prefix, second_column=input_file)
        elif running_type == 'assembly':
            self.prefix = self.prefix if self.prefix != '' else f'abritamr'
            self.from_contigs = True
            self.make_links(first_column=self.prefix, second_column=input_file)

        return [self.prefix] if running_type != 'batch' else isos

    def check_singularity(self):
        '''
        if singularity path check it exists.
        '''
        self.logger.info(f"Checking singularity.")
        if self.singularity_path != f"shub://phgenomics-singularity/amrfinderplus":
            return self.file_present(self.singularity_path)

    def generate_workflow_files(self, samples):
        # varaiables for config.yaml
        self.logger.info(f"Setting up workflow files. ")
        script_path = self.resources / "utils"
        amrfinder = " " if not self.from_contigs else "run_amrfinder"
        mduqc = "mduqc" if self.mduqc else ""
        species = "species" if self.species_detect else ""
        singularity_path = (
            f"singularity:'{self.singularity_path}'"
            if self.run_singulairty
            else " ")
        if self.mduqc:
            self.file_present(self.qc)
        # if running singleton put summary files in prefix dir
        config_source = self.resources / "templates" / "config.j2"
        self.logger.info(f"Writing config file")
        config_template = jinja2.Template(config_source.read_text())
        config_target = self.workdir / "config.yaml"
        config_target.write_text(
            config_template.render(
                script_path=script_path,
                amrfinder=amrfinder,
                mduqc=mduqc,
                samples=' '.join(samples),
                qc=self.qc,
                final_output=self.finaloutput,
                workdir=self.workdir,
                singularity_path=singularity_path,
                database_version = self.db,
                species_detect = species
            )
        )
        self.logger.info(f"Written config.yaml to {self.workdir}")

    def run_snakemake(self):

        singularity = "--use-singularity --singularity-args '--bind /home'" if self.run_singulairty else ""
        cmd = f"snakemake -s \"{self.snakefile}\" -j {self.jobs} -d {self.workdir} {singularity} 2>&1"
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

    def clean(self):

        logs = self.workdir / ".snakemake" / "log"
        if logs.exists():
            rmlog = subprocess.run(f"rm -rf {logs}", shell=True, capture_output=True)
            if rmlog.returncode == 0:
                self.logger.info("Removed old log files")
        conda = self.workdir / ".snakemake" / "conda"
        if conda.exists():
            cleanconda = subprocess.run(
                f"snakemake --cleanup-conda", shell=True, capture_output=True
            )
            if cleanconda.returncode == 0:
                self.logger.info("Cleaned unused conda environments")

    def run_amr(self):
        # setup the pipeline
        samples = self.link_input_files()
        # write snakefile
        self.generate_workflow_files(samples)
        # run snakefile
        wkflow = self.run_snakemake()
        if wkflow:
            self.logger.info(f"Pipeline completed")
            if self.mduqc:
                outfile = self.workdir / "MMS118.xlsx"
                if outfile.exists():
                    self.logger.info(f"MMS118.xlsx found, pipeline successfully completed. Come again soon.")
                else:
                    self.logger.warning(f"MMS118.xlsx is not present. Please check logs and try again.")
                    raise SystemExit
            else:
                outfile_matches = self.workdir / "summary_matches.csv"
                outfile_partials = self.workdir / "summary_partials.csv"
                if outfile_matches.exists() and outfile_partials.exists():
                    self.logger.info(f"'summary_matches.csv' and 'summary_partials.csv' found, pipeline successfully completed. Come again soon.")
                else:
                    self.logger.warning(f"'summary_matches.csv' and 'summary_partials.csv' are not present. Please check logs and try again.")
                    raise SystemExit
            if not self.keep:
                self.logger.info(f"Cleaning up the working directory.")
                self.clean()
            self.logger.info("Thank you and come again!")
        else:
            self.logger.info(f"Pipeline did not complete successfully. Check logs and try again")
