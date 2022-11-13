import pathlib, pandas, datetime, subprocess, os, logging,subprocess,collections
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
        self.species_list = ["Burkholderia_cepacia","Acinetobacter_baumannii","Streptococcus_pyogenes","Streptococcus_agalactiae","Streptococcus_pneumoniae","Enterococcus_faecium","Pseudomonas_aeruginosa","Staphylococcus_pseudintermedius","Clostridioides_difficile","Klebsiella","Neisseria","Campylobacter","Salmonella","Escherichia","Staphylococcus_aureus","Burkholderia_pseudomallei","Enterococcus_faecalis"]
        
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

        
        self.jobs = args.jobs # number of amrfinderplus to run at a time
        self.contigs = args.contigs
        self.from_contigs = True if args.contigs != '' else False
        self.prefix = args.prefix
        # amr
        self.species = args.species if args.species in self.species_list else ""
        self.identity = args.identity
        self.amrfinder_db = args.amrfinder_db

        

    def _check_prefix(self):
        """
        If running type is not batch, then check that prefix is present and a string
        """

        if self.prefix == '':
            self.logger.critical(f"You must supply a sample or sequence id.")
            raise SystemExit
        else:
            return True

            
    
    def _get_input_shape(self):
        """
        determine shape of file
        """
        run_type = 'assembly'
        with open(self.contigs, 'r') as c:
            data = c.read().strip().split('\n')
            firstline = data[0]
            if not firstline.startswith('>'):
                for line in data:
                    if len(line.split('\t')) != 2:
                        self.logger.critical("Your input file should either be a tab delimited file with two columns or the path to contigs. Please check your input and try again.")
                        raise SystemExit
                run_type = 'batch'
        self.logger.info(f"The input file seems to be in the correct format. Thank you.")
        return run_type
    

    def _input_files(self):
        """
        Ensure that the files (either contigs or amrfinder output) exist and return running type
        """
        
        running_type = self._get_input_shape()
        if running_type == 'batch':
            self.logger.info(f"Checking that the input data is present.")
            with open(self.contigs, 'r') as c:
                data = c.read().strip().split('\n')
                for line in data:
                    row = line.split('\t')
                    if not self.file_present(row[1]):
                        self.logger.critical(f"{row[1]} is not a valid file path. Please check your input and try again.")
                        raise SystemExit
        elif running_type == 'assembly' and self.file_present(self.contigs):
            self.logger.info(f"{self.contigs} is present. abritamr can proceed.")
        else:
            self.logger.critical(f"Something has gone wrong with your inputs. Please try again.")
            raise SystemExit
        
        return running_type
   

    def setup(self):
        # check that inputs are correct and files are present
        running_type = self._input_files()
        # check that prefix is present (if needed)
        if running_type == 'assembly':
            self._check_prefix()
        
        Data = collections.namedtuple('Data', ['run_type', 'input', 'prefix', 'jobs', 'organism', 'identity','amrfinder_db'])
        input_data = Data(running_type, self.contigs, self.prefix, self.jobs, self.species, self.identity, self.amrfinder_db)
        
        return input_data


class SetupMDU(Setup):
    """
    Setup MDUify of abritamr results
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
        self.db = db
        self.qc = args.qc
        self.runid = args.runid
        self.matches = args.matches
        self.partials = args.partials   
        self.sop = args.sop
        self.sop_name = args.sop_name

    def _check_runid(self):
        if self.runid == '':
            self.logger.critical(f"Run ID can not be empty, please try again.")
            raise SystemExit
        else:
            return True

    def setup(self):
        """
        Check the inputs for MDU - ensure all files are present for collation.
        """
        
        file_dict = {
            'QC': self.qc, 
            'summary_matches': self.matches,
            'summary_partials':self.partials
            } if self.sop == 'general' else {
            'QC': self.qc, 
            'summary_matches': self.matches,
            # 'summary_partials':self.partials
            }

        if self._check_runid():
            self.logger.info(f"You are generating a {'general report' if self.sop == 'general' else 'species specific report'}")
            self.logger.info(f"Now checking all input files are present.")
            for _file in file_dict:
                self.logger.info(f"Checking that {_file} is present.")
                if self.file_present(file_dict[_file]):
                    self.logger.info(f"{file_dict[_file]} is present.")
                else:
                    self.logger.critical(f"The {_file} file supplied ({file_dict[_file]}) does not exist. Please check your inputs and try again.")
                    raise SystemExit

            Data = collections.namedtuple('Data', ['qc', 'matches', 'partials', 'db', 'runid', 'sop','sop_name'])
        
            return Data(self.qc, self.matches, self.partials, self.db, self.runid, self.sop, self.sop_name)
        