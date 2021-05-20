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

    def _check_prefix(self):
        """
        If running type is not batch, then check that prefix is present and a string
        """

        if self.prefix == '':
            self.logger.critical(f"You must supply a sample or sequence id.")
            raise SystemExit
        else:
            return True

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
   

    def setup(self):
        # check that inputs are correct and files are present
        running_type = self._input_files()
        # check that prefix is present (if needed)
        if running_type == 'assembly':
            self._check_prefix()
        Data = collections.namedtuple('Data', ['run_type', 'input', 'prefix', 'jobs', 'organism'])
        input_data = Data(running_type, self.contigs, self.prefix, self.jobs, self.species)
        
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
        self.matches = args.matches
        self.partials = args.matches
        
    def setup(self):
        """
        Check the inputs for MDU - ensure all files are present for collation.
        """
        Data = collections.namedtuple('Data', ['qc', 'matches', 'partials', 'db'])
        if self.file_present(self.qc) and self.file_present(self.matches) and self.file_present(self.partials):
            return Data(self.qc, self.matches, self.partials, self.db)
        else:
            self.logger.critical(f"Something has gone wrong with your inputs. Please try again!")
            raise SystemExit
