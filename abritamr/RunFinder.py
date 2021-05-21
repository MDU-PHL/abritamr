import pathlib, pandas, datetime, subprocess, os, logging,subprocess,collections, re
from abritamr.version import db
from abritamr.CustomLog import CustomFormatter


class RunFinder(object):
    """
    A class to run amrfinderplus
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
        self.organism = args.organism
        self.input = args.input
        self.run_type = args.run
        self.jobs = args.jobs
        self.prefix = args.prefix
    
    def _batch_cmd(self):
        """
        generate cmd with parallel
        """
        org = f"--plus --organism {self.organism}" if self.organism != '' else ''
        cmd = f"parallel -j {self.jobs} --colsep '\t' mkdir {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out {org} :::: {self.input}"
        return cmd
    
    def _single_cmd(self):
        """
        generate a single amrfinder command
        """
        org = f"--plus --organism {self.organism}" if self.organism != '' else ''
        cmd = f"mkdir {self.prefix} && amrfinder -n {self.input} -o {self.prefix}/amrfinder.out {org}"
        return cmd
    
    def _check_amrfinder(self):
        """
        Check that amrfinder is installed and db setup properly.
        """
        cmd = f"amrfinder --debug"
        pat = re.compile(r'(?P<id>[0-9]{4}-[0-9]{5,6})-?(?P<itemcode>.{1,2})?')
        p = subprocess.run("amrfinder --debug", shell = True, encoding = "utf-8")
        m = re.search(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', p.stderr)
        if m:
            return True
        else:
            return False


    def _generate_cmd(self):
        """
        Generate a command to run amrfinder
        """
        cmd = self._batch_cmd() if self.run_type == 'batch' else self._single_cmd()
        return cmd
        
    def _run_cmd(self, cmd):
        """
        Use subprocess to run the command for amrfinder
        """

        p = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
        if p.returncode == 0:
            self.logger.info(f"AMRfinder completed successfully. Will now move on to collation.")
            return True
        else:
            self.logger.critical(f"There appears to have been a problem with running amrfinder plus. The following erro has been reported : \n {p.stderr}")

    def _check_output_file(self, path):
        """
        check that amrfinderplus outputs are present
        """
        if not pathlib.Path(path).exists():
            self.logger.critical(f"The amrfinder output : {path} is missing. Something has gone wrong with AMRfinder plus. Please check all inputs and try again.")
            raise SystemExit
        else:
            return True

    def _check_outputs(self):
        """
        use inputs to check if files made
        """
        if self.run_type != 'batch':
            self._check_output_file(f"{self.prefix}/amrfinder.out")
        else:
            tab = pandas.read_csv(self.input, sep = '\t', header = None)
            for row in tab.iterrows():
                self._check_output_file(f"{row[1][0]}/amrfinder.out"):

    def run(self):
        """
        run amrfinder
        """
        if self._check_amrfinder():
            cmd = self._generate_cmd()
            self._run_cmd(cmd)
            self._check_outputs()
            
        else:
            self.logger.critical(f"It seems that amrfinder is not properly configured. Please check amrfinder documentation (https://github.com/ncbi/amr/wiki/AMRFinderPlus-database) and try again.")
            raise SystemExit
        Data = collections.namedtuple('Data', ['run_type', 'input', 'prefix'])
        amr_data = Data(self.run_type, self.input, self.prefix)

        return amr_data
