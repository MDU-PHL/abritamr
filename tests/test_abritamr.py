import sys, pathlib, pandas, pytest, numpy, logging, collections

from unittest.mock import patch, PropertyMock

from abritamr.AmrSetup import Setup, SetupAMR, SetupMDU
from abritamr.RunFinder import RunFinder
from abritamr.Collate import Collate, MduCollate

test_folder = pathlib.Path(__file__).parent


def test_file_present():
    """
    assert true when the input file is true
    """
    with patch.object(Setup, "__init__", lambda x: None):
        amr_obj = Setup()
        amr_obj.logger = logging.getLogger(__name__)
        p = test_folder / "contigs.fa"
    
        assert amr_obj.file_present(p)

# Test SetupAMR
def test_check_run_type_batch():
    """
    assert true when a tab file is provided with two columns
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'batch.txt'}"
        amr_obj.prefix = ''
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._get_input_shape() == 'batch'

def test_check_run_type_single():
    """
    assert true when a contigs file is provided and batch recorded as assembly
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'contigs.fa'}"
        amr_obj.prefix = 'somename'
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._get_input_shape() == 'assembly'

def test_check_run_type_wrong_data():
    """
    assert true when error raised because input is wrong shape
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'batch_fail.txt'}"
        amr_obj.prefix = ''
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj._get_input_shape()

def test_prefix_string():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.prefix = "some_isolate"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_prefix()

def test_prefix_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.prefix = ""
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj._check_prefix()
       

def test_setup_contigs():
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'contigs.fa'}"
        amr_obj.prefix = 'somename'
        amr_obj.jobs  = 16
        amr_obj.species = ''
        amr_obj.logger = logging.getLogger(__name__)
        T = collections.namedtuple('T', ['run_type', 'input', 'prefix', 'jobs', 'organism'])
        input_data = T('assembly', amr_obj.contigs, amr_obj.prefix, amr_obj.jobs, amr_obj.species)
        assert amr_obj.setup() == input_data


def test_setup_fail():
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'batch.txt'}"
        amr_obj.prefix = ''
        amr_obj.jobs  = 16
        amr_obj.species = ''
        amr_obj.logger = logging.getLogger(__name__)
        T = collections.namedtuple('T', ['run_type', 'input', 'prefix', 'jobs', 'organism'])
        input_data = T('batch', amr_obj.contigs, amr_obj.prefix, amr_obj.jobs, amr_obj.species)
        assert amr_obj.setup() == input_data
 
def test_setup_fail():
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'batch_fail.txt'}"
        amr_obj.prefix = ''
        amr_obj.jobs  = 16
        amr_obj.species = ''
        amr_obj.logger = logging.getLogger(__name__)
        T = collections.namedtuple('T', ['run_type', 'input', 'prefix', 'jobs', 'organism'])
        input_data = T('batch', amr_obj.contigs, amr_obj.prefix, amr_obj.jobs, amr_obj.species)
        with pytest.raises(SystemExit):
            amr_obj.setup()

# Test SetupMDU
MDU = collections.namedtuple('MDU', ['runid', 'matches', 'partial', 'qc'])
def test_prefix_string():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("RUNID", 'tests/summary_matches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv')
        amr_obj = SetupMDU(args)
        # amr_obj.runid= 
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_runid()

def test_prefix_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("", 'tests/summary_matches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv')
        amr_obj = SetupMDU(args)
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj._check_runid()


def test_mdu_setup_success():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("RUNID", 'tests/summary_matches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv')
        amr_obj = SetupMDU(args)
        Data = collections.namedtuple('Data', ['qc', 'matches', 'partials', 'db', 'runid'])
        d = Data(args.qc, args.matches, args.partial, amr_obj.db, args.runid)
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.setup() == d


def test_mdu_setup_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("RUNID", 'tests/summarymatches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv')
        amr_obj = SetupMDU(args)
        Data = collections.namedtuple('Data', ['qc', 'matches', 'partials', 'db', 'runid'])
        d = Data(args.qc, args.matches, args.partial, amr_obj.db, args.runid)
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj.setup()

# test RunFInder

Data = collections.namedtuple('Data', ['run_type', 'input', 'prefix', 'jobs', 'organism'])
def test_batch_cmd_no_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', '', 9, '')
        amr_obj = RunFinder(args)
        cmd = f"parallel -j {args.jobs} --colsep '\t' mkdir {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out  :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._batch_cmd() == cmd

def test_batch_cmd_with_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', '', 9, 'Salmonella')
        amr_obj = RunFinder(args)
        cmd = f"parallel -j {args.jobs} --colsep '\t' mkdir {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out --plus --organism {args.organism} :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._batch_cmd() == cmd

def test_single_cmd_with_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/contigs.fa', 'somename', 9, 'Salmonella')
        amr_obj = RunFinder(args)
        cmd = f"mkdir {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out --plus --organism {args.organism}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._single_cmd() == cmd


def test_single_cmd_no_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/contigs.fa', 'somename', 9, '')
        amr_obj = RunFinder(args)
        cmd = f"mkdir {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out "
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._single_cmd() == cmd


def test_generate_cmd_single():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'somename', 9, '')
        amr_obj = RunFinder(args)
        cmd = f"mkdir {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out "
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._generate_cmd() == cmd


def test_generate_cmd_single_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/contigs.fa', 'somename', 9, '')
        amr_obj = RunFinder(args)
        cmd = f"mkdir {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out "
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._generate_cmd() != cmd

def test_generate_cmd_batch():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', '', 9, '')
        amr_obj = RunFinder(args)
        cmd = f"parallel -j {args.jobs} --colsep '\t' mkdir {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out  :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._generate_cmd() == cmd


def test_generate_cmd_batch_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("assembly", 'tests/batch.txt', '', 9, '')
        amr_obj = RunFinder(args)
        cmd = f"parallel -j {args.jobs} --colsep '\t' mkdir {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out  :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._generate_cmd() != cmd


def test_check_output_file():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'tests', 9, '')
        amr_obj = RunFinder(args)
        p = 'tests/amrfinder.out'
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_output_file(p)


def test_check_output_file_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'tests', 9, '')
        amr_obj = RunFinder(args)
        p = 'tests/amrfinders.out'
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj._check_output_file(p)


def test_check_outputs_single():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'tests', 9, '')
        amr_obj = RunFinder(args)
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_outputs()



def test_check_outputs_batch():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', 'tests', 9, '')
        amr_obj = RunFinder(args)
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_outputs()

# # def test_input_exists_amrfinder():

# #     with patch.object(Setupamr, "__init__", lambda x: None):
# #         amr_obj = Setupamr()
# #         amr_obj.contigs = ""
# #         amr_obj.mduqc = False
# #         amr_obj.amrfinder_output = test_folder / "amrfinder_input_file.txt"
# #         amr_obj.logger = logging.getLogger(__name__)
# #         assert amr_obj.check_input_exists() == True


# # def test_input_exists_amrfinder_contigs_fail():

# #     with patch.object(Setupamr, "__init__", lambda x: None):
# #         amr_obj = Setupamr()
# #         amr_obj.contigs = test_folder / "contigs_input_file.txt"
# #         amr_obj.mduqc = False
# #         amr_obj.amrfinder_output = test_folder / "amrfinder_input_file.txt"
# #         amr_obj.logger = logging.getLogger(__name__)
# #         with pytest.raises(SystemExit):
# #             amr_obj.check_input_exists()


# # def test_input_exists_fail():

# #     with patch.object(Setupamr, "__init__", lambda x: None):
# #         amr_obj = Setupamr()
# #         amr_obj.contigs = f""
# #         amr_obj.mduqc = False
# #         amr_obj.amrfinder_output = f""
# #         amr_obj.logger = logging.getLogger(__name__)
# #         with pytest.raises(SystemExit):
# #             amr_obj.check_input_exists()


# # def test_input_exists_mduqc_fail():

# #     with patch.object(Setupamr, "__init__", lambda x: None):
# #         amr_obj = Setupamr()
# #         amr_obj.contigs = ""
# #         amr_obj.mduqc = True
# #         amr_obj.amrfinder_output = test_folder / "amrfinder_input_file.txt"
# #         amr_obj.logger = logging.getLogger(__name__)
# #         with pytest.raises(SystemExit):
# #             amr_obj.check_input_exists()
