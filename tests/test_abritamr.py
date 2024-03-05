import sys, pathlib, pandas, pytest, numpy, logging, collections

from unittest.mock import patch, PropertyMock

from abritamr.AmrSetup import Setup, SetupAMR, SetupMDU
from abritamr.RunFinder import RunFinder
from abritamr.Collate import Collate, MduCollate



test_folder = pathlib.Path(__file__).parent
REFGENES = f"{pathlib.Path(__file__).parent.parent /'abritamr' /'db' / 'refgenes_latest.csv'}"
CONTROLS = pathlib.Path(__file__).parent.parent /'abritamr' / 'control'

def test_file_present():
    """
    assert true when the input file is true
    """
    with patch.object(Setup, "__init__", lambda x: None):
        amr_obj = Setup()
        amr_obj.logger = logging.getLogger(__name__)
        p = CONTROLS / "contigs.fa"
    
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
        amr_obj.contigs = f"{CONTROLS / 'contigs.fa'}"
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
        amr_obj.contigs = f"{CONTROLS / 'contigs.fa'}"
        amr_obj.prefix = 'somename'
        amr_obj.jobs  = 16
        amr_obj.species = ''
        amr_obj.identity = ''
        amr_obj.amrfinder_db = f"{pathlib.Path(__file__).parent.parent /'abritamr' /'db' / 'amrfinderplus'/ 'data'/ '2022-08-09.1'}"
        amr_obj.logger = logging.getLogger(__name__)
        T = collections.namedtuple('T', ['run_type', 'input', 'prefix', 'jobs', 'organism', 'identity','amrfinder_db'])
        input_data = T('assembly', amr_obj.contigs, amr_obj.prefix, amr_obj.jobs, amr_obj.species, amr_obj.identity, amr_obj.amrfinder_db)
        assert amr_obj.setup() == input_data

def test_species():
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{CONTROLS / 'contigs.fa'}"
        amr_obj.prefix = 'somename'
        amr_obj.jobs  = 16
        amr_obj.species = 'Neiserria'
        amr_obj.identity = ''
        amr_obj.amrfinder_db = f"{pathlib.Path(__file__).parent.parent /'abritamr' /'db' / 'amrfinderplus'/ 'data'/ '2022-08-09.1'}"
        amr_obj.logger = logging.getLogger(__name__)
        T = collections.namedtuple('T', ['run_type', 'input', 'prefix', 'jobs', 'organism', 'identity','amrfinder_db'])
        input_data = T('assembly', amr_obj.contigs, amr_obj.prefix, amr_obj.jobs, amr_obj.species, amr_obj.identity, amr_obj.amrfinder_db)
        assert amr_obj.setup() == input_data


def test_setup_batch():
    with patch.object(SetupAMR, "__init__", lambda x: None):
        amr_obj = SetupAMR()
        amr_obj.contigs = f"{test_folder / 'batch.txt'}"
        amr_obj.prefix = ''
        amr_obj.jobs  = 16
        amr_obj.species = ''
        amr_obj.identity = ''
        amr_obj.amrfinder_db = f"{pathlib.Path(__file__).parent.parent /'abritamr' /'db' / 'amrfinderplus'/ 'data'/ '2022-08-09.1'}"
        amr_obj.logger = logging.getLogger(__name__)
        T = collections.namedtuple('T', ['run_type', 'input', 'prefix', 'jobs', 'organism','identity', 'amrfinder_db'])
        input_data = T('batch', amr_obj.contigs, amr_obj.prefix, amr_obj.jobs, amr_obj.species, amr_obj.identity,amr_obj.amrfinder_db)
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

# # Test SetupMDU
MDU = collections.namedtuple('MDU', ['runid', 'matches', 'partials', 'qc', 'sop', 'sop_name'])

def test_prefix_string():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("RUNID", 'tests/summary_matches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv', 'general', 'sop_name')
        amr_obj = SetupMDU(args)
        # amr_obj.runid= 
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_runid()

def test_prefix_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("", 'tests/summary_matches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv', 'general', 'sop_name')
        amr_obj = SetupMDU(args)
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj._check_runid()


def test_mdu_setup_success():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("RUNID", 'tests/summary_matches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv', 'general', 'sop_name')
        amr_obj = SetupMDU(args)
        Data = collections.namedtuple('Data', ['qc', 'matches', 'partials', 'db', 'runid','sop', 'sop_name'])
        d = Data(args.qc, args.matches, args.partials, amr_obj.db, args.runid, args.sop, args.sop_name)
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.setup() == d


def test_mdu_setup_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(SetupAMR, "__init__", lambda x: None):
        args = MDU("RUNID", 'tests/summarymatches.txt', 'tests/summary_matches.txt', 'tests/mdu_qc_checked.csv', 'general', 'sop_name')
        amr_obj = SetupMDU(args)
        Data = collections.namedtuple('Data', ['qc', 'matches', 'partials', 'db', 'runid'])
        d = Data(args.qc, args.matches, args.partials, amr_obj.db, args.runid)
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj.setup()

# # test RunFinder

Data = collections.namedtuple('Data', ['run_type', 'input', 'prefix', 'jobs', 'organism'])
def test_batch_cmd_no_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        
        args = Data("batch", 'tests/batch.txt', '', 9, '')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"parallel -j {args.jobs} --colsep '\\t' 'mkdir -p {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out --plus  --threads 1 -d {amr_obj.amrfinder_db}' :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._batch_cmd() == cmd

def test_batch_cmd_with_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', '', 9, 'Salmonella')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"parallel -j {args.jobs} --colsep '\\t' 'mkdir -p {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out --plus --organism {args.organism} --threads 1 -d {amr_obj.amrfinder_db}' :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._batch_cmd() == cmd

def test_batch_cmd_with_neiserria():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', '', 9, 'Neiserria')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"parallel -j {args.jobs} --colsep '\\t' 'mkdir -p {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out --plus --organism {args.organism} --threads 1 -d {amr_obj.amrfinder_db}' :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._batch_cmd() == cmd

def test_single_cmd_with_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("batch", 'tests/contigs.fa', 'somename', 9, 'Salmonella')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"mkdir -p {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out --plus --organism {args.organism} --threads {args.jobs} -d {amr_obj.amrfinder_db}"
        amr_obj.logger = logging.getLogger()
        assert amr_obj._single_cmd() == cmd


def test_single_cmd_no_org():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("batch", 'tests/contigs.fa', 'somename', 9, '')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"mkdir -p {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out --plus  --threads {args.jobs} -d {amr_obj.amrfinder_db}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._single_cmd() == cmd



def test_generate_cmd_single_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("batch", 'tests/contigs.fa', 'somename', 9, '')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"mkdir -p {args.prefix} && amrfinder -n {args.input} -o {args.prefix}/amrfinder.out -d {amr_obj.amrfinder_db}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._generate_cmd() != cmd


def test_generate_cmd_batch_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("assembly", 'tests/batch.txt', '', 9, '')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.amrfinder_db = "2021-06-01.1"
        amr_obj.identity = ''
        cmd = f"parallel -j {args.jobs} --colsep '\t' mkdir -p {{1}} && amrfinder -n {{2}} -o {{1}}/amrfinder.out -d {amr_obj.amrfinder_db} :::: {args.input}"
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._generate_cmd() != cmd


def test_check_output_file():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'tests', 9, '')
        amr_obj = RunFinder()
        p = 'tests/amrfinder.out'
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_output_file(p)


def test_check_output_file_fail():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'tests', 9, '')
        amr_obj = RunFinder()
        p = 'tests/amrfinders.out'
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            amr_obj._check_output_file(p)


def test_check_outputs_single():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("assembly", 'tests/contigs.fa', 'tests', 9, '')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_outputs()



def test_check_outputs_batch():
    """
    assert True when non-empty string is given
    """
    with patch.object(RunFinder, "__init__", lambda x: None):
        args = Data("batch", 'tests/batch.txt', 'tests', 9, '')
        amr_obj = RunFinder()
        amr_obj.organism = args.organism
        amr_obj.run_type = args.run_type
        amr_obj.prefix = args.prefix
        amr_obj.jobs = args.jobs
        amr_obj.input = args.input
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj._check_outputs()

# # test Collate
# # 
Colls = collections.namedtuple('Data', ['run_type', 'input', 'prefix'])
def test_get_drugclass_allele_1():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        reftab = pandas.read_csv(REFGENES)
        df= pandas.read_csv('tests/amrfinder.out', sep = '\t')
        for i in df.iterrows():
            if i[1]['Gene symbol'] == 'blaSHV-11':
                row = i
        colname = 'allele'
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.get_drugclass(reftab, row, colname) == "Beta-lactam"

def test_extract_gene_name_1():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        reftab = pandas.read_csv(REFGENES)
        protein = 'WP_004176269.1'
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.extract_gene_name(protein, reftab) == 'blaSHV-11'

def test_extract_gene_name_2():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        reftab = pandas.read_csv(REFGENES)
        reftab = reftab.fillna('-')
        protein = 'WP_063839881.1'
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.extract_gene_name(protein, reftab) == "aac(2')-IIa"

def test_setup_dict():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        reftab = pandas.read_csv(REFGENES)
        reftab = reftab.fillna('-')
        drugclass_dict = {}
        df= pandas.read_csv('tests/amrfinder.out', sep = '\t')
        for i in df.iterrows():
            if i[1]['Gene symbol'] == 'blaSHV-11':
                row = i
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.setup_dict(drugclass_dict, reftab, row) == {"Beta-lactam":['blaSHV-11']}


def test_get_per_isolate():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        reftab = pandas.read_csv(REFGENES)
        reftab = reftab.fillna('-')
        
        df= pandas.read_csv('tests/amrfinder.out', sep = '\t')
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.get_per_isolate(reftab=reftab, df=df, isolate=isolate) == ({"Isolate":isolate, "Beta-lactam":'blaSHV-11','ESBL':'blaCTX-M-15'},{"Isolate":isolate,'ESBL':'blaCTX-M-15'},{"Isolate":isolate,'Metal':'qnrB1'})



def test_collate():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        reftab = pandas.read_csv(REFGENES)
        reftab = reftab.fillna('-')
        isolate = 'tests'
        drugs = pandas.DataFrame({"Isolate":isolate,'ESBL': 'blaCTX-M-15', "Beta-lactam":'blaSHV-11'}, index = [0])
        partial = pandas.DataFrame({"Isolate":isolate,'ESBL':'blaCTX-M-15'}, index = [0])
        virulence = pandas.DataFrame({"Isolate":isolate,'Metal':'qnrB1'}, index = [0])
        # print(drugs)
        
        amr_obj.logger = logging.getLogger(__name__)
        assert amr_obj.collate(isolate)[0].equals(drugs)
        assert amr_obj.collate(isolate)[1].equals(partial)
        assert amr_obj.collate(isolate)[2].equals(virulence)


def test_save():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate,'ESBL': 'blaCTX-M-15', "Beta-lactamase (not ESBL or carbapenemase)":'blaSHV-11'}, index = [0])
        summary_partial = pandas.DataFrame({"Isolate":isolate,'ESBL':'blaCTX-M-15'}, index = [0])
        virulence = pandas.DataFrame({"Isolate":isolate,'METAL':'qnrB1'}, index = [0])
        assert amr_obj.save_files('',summary_drugs,summary_partial, virulence)



def test_save_one_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate}, index = [0])
        summary_partial = pandas.DataFrame({"Isolate":isolate,'ESBL':'blaCTX-M-15'}, index = [0])
        virulence = pandas.DataFrame({"Isolate":isolate,'METAL':'qnrB1'}, index = [0])
        assert amr_obj.save_files('',summary_drugs,summary_partial, virulence)



def test_save_all_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate}, index = [0])
        summary_partial = pandas.DataFrame({"Isolate":isolate}, index = [0])
        virulence = pandas.DataFrame({"Isolate":isolate}, index = [0])
        assert amr_obj.save_files('',summary_drugs,summary_partial, virulence)


def test_get_cols_works():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate,'ESBL': 'blaCTX-M-15', "Beta-lactamase (not ESBL or carbapenemase)":'blaSHV-11'}, index = [0])
        # summary_partial = pandas.DataFrame({"Isolate":isolate,'ESBL':'blaCTX-M-15'}, index = [0])
        # virulence = pandas.DataFrame({"Isolate":isolate,'METAL':'qnrB1'}, index = [0])
        assert amr_obj._get_cols(df= summary_drugs) == ["ESBL","Beta-lactamase (not ESBL or carbapenemase)"]

def test_get_cols_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame()
        # summary_partial = pandas.DataFrame({"Isolate":isolate,'ESBL':'blaCTX-M-15'}, index = [0])
        # virulence = pandas.DataFrame({"Isolate":isolate,'METAL':'qnrB1'}, index = [0])
        assert amr_obj._get_cols(df= summary_drugs) == []


def test_commbine_all_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame()
        summary_partial = pandas.DataFrame()
        virulence = pandas.DataFrame()
        assert amr_obj._combine_dfs(summary_drugs,summary_partial, virulence).empty



def test_save_merge():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate,'ESBL': 'blaCTX-M-15', "Beta-lactamase (not ESBL or carbapenemase)":'blaSHV-11'}, index = [0])
        summary_partial = pandas.DataFrame({"Isolate":isolate,'ESBL':'blaCTX-M-15'}, index = [0])
        virulence = pandas.DataFrame()
        assert not amr_obj._merge(summary_drugs,summary_partial).empty



def test_save_merge_one_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate,'ESBL': 'blaCTX-M-15', "Beta-lactamase (not ESBL or carbapenemase)":'blaSHV-11'}, index = [0])
        summary_partial = pandas.DataFrame({"Isolate":isolate}, index = [0])
        virulence = pandas.DataFrame()
        assert not amr_obj._merge(summary_drugs,summary_partial).empty



def test_save_merge_both_empty():
    """
    assert True when non-empty string is given
    """
    with patch.object(Collate, "__init__", lambda x: None):
        args = Colls("assembly", 'tests/contigs.fa', '')
        amr_obj = Collate()
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate}, index = [0])
        summary_partial = pandas.DataFrame({"Isolate":isolate}, index = [0])
        virulence = pandas.DataFrame()
        assert not amr_obj._merge(summary_drugs,summary_partial).empty



MduColls = collections.namedtuple('Data', ['sop', 'sop_name', 'qc','db','partials','matches','runid'])
def test_reporting_efaec_optra():

    with patch.object(Collate, "__init__", lambda x: None):
        args = MduColls('sop', 'sop_name', 'mduqc','db','partials','match','runid')
        amr_obj = MduCollate(args)
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate,'Florfenicol/Oxazolidinone': 'optrA'}, index = [0])
        
        for row in summary_drugs.iterrows():
            assert amr_obj.reporting_logic_general(row = row, species = 'Enterococcus faecium') == (['optrA'],["No non-reportable genes found."])


def test_reporting_ni_optra():

    with patch.object(Collate, "__init__", lambda x: None):
        args = MduColls('sop', 'sop_name', 'mduqc','db','partials','match','runid')
        amr_obj = MduCollate(args)
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        summary_drugs = pandas.DataFrame({"Isolate":isolate,'Florfenicol/Oxazolidinone': 'optrA'}, index = [0])
        
        for row in summary_drugs.iterrows():
            assert amr_obj.reporting_logic_general(row = row, species = 'No identification') == (['CPase_16S_mcr_NEG'],['optrA'])



def test_reporting_df_efaec():

    with patch.object(Collate, "__init__", lambda x: None):
        args = MduColls('sop', 'sop_name', f"{test_folder / 'mdu_qc_checked.csv'}",'db',f"{test_folder / 'summary_partials.txt'}",f"{test_folder / 'summary_matches.txt'}",'runid')
        amr_obj = MduCollate(args)
        isolate = 'tests'
        amr_obj.logger = logging.getLogger(__name__)
        
        df = pandas.DataFrame(
            {"MDU sample ID":isolate,
            'Item code':"",
            'Resistance genes (alleles) detected':"optrA",
            'Resistance genes (alleles) det (non-rpt)':"No non-reportable genes found.",
            'Species_obs':'Enterococcus faecium',
            'Species_exp':'No identification', 
            'db_version':'db'},
            index = [0]
        )
        df = df.set_index("MDU sample ID")
        
        assert amr_obj.mdu_reporting_general(match=f"{test_folder / 'summary_matches.txt'}").equals(df)



def test_reporting_df_shig():

    with patch.object(Collate, "__init__", lambda x: None):
        args = MduColls('sop', 'sop_name', f"{test_folder / 'mdu_qc_checked.csv'}",'db',f"{test_folder / 'summary_partials_shig.txt'}",f"{test_folder / 'summary_matches_shig.txt'}",'runid')
        amr_obj = MduCollate(args)
        isolate = 'tests2'
        amr_obj.logger = logging.getLogger(__name__)
        
        df = pandas.DataFrame(
            {"MDU sample ID":isolate,
            'Item code':"",
            'Resistance genes (alleles) detected':"CTX-M-15",
            'Resistance genes (alleles) det (non-rpt)':"No non-reportable genes found.",
            'Species_obs':'Eschericia coli', 
            'Species_exp':'Shigella flexneri',
            'db_version':'db'},
            index = [0]
        )
        df = df.set_index("MDU sample ID")
        
        assert amr_obj.mdu_reporting_general(match=f"{test_folder / 'summary_matches_shig.txt'}").equals(df)

