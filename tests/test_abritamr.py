import sys, pathlib, pandas, pytest, numpy

from unittest.mock import patch

from abritamr.AmrSetup import Setupamr


def test_file_present():
        '''
        assert true when the input file is true
        '''
        with patch.object(Setupamr, "__init__", lambda x: None):
                amr_obj = Setupamr()
                p = pathlib.Path('abritamr','test', 'contigs_input_file.txt')
                assert amr_obj.file_present(p)


def test_structure():
        with patch.object(Setupamr, "__init__", lambda x: None):
                amr_obj = Setupamr()
                tab = pandas.DataFrame({'A':[1,2,3,4], 'B':[5,6,7,8]})
                assert amr_obj.check_input_tab(tab) == True

def test_structure_fail():
        with patch.object(Setupamr, "__init__", lambda x: None):
                amr_obj = Setupamr()
                tab = pandas.DataFrame({'A':[1,2,3,4], 'B':[5,6,7,8], 'C':[1,2,5,4]})
                with pytest.raises(SystemExit):
                        amr_obj.check_input_tab(tab)


def test_input_exists_conitgs():

        with patch.object(Setupamr, "__init__", lambda x:None):
                amr_obj = Setupamr()     
                amr_obj.contigs = f"{pathlib.Path('abritamr', 'test', 'contigs_input_file.txt')}"
                print(amr_obj.contigs)
                amr_obj.mduqc = False
                amr_obj.amrfinder_output = ''
                assert amr_obj.check_input_exists() == True

def test_input_exists_amrfinder():

        with patch.object(Setupamr, "__init__", lambda x:None):
                amr_obj = Setupamr()     
                amr_obj.contigs = ''
                amr_obj.mduqc = False
                amr_obj.amrfinder_output = f"{pathlib.Path('abritamr', 'test', 'amrfinder_input_file.txt')}"
                assert amr_obj.check_input_exists() == True

def test_input_exists_amrfinder_contigs_fail():

        with patch.object(Setupamr, "__init__", lambda x:None):
                amr_obj = Setupamr()     
                amr_obj.contigs = f"{pathlib.Path('abritamr', 'test', 'contigs_input_file.txt')}"
                amr_obj.mduqc = False
                amr_obj.amrfinder_output = f"{pathlib.Path('abritamr', 'test', 'amrfinder_input_file.txt')}"
                with pytest.raises(SystemExit):
                        amr_obj.check_input_exists()

def test_input_exists_fail():

        with patch.object(Setupamr, "__init__", lambda x:None):
                amr_obj = Setupamr()     
                amr_obj.contigs = f""
                amr_obj.mduqc = False
                amr_obj.amrfinder_output = f""
                with pytest.raises(SystemExit):
                        amr_obj.check_input_exists()

def test_input_exists_mduqc_fail():

        with patch.object(Setupamr, "__init__", lambda x:None):
                amr_obj = Setupamr()     
                amr_obj.contigs = ''
                amr_obj.mduqc = True
                amr_obj.amrfinder_output = f"{pathlib.Path('abritamr', 'test', 'amrfinder_input_file.txt')}"
                with pytest.raises(SystemExit):
                        amr_obj.check_input_exists()
