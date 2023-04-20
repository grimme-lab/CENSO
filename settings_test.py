import unittest
import os
from functools import reduce

from censo.inputhandling import cml
from censo.cfg import DESCR
from censo.settings import CensoSettings


def parentdir(path):
    tmp = os.path.split(path)
    tmp = tmp[:len(tmp)-1]
    return reduce(lambda x, y: os.path.join(x, y), tmp)


test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz -solvent water".split())
test_dir = os.getcwd()
test_path = os.path.join(parentdir(__file__), "censorc_new")

class SettingsTest(unittest.TestCase):
    def setUp(self):
        self.test = CensoSettings()
        
        
    def test_write_config(self):
        self.test.censorc_path = self.test.write_config(test_args, test_dir)
        self.assertEqual(self.test.censorc_path, test_path)
        os.remove(test_path)
        
        
    def test_settings_current(self):
        self.test.censorc_path = self.test.write_config(test_args, test_dir)
        self.test.censorc_path = test_path
        self.test.settings_current = test_args
        os.remove(test_path)
        
        # setting for solvent should be equal to solvent given in cml args
        self.assertEqual(
            self.test.settings_current.get_setting(
                type_t=str, 
                part="general", 
                name="solvent"
            ).value, 
            "water"
        )
        
        # should not find a setting for this combination of parameters
        self.assertIsNone(self.test.settings_current.get_setting(type_t=bool, part="general", name="solvent"))
        
        
    def test_read_program_paths(self):
        self.test.censorc_path = self.test.write_config(test_args, test_dir)
        self.test.censorc_path = test_path
        self.test.settings_current = test_args
        
        self.test.read_program_paths()
        os.remove(test_path)
        print(self.test.external_paths)