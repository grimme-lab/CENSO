import unittest
import os
import numpy as np
from copy import deepcopy


from censo_test.inputhandling import cml
from censo_test.cfg import DESCR
from censo_test.settings import InternalSettings

test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz -solvent water".split())
test_dir = os.getcwd()
test_path = "/home/polt/prog/CENSO/debug/censorc_new"

class SettingsTest(unittest.TestCase):
    def setUp(self):
        self.test = InternalSettings()
        
        
    def test_write_config(self):
        self.test.censorc_path = self.test.write_config(test_args, test_dir)
        self.assertEqual(self.test.censorc_path, test_path)
        
        
    def test_settings_current(self):
        self.test.censorc_path = test_path
        self.test.settings_current = test_args
        
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