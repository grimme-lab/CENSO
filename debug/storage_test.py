import unittest
import os

from censo_test.inputhandling import cml
from censo_test.cfg import DESCR
from censo_test.storage import CensoStorage

test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz".split())
test_dir = os.getcwd()

class StorageTest(unittest.TestCase):
    def setUp(self):
        self.test = CensoStorage(test_args, test_dir)
        
        
    def test_assets_path(self):
        self.assertEqual(
            self.test.assets_path, 
            os.path.expanduser("~/.censo_assets")
        )
        
        
    def test_ensemble_path(self):
        self.assertEqual(
            self.test.ensemble_path,
            os.path.abspath(os.path.join("testfiles", "crest_conformers.xyz"))
        )