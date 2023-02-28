import unittest
import os
import numpy as np

from censo_test.inputhandling import cml
from censo_test.cfg import DESCR
from censo_test.core import CensoCore
from censo_test.storage import CensoStorage

test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz".split())
test_dir = os.getcwd()

class CoreTest(unittest.TestCase):
    def setUp(self):
        self.test = CensoCore.factory(CensoStorage(test_args, test_dir))
        
    
    @CensoCore.check_instance
    def test_instances(self):
        self.assertEqual(np.sum([1, 5, 6]), 12)
        
        
    def test_core(self):
        self.assertEqual(CensoCore.core(), self.test)
        
        
    def test_read_input(self):
        self.test.read_input()