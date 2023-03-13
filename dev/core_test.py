import unittest
import os
import numpy as np

from censo.inputhandling import cml
from censo.cfg import DESCR
from censo.core import CensoCore

test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz".split())
test_dir = os.getcwd()

class CoreTest(unittest.TestCase):
    def setUp(self):
        self.test = CensoCore(test_args, test_dir)

        
    def test_read_input(self):
        self.test.read_input()