import unittest
import os
import numpy as np

from censo_test.inputhandling import cml
from censo_test.cfg import DESCR
from censo_test.storage import CensoStorage
from censo_test.settings import InternalSettings

test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz".split())
test_dir = os.getcwd()

class SettingsTest(unittest.TestCase):
    def setUp(self):
        self.test = InternalSettings(CensoStorage(test_args, test_dir))
        
        
    def test_settings_current(self):
        print(self.test.settings_current())