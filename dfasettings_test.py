import unittest
import json

from censo.configuration import DfaHelper

dfas = json.load(open("./censo/assets/censo_dfa_settings.json", "r"))

class DfaSettingsTest(unittest.TestCase):
    def setUp(self):
        self.test = DfaHelper(dfas)
        
        
    def test_find_func(self):
        self.assertEqual(self.test.find_func())