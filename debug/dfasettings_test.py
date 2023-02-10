import unittest
import json

from censo_test.settings import DfaSettings

dfas = json.load(open("./censo_test/assets/censo_dfa_settings.json", "r"))

class DfaSettingsTest(unittest.TestCase):
    def setUp(self):
        self.test = DfaSettings(dfas)
        
        
    def test_find_func(self):
        self.assertEqual(self.test.find_func())