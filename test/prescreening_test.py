import unittest

from src import Prescreening
from src.censo import startup

core, settings = startup("-inp testfiles/crest_conformers.xyz -solvent water".split())

class PartTest(unittest.TestCase):
    def setUp(self):
        self.test = Prescreening(core, settings)


    def something(self):
        print("do something")