import unittest

from censo.prescreening import Prescreening
from censo.censo import startup

core, settings = startup("-inp testfiles/crest_conformers.xyz -solvent water".split())

class PartTest(unittest.TestCase):
    def setUp(self):
        self.test = Prescreening(core, settings)


    def something(self):
        print("do something")