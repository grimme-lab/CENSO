import unittest

from src import DESCR
from src import cml


class CmlTest(unittest.TestCase):
    def setUp(self):
        self.test = "-inp crest_conformers.xyz".split()

    def test_cml(self):
        try:
            self.assertIsNotNone(cml(DESCR, self.test))
        except Exception:
            raise AssertionError

        print(cml(DESCR, self.test))
