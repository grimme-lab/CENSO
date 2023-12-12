import shutil
import unittest
import os
from collections import OrderedDict

from src import OrcaParser

test_dir = os.getcwd()


class OParserTest(unittest.TestCase):
    def setUp(self):
        self.test = OrcaParser()

    def test_read(self):
        inp = self.test.read_input(os.path.join(test_dir, "testfiles", "inp"))

        should = OrderedDict()
        should["main"] = ["RHF", "CCSD(T)", "def2-TZVP", "TightSCF"]
        should["paras"] = {"R=": ["4.0,0.5,35"]}
        should["geom"] = {}
        should["geom"]["def"] = ["xyz", "0", "1"]
        should["geom"]["coord"] = [
            ["H", "0", "0", "0"],
            ["F", "0", "0", "{R}"],
        ]

        self.assertDictEqual(inp, should)

    def test_write(self):
        towrite = OrderedDict()
        towrite["main"] = ["RHF", "CCSD(T)", "def2-TZVP", "TightSCF"]
        towrite["paras"] = {"R=": ["4.0,0.5,35"]}
        towrite["geom"] = {}
        towrite["geom"]["def"] = ["xyz", "0", "1"]
        towrite["geom"]["coord"] = [
            ["H", "0", "0", "0"],
            ["F", "0", "0", "{R}"],
        ]

        self.test.write_input(os.path.join(test_dir, "testfiles", "testinp"), towrite)
        with open(os.path.join(test_dir, "testfiles", "testinp"), "r") as file:
            written = file.readlines()

        with open(os.path.join(test_dir, "testfiles", "inp2"), "r") as file:
            should = file.readlines()

        self.assertListEqual(written, should)

    def test_read_template(self):
        inp = self.test.read_input(os.path.join(test_dir, "testfiles", "test.template"))

        should = OrderedDict()
        should["main"] = ["OPT"]
        should["geom"] = {}
        should["mp2"] = {"bla": ["bla"]}

        self.assertDictEqual(inp, should)

    def doCleanups(self):
        # perform cleanup
        delete = [
            "censo.log",
        ]
        for f in delete:
            f = os.path.join(os.getcwd(), f)
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)
