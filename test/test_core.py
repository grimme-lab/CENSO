import shutil
import unittest
import os

from ..src.censo.inputhandling import cml
from ..src.censo.params import DESCR
from ..src.censo.core import CensoCore

test_args = cml(DESCR, "-inp testfiles/crest_conformers.xyz -chrg 0 -u 0".split())
test_dir = os.getcwd()


def getconfcount(path: str) -> int:
    with open(path, "r") as file:
        lines = file.readlines()

    nat = int(lines[0])
    return len(lines) // nat


class CoreTest(unittest.TestCase):

    def test_read_input_args(self):
        core = CensoCore(test_dir, args=test_args)
        core.read_input(test_args.inp)
        nconf = getconfcount("testfiles/crest_conformers.xyz")
        self.assertEqual(nconf, len(core.conformers))
        self.assertEqual(0, core.runinfo["charge"])
        self.assertEqual(0, core.runinfo["unpaired"])

    def test_read_input_script(self):
        core = CensoCore(test_dir)
        core.read_input(test_args.inp, charge=0, unpaired=0)
        nconf = getconfcount("testfiles/crest_conformers.xyz")
        self.assertEqual(nconf, len(core.conformers))
        self.assertEqual(0, core.runinfo["charge"])
        self.assertEqual(0, core.runinfo["unpaired"])

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
