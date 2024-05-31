from censo.cli.cml_parser import parse
from censo.params import DESCR
from censo.ensembledata import EnsembleData
import shutil
import unittest
import os

os.chdir(os.path.split(__file__)[0])


test_args = parse(
    DESCR, "-inp testfiles/crest_conformers.xyz -chrg 0 -u 0".split())
test_dir = os.getcwd()


def getconfcount(path: str) -> int:
    with open(path, "r") as file:
        lines = file.readlines()

    nat = int(lines[0])
    return len(lines) // nat


class CoreTest(unittest.TestCase):
    def test_read_input_args(self):
        ensemble = EnsembleData(test_dir, args=test_args)
        ensemble.read_input(test_args.inp)
        nconf = getconfcount("testfiles/crest_conformers.xyz")
        self.assertEqual(nconf, len(ensemble.conformers))
        self.assertEqual(0, ensemble.runinfo["charge"])
        self.assertEqual(0, ensemble.runinfo["unpaired"])

    def test_read_input_script(self):
        ensemble = EnsembleData(test_dir)
        ensemble.read_input(test_args.inp, charge=2, unpaired=7)
        nconf = getconfcount("testfiles/crest_conformers.xyz")
        self.assertEqual(nconf, len(ensemble.conformers))
        self.assertEqual(2, ensemble.runinfo["charge"])
        self.assertEqual(7, ensemble.runinfo["unpaired"])

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
