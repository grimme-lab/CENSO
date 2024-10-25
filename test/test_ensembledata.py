from censo.cli.cml_parser import parse
from censo.params import DESCR
from censo.ensembledata import EnsembleData
import pytest


def test_read_input(self):
    # Read input via python instruction

    # Read input passed via cml args
    test_args = parse(argv="-i fixtures/crest_conformers.xyz".split())
    ensemble = EnsembleData(test_dir, args=test_args)
    ensemble.read_input(test_args.inp)
    nconf = 7
    assert nconf == len(ensemble.conformers)
    assert 0 == ensemble.runinfo["charge"]
    assert 0 == ensemble.runinfo["unpaired"]
