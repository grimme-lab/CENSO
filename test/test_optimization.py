import os
import random
import shutil
import unittest
from unittest.mock import patch
from pprint import pprint

from censo.ensembledata import EnsembleData
from censo.ensembleopt.optimization import Optimization


class TestOptimization(unittest.TestCase):
    @patch("censo.ensembleopt.optimization.execute")
    def test_run(self, mock_execute):
        ensemble = EnsembleData(os.getcwd())
        ensemble.read_input(
            "testfiles/crest_conformers.xyz", charge=2, unpaired=7)

        # Create an instance of the Optimization class
        optimization = Optimization(ensemble)

        # Mock execution
        mock_results = {
            id(conf): {"xtb_opt": {}, "xtb_rrho": {}} for conf in ensemble.conformers
        }
        mock_results_2 = []
        for conf in ensemble.conformers:
            mock_results[id(conf)]["xtb_opt"][
                "energy"
            ] = -1396.397775 + random.normalvariate(0, 0.001)
            mock_results[id(conf)]["xtb_opt"]["ecyc"] = [-1396.387775 -
                                                         i * random.normalvariate(0, 0.0001) for i in range(7)]
            mock_results[id(conf)]["xtb_opt"]["converged"] = True
            mock_results[id(conf)]["xtb_opt"]["cycles"] = 7
            mock_results[id(conf)]["xtb_opt"]["geom"] = conf.geom.xyz
            mock_results[id(conf)]["xtb_opt"][
                "grad_norm"
            ] = 0.0001 + random.normalvariate(0, 0.0003)
            mock_results[id(conf)]["xtb_opt"]["gncyc"] = [
                0.001 - i * random.normalvariate(0, 0.0001) for i in range(7)]
            mock_results[id(conf)]["xtb_rrho"]["energy"] = 0.001
            mock_results[id(conf)]["xtb_rrho"].setdefault("gibbs", {}).setdefault(
                optimization.get_general_settings()["temperature"], 0.001
            )

        mock_execute.return_value = (mock_results, mock_results_2)

        # Call the run method
        optimization.run()

    def doCleanups(self):
        # perform cleanup
        delete = [
            "optimization",
            "censo.log",
            "censo_ensemble_optimization.xyz",
            "optimization.out",
            "optimization.json",
        ]
        for f in delete:
            f = os.path.join(os.getcwd(), f)
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)


if __name__ == "__main__":
    unittest.main()
