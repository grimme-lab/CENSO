import os
import random
import shutil
import unittest
from unittest.mock import patch

from src.core import CensoCore
from src.ensembleopt.optimization import Optimization


class TestOptimization(unittest.TestCase):

    @patch("src.ensembleopt.optimization.execute")
    def test_run(self, mock_execute):
        core = CensoCore(os.getcwd())
        core.read_input("testfiles/crest_conformers.xyz", charge=2, unpaired=7)

        # Create an instance of the Optimization class
        optimization = Optimization(core)

        # Mock execution
        mock_results = {id(conf): {"xtb_opt": {}, "xtb_rrho": {}} for conf in core.conformers}
        for conf in core.conformers:
            mock_results[id(conf)]["xtb_opt"]["energy"] = -1396.397775 + random.normalvariate(0, 0.1)
            mock_results[id(conf)]["xtb_opt"]["converged"] = True
            mock_results[id(conf)]["xtb_opt"]["cycles"] = 8
            mock_results[id(conf)]["xtb_opt"]["geom"] = conf.geom.xyz
            mock_results[id(conf)]["xtb_opt"]["grad_norm"] = 0.0001 + random.normalvariate(0, 0.001)
            mock_results[id(conf)]["xtb_rrho"]["energy"] = 0.001
            mock_results[id(conf)]["xtb_rrho"].setdefault("gibbs", {}).setdefault(optimization._instructions["temperature"], 0.001)

        mock_execute.return_value = mock_results

        # Call the run method
        optimization.run()

    def doCleanups(self):
        # perform cleanup
        delete = [
            "optimization",
            "censo.log",
            "censo_ensemble_optimization.xyz",
            "optimization.out",
        ]
        for f in delete:
            f = os.path.join(os.getcwd(), f)
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)


if __name__ == '__main__':
    unittest.main()
