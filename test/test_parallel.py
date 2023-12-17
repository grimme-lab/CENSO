import os
import random
import shutil
import unittest
from unittest.mock import patch

from censo.core import CensoCore
from censo.datastructure import ParallelJob
from censo.parallel import execute


class TestParallel(unittest.TestCase):
    def setUp(self):
        self.core = CensoCore(os.getcwd())
        self.core.read_input("testfiles/crest_conformers.xyz", charge=2, unpaired=7)

    @patch("censo.parallel.dqp")
    def test_execute(self, mock_dqp):
        mock_instructions = {
            "prog": "orca",
            "balance": True,
            "maxcores": 12,
            "jobtype": ["sp"],
            "omp": 4,
            "copy_mo": False,
            "retry_failed": True,
        }
        mock_dqp_results = self.__mock_dqp(mock_instructions)
        mock_dqp.return_value = mock_dqp_results
        print(
            f"Failed jobs: "
            f"{[job.conf.name for job in mock_dqp_results if not all(job.meta[jt]['success'] for jt in job.jobtype)]}"
        )
        failed_jobs = [
            i
            for i, job in enumerate(mock_dqp_results)
            if any(not job.meta[jt]["success"] for jt in job.jobtype)
        ]
        print(f"Failed jobs (indices): {failed_jobs}")

        execute(self.core.conformers, mock_instructions, os.getcwd())

    def __mock_dqp(self, instructions: dict) -> list[ParallelJob]:
        mock_dqp_results = []
        for conf in self.core.conformers:
            mock_dqp_results.append(
                ParallelJob(conf.geom, instructions["jobtype"], instructions["omp"])
            )

        for job in mock_dqp_results:
            for jt in job.jobtype:
                job.meta[jt] = {
                    "success": True if random.randint(0, 1) == 1 else False,
                    "error": "SCF not converged",
                }

        return mock_dqp_results

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
