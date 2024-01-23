from functools import reduce

from .screening import Screening
from .prescreening import Prescreening
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
)
from ..utilities import DfaHelper, setup_logger
from ..utilities import print, timeit, format_data

logger = setup_logger(__name__)


class Refinement(Screening):
    alt_name = "part3"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    # __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "threshold": {"default": 0.95, "range": [0.01, 0.99]},
        "func": {"default": "wb97x-d3", "options": DfaHelper.find_func("refinement")},
        "basis": {"default": "def2-TZVP(-f)", "options": BASIS_SETS},
        "prog": {"default": "orca", "options": PROGS},
        "sm": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
        "grid": {"default": "high+", "options": GRIDOPTIONS},
        "run": {"default": True},
        "gcp": {"default": True},
        "implicit": {"default": True},
        "template": {"default": False},
    }

    _settings = {}

    @timeit
    # @CensoPart._create_dir - not required here because Prescreening.run(self) already does this
    def run(self) -> None:
        """
        """
        Prescreening.run(self)

        if self.get_general_settings()["evaluate_rrho"]:
            # Check if evaluate_rrho, then check if optimization was run and use that value, otherwise do xtb_rrho
            if not all("optimization" in conf.results.keys() for conf in self.ensemble.conformers):
                jobtype = ["xtb_rrho"]
                prepinfo = self.setup_prepinfo(jobtype)

                # append results to previous results
                results, failed = execute(
                    self.ensemble.conformers,
                    self.dir,
                    self.get_settings()["prog"],
                    prepinfo,
                    jobtype,
                    copy_mo=self.get_general_settings()["copy_mo"],
                    balance=self.get_general_settings()["balance"],
                    omp=self.get_general_settings()["omp"],
                    maxcores=self.get_general_settings()["maxcores"],
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Remove failed conformers
                self.ensemble.remove_conformers(failed)

                for conf in self.ensemble.conformers:
                    # update results for each conformer
                    conf.results[self._name].update(results[id(conf)])

                    # calculate new gtot including RRHO contribution
                    conf.results[self._name]["gtot"] = self.grrho(conf)
            else:
                # Use values from optimization rrho
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["gtot"] = conf.results["optimization"]["gtot"]

        # sort conformers list
        self.ensemble.conformers.sort(
            key=lambda conf: conf.results[self._name]["gtot"])

        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        # Get Boltzmann population threshold from settings
        threshold = self.get_settings()["threshold"]

        # Update ensemble using Boltzman population threshold
        self.ensemble.update_conformers(
            lambda conf: conf.results[self._name]["gtot"], threshold, boltzmann=True)

        # second 'write_results' for the updated sorting with RRHO contributions
        self.write_results2()

        # dump ensemble
        self.ensemble.dump_ensemble(self._name)

        # DONE
