from functools import reduce
import os

from .screening import Screening
from .prescreening import Prescreening
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
    AU2KCAL,
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
    def run(self, cut: bool = True) -> None:
        """
        """
        Prescreening.run(self, cut=False)

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
                    conf.results[self._name]["xtb_rrho"] = conf.results["optimization"]["xtb_rrho"]
                    conf.results[self._name]["gtot"] = conf.results["optimization"]["gtot"]

        # sort conformers list
        self.ensemble.conformers.sort(
            key=lambda conf: conf.results[self._name]["gtot"])

        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        if cut:
            # Get Boltzmann population threshold from settings
            threshold = self.get_settings()["threshold"]

            # Update ensemble using Boltzman population threshold
            for confname in self.ensemble.update_conformers(lambda conf: conf.results[self._name]["bmw"], threshold, boltzmann=True):
                print(f"No longer considering {confname}.")

        # second 'write_results' for the updated sorting with RRHO contributions
        self.write_results2()

        # dump ensemble
        self.ensemble.dump_ensemble(self._name)

        # DONE

    def write_results2(self) -> None:
        """
        Additional write function in case RRHO is used.
        Write the results to a file in formatted way. This is appended to the first file.
        writes (2):
            G (xtb),
            δG (xtb),
            E (DFT),
            δGsolv (DFT),
            Grrho,
            Gtot,
            δGtot

        Also writes them into an easily digestible format.
        """
        # column headers
        headers = [
            "CONF#",
            "E (DFT)",
            "ΔGsolv",
            "GmRRHO",
            "Gtot",
            "ΔGtot",
            "Boltzmann weight"
        ]

        # column units
        units = [
            "",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
            f"% at {self.get_general_settings().get('temperature', 298.15)} K",
        ]

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(conf.results[self._name]["gtot"]
                      for conf in self.ensemble.conformers)

        # collect all dft single point energies
        dft_energies = (
            {
                id(conf): conf.results[self._name]["sp"]["energy"]
                for conf in self.ensemble.conformers
            }
            if not all(
                "gsolv" in conf.results[self._name].keys()
                for conf in self.ensemble.conformers
            )
            else {
                id(conf): conf.results[self._name]["gsolv"]["energy_gas"]
                for conf in self.ensemble.conformers
            }
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT)": lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv": lambda conf: f"{self.gtot(conf) - dft_energies[id(conf)]:.6f}"
            if not self.get_settings().get("implicit", False)
            else "---",
            "GmRRHO": lambda conf: f"{conf.results[self._name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
            if self.get_general_settings()["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{conf.results[self._name]['gtot']:.6f}",
            "ΔGtot": lambda conf: f"{(conf.results[self._name]['gtot'] - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # append lines to already existing file
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, f'{self._name}.out')}."
        )
        with open(
            os.path.join(self.ensemble.workdir, f"{self._name}.out"), "a", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results to a json file
        self.write_json()
