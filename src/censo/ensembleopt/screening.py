"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv implicitly and include the RRHO contribution.
"""
import os
from functools import reduce
from math import isclose, exp
from statistics import stdev

from ..datastructure import MoleculeData
from .prescreening import Prescreening
from ..parallel import execute
from ..params import (
    AU2KCAL,
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
)
from ..utilities import DfaHelper, setup_logger
from ..utilities import print, timeit, format_data

logger = setup_logger(__name__)


class Screening(Prescreening):
    alt_name = "part1"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    # __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "threshold": {"default": 3.5, "range": [0.75, 7.5]},
        "func": {"default": "r2scan-3c", "options": DfaHelper.find_func("screening")},
        "basis": {"default": "def2-TZVP", "options": BASIS_SETS},
        "prog": {"default": "orca", "options": PROGS},
        "sm": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
        "grid": {"default": "low+", "options": GRIDOPTIONS},
        "run": {"default": True},
        "gcp": {"default": True},
        "implicit": {"default": True},
        "template": {"default": False},
    }

    _settings = {}

    @timeit
    # @CensoPart._create_dir - not required here because super().run() already does this
    def run(self, cut: bool = True) -> None:
        """
        Advanced screening of the ensemble by doing single-point calculations on the input geometries,
        but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

        Basically consists of two parts:
            - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
            - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
        """
        super().run(cut=cut)

        # NOTE: the following is only needed if 'evaluate_rrho' is enabled, since 'screening' runs the same procedure as prescreening before
        # therefore the sorting and filtering only needs to be redone if the rrho contributions are going to be included
        if self.get_general_settings()["evaluate_rrho"]:
            # PART (2)
            threshold = self.get_settings()["threshold"] / AU2KCAL

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

            # sort conformers list
            self.ensemble.conformers.sort(
                key=lambda conf: conf.results[self._name]["gtot"])

            if cut:
                # calculate fuzzyness of threshold (adds 1 kcal/mol at max to the threshold)
                fuzzy = (1 / AU2KCAL) * (1 - exp(-AU2KCAL * stdev(
                    [conf.results[self._name]["xtb_rrho"]["energy"]
                        for conf in self.ensemble.conformers]
                )))
                threshold += fuzzy
                print(
                    f"Updated fuzzy threshold: {threshold * AU2KCAL:.2f} kcal/mol.")

                # update the conformer list in ensemble (remove conf if below threshold)
                for confname in self.ensemble.update_conformers(self.grrho, threshold):
                    print(f"No longer considering {confname}.")

            # calculate boltzmann weights from gtot values calculated here
            # trying to get temperature from instructions, set it to room temperature if that fails for some reason
            self.ensemble.calc_boltzmannweights(
                self.get_general_settings().get("temperature", 298.15), self._name
            )

            # if no conformers are filtered basically nothing happens

            # second 'write_results' for the updated sorting with RRHO contributions
            self.write_results2()

        # dump ensemble
        self.ensemble.dump_ensemble(self._name)

        # DONE

    def gtot(self, conf: MoleculeData) -> float:
        """
        Override of the function from Prescreening.
        """
        # If solvation contributions should be included and the solvation free enthalpy
        # should not be included in the single-point energy the 'gsolv' job should've been run
        if (
                not self.get_general_settings()["gas-phase"]
                and not self.get_settings()["implicit"]
        ):
            return conf.results[self._name]["gsolv"]["energy_gas"]
            + conf.results[self._name]["gsolv"]["gsolv"]
        # Otherwise, return just the single-point energy
        else:
            return conf.results[self._name]["sp"]["energy"]

    def grrho(self, conf: MoleculeData) -> float:
        """
        Calculate the total Gibbs free energy (Gtot) of a given molecule using DFT single-point and gsolv (if included) and RRHO contributions.

        Parameters:
            conf (MoleculeData): The MoleculeData object containing the necessary information for the calculation.

        Returns:
            float: The total Gibbs free energy (Gtot) of the molecule.
        """
        # Gtot = E(DFT) + Gsolv + Grrho
        # NOTE: grrho should only be called if evaluate_rrho is True
        return self.gtot(conf) + conf.results[self._name]["xtb_rrho"]["energy"]

    def write_results(self) -> None:
        """
        Overrides the write_results function of Prescreening.
        Write the results to a file in formatted way.
        writes (1):
            E (xtb),
            δE (xtb),
            E (DFT),
            δGsolv (DFT),
            Gtot,
            δGtot

        Generates NO csv file. All info is included in the file written in write_results2.
        """
        # PART (1) of writing
        # column headers
        headers = [
            "CONF#",
            "E (xTB)",
            "ΔE (xTB)",
            "E (DFT)",
            "ΔGsolv",
            "Gtot",
            "ΔGtot",
        ]

        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
        ]

        # variables for printmap
        # minimal xtb single-point energy (taken from prescreening)
        if (
            all("prescreening" in conf.results.keys()
                for conf in self.ensemble.conformers)
            and not self.get_general_settings()["gas-phase"]
        ):
            xtb_energies = {
                id(conf): conf.results["prescreening"]["xtb_gsolv"]["energy_xtb_gas"]
                for conf in self.ensemble.conformers
            }
            xtbmin = min(xtb_energies.values())
        else:
            xtb_energies = None

        # minimal total free enthalpy (single-point and potentially gsolv)
        gtotmin = min(self.gtot(conf) for conf in self.ensemble.conformers)

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

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: f"{xtb_energies[id(conf)]:.6f}"
            if xtb_energies is not None
            else "---",
            "ΔE (xTB)": lambda conf: f"{(xtb_energies[id(conf)] - xtbmin) * AU2KCAL:.2f}"
            if xtb_energies is not None
            else "---",
            "E (DFT)": lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv": lambda conf: f"{self.gtot(conf) - dft_energies[id(conf)]:.6f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            or "gsolv" in conf.results[self._name].keys()
            else "---",
            "Gtot": lambda conf: f"{self.gtot(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.gtot(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # write everything to a file
        with open(
            os.path.join(self.ensemble.workdir, f"{self._name}.out"), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

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
            "G (xTB)",
            "ΔG (xTB)",
            "E (DFT)",
            "ΔGsolv",
            "GmRRHO",
            "Gtot",
            "ΔGtot",
        ]

        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
        ]

        # minimal xtb energy from single-point (and mRRHO)
        if (
            all("prescreening" in conf.results.keys()
                for conf in self.ensemble.conformers)
            and not self.get_general_settings()["gas-phase"]
        ):
            gxtb = {
                id(conf): conf.results["prescreening"]["xtb_gsolv"]["energy_xtb_gas"]
                for conf in self.ensemble.conformers
            }
            if self.get_general_settings()["evaluate_rrho"]:
                for conf in self.ensemble.conformers:
                    gxtb[id(conf)] += conf.results[self._name]["xtb_rrho"]["gibbs"][
                        self.get_general_settings()["temperature"]
                    ]
            gxtbmin = min(gxtb.values())
        else:
            gxtb = None

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        if self.get_general_settings()["evaluate_rrho"]:
            gtotmin = min(self.grrho(conf)
                          for conf in self.ensemble.conformers)
        else:
            gtotmin = min(self.gtot(conf) for conf in self.ensemble.conformers)

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
            "G (xTB)": lambda conf: f"{gxtb[id(conf)]:.6f}"
            if gxtb is not None
            else "---",
            "ΔG (xTB)": lambda conf: f"{(gxtb[id(conf)] - gxtbmin) * AU2KCAL:.2f}"
            if gxtb is not None
            else "---",
            "E (DFT)": lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv": lambda conf: f"{self.gtot(conf) - dft_energies[id(conf)]:.6f}"
            if not self.get_settings().get("implicit", False)
            else "---",
            "GmRRHO": lambda conf: f"{conf.results[self._name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
            if self.get_general_settings()["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

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
