"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv implicitly and include the RRHO contribution.
"""
import os
from functools import reduce
from math import exp
from statistics import stdev

from ..datastructure import MoleculeData
from ..logging import setup_logger
from ..parallel import execute
from ..params import (AU2KCAL, GFNOPTIONS, GRIDOPTIONS, PLENGTH, PROGS,
                      SOLV_MODS)
from ..utilities import format_data, h1, print
from .prescreening import Prescreening

logger = setup_logger(__name__)


class Screening(Prescreening):
    _part_no = "1"

    _grid = "low+"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    # __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "threshold": {
            "default": 3.5
        },
        "func": {
            "default": "r2scan-3c"
        },
        "basis": {
            "default": "def2-TZVP"
        },
        "prog": {
            "default": "orca",
            "options": PROGS
        },
        "sm": {
            "default": "smd",
            "options": __solv_mods
        },
        "gfnv": {
            "default": "gfn2",
            "options": GFNOPTIONS
        },
        "run": {
            "default": True
        },
        "implicit": {
            "default": True
        },
        "template": {
            "default": False
        },
    }

    _settings = {}

    def optimize(self, cut: bool = True) -> None:
        """
        Advanced screening of the ensemble by doing single-point calculations on the input geometries,
        but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

        Basically consists of two parts:
            - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
            - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
        """
        super().optimize(cut=False)

        # NOTE: the following is only needed if 'evaluate_rrho' is enabled, since 'screening' runs the same procedure as prescreening before
        # therefore the sorting and filtering only needs to be redone if the rrho contributions are going to be included
        if self.get_general_settings()["evaluate_rrho"]:
            # PART (2)
            threshold = self.get_settings()["threshold"] / AU2KCAL

            jobtype = ["xtb_rrho"]
            prepinfo = self.setup_prepinfo(jobtype)

            # append results to previous results
            success, _, failed = execute(
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
                # calculate new gtot including RRHO contribution
                conf.results[self._name]["gtot"] = self.grrho(conf)

            # sort conformers list
            self.ensemble.conformers.sort(
                key=lambda conf: conf.results[self._name]["gtot"])

            if cut and len(self.ensemble.conformers) > 1:
                # calculate fuzzyness of threshold (adds 1 kcal/mol at max to the threshold)
                fuzzy = (1 / AU2KCAL) * (1 - exp(-AU2KCAL * stdev([
                    conf.results[self._name]["xtb_rrho"]["energy"]
                    for conf in self.ensemble.conformers
                ])))
                threshold += fuzzy
                print(
                    f"Updated fuzzy threshold: {threshold * AU2KCAL:.2f} kcal/mol."
                )

                # update the conformer list in ensemble (remove confs if below threshold)
                for confname in self.ensemble.update_conformers(
                        self.grrho, threshold):
                    print(f"No longer considering {confname}.")

            # calculate boltzmann weights from gtot values calculated here
            # trying to get temperature from instructions, set it to room temperature if that fails for some reason
            self.ensemble.calc_boltzmannweights(
                self.get_general_settings().get("temperature", 298.15),
                self._name)

            # if no conformers are filtered basically nothing happens

            # second 'write_results' for the updated sorting with RRHO contributions
            self.write_results2()

    def gsolv(self, conf: MoleculeData) -> float:
        """
        Override of the function from Prescreening.
        """
        # If solvation contributions should be included and the solvation free enthalpy
        # should not be included in the single-point energy the 'gsolv' job should've been run
        if not self.get_general_settings(
        )["gas-phase"] and not self.get_settings()["implicit"]:
            return conf.results[self._name]["gsolv"]["energy_solv"]
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
        return self.gsolv(conf) + conf.results[
            self._name]["xtb_rrho"]["energy"]

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
        print(h1(f"{self._name.upper()} SINGLE-POINT RESULTS"))
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
        if (all("prescreening" in conf.results.keys()
                for conf in self.ensemble.conformers)
                and not self.get_general_settings()["gas-phase"]):
            xtb_energies = {
                id(conf):
                conf.results["prescreening"]["xtb_gsolv"]["energy_xtb_gas"]
                for conf in self.ensemble.conformers
            }
            xtbmin = min(xtb_energies.values())
        else:
            xtb_energies = None

        # minimal total free enthalpy (single-point and potentially gsolv)
        gsolvmin = min(self.gsolv(conf) for conf in self.ensemble.conformers)

        # collect all dft single point energies
        dft_energies = ({
            id(conf): conf.results[self._name]["sp"]["energy"]
            for conf in self.ensemble.conformers
        } if not all("gsolv" in conf.results[self._name].keys()
                     for conf in self.ensemble.conformers) else {
                         id(conf):
                         conf.results[self._name]["gsolv"]["energy_gas"]
                         for conf in self.ensemble.conformers
                     })

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#":
            lambda conf: conf.name,
            "E (xTB)":
            lambda conf: f"{xtb_energies[id(conf)]:.6f}"
            if xtb_energies is not None else "---",
            "ΔE (xTB)":
            lambda conf: f"{(xtb_energies[id(conf)] - xtbmin) * AU2KCAL:.2f}"
            if xtb_energies is not None else "---",
            "E (DFT)":
            lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv":
            lambda conf: f"{self.gsolv(conf) - dft_energies[id(conf)]:.6f}"
            if "xtb_gsolv" in conf.results[self._name].keys() or "gsolv" in
            conf.results[self._name].keys() else "---",
            "Gtot":
            lambda conf: f"{self.gsolv(conf):.6f}",
            "ΔGtot":
            lambda conf: f"{(self.gsolv(conf) - gsolvmin) * AU2KCAL:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers]
                for conf in self.ensemble.conformers]

        lines = format_data(headers, rows, units=units)

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        print("".ljust(PLENGTH, "-"))

        # write everything to a file
        filename = f"{self._part_no}_{self._name.upper()}.out"
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, filename)}.")
        with open(os.path.join(self.ensemble.workdir, filename),
                  "w",
                  newline=None) as outfile:
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
        print(h1(f"{self._name.upper()} RRHO RESULTS"))

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
        if (all("prescreening" in conf.results.keys()
                for conf in self.ensemble.conformers)
                and not self.get_general_settings()["gas-phase"]):
            gxtb = {
                id(conf):
                conf.results["prescreening"]["xtb_gsolv"]["energy_xtb_gas"]
                for conf in self.ensemble.conformers
            }
            if self.get_general_settings()["evaluate_rrho"]:
                for conf in self.ensemble.conformers:
                    gxtb[id(conf)] += conf.results[self._name]["xtb_rrho"][
                        "gibbs"][self.get_general_settings()["temperature"]]
            gxtbmin = min(gxtb.values())
        else:
            gxtb = None

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(conf.results[self._name]["gtot"]
                      for conf in self.ensemble.conformers)

        # collect all dft single point energies
        dft_energies = ({
            id(conf): conf.results[self._name]["sp"]["energy"]
            for conf in self.ensemble.conformers
        } if not all("gsolv" in conf.results[self._name].keys()
                     for conf in self.ensemble.conformers) else {
                         id(conf):
                         conf.results[self._name]["gsolv"]["energy_gas"]
                         for conf in self.ensemble.conformers
                     })

        printmap = {
            "CONF#":
            lambda conf: conf.name,
            "G (xTB)":
            lambda conf: f"{gxtb[id(conf)]:.6f}"
            if gxtb is not None else "---",
            "ΔG (xTB)":
            lambda conf: f"{(gxtb[id(conf)] - gxtbmin) * AU2KCAL:.2f}"
            if gxtb is not None else "---",
            "E (DFT)":
            lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv":
            lambda conf: f"{self.gsolv(conf) - dft_energies[id(conf)]:.6f}"
            if not self.get_settings().get("implicit", False) else "---",
            "GmRRHO":
            lambda conf:
            f"{conf.results[self._name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
            if self.get_general_settings()["evaluate_rrho"] else "---",
            "Gtot":
            lambda conf: f"{conf.results[self._name]['gtot']:.6f}",
            "ΔGtot":
            lambda conf:
            f"{(conf.results[self._name]['gtot'] - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers]
                for conf in self.ensemble.conformers]

        lines = format_data(headers, rows, units=units)

        # list the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):\n"
        )
        lines.append(
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n"
        )

        # calculate averaged free enthalpy
        avG = sum([
            conf.results[self._name]["bmw"] * conf.results[self._name]["gtot"]
            for conf in self.ensemble.conformers
        ])

        # calculate averaged free energy
        avE = sum([
            conf.results[self._name]["bmw"] *
            conf.results[self._name]["sp"]["energy"]
            for conf in self.ensemble.conformers
        ])

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part1==\n"
        )
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # append lines to already existing file
        filename = f"{self._part_no}_{self._name.upper()}.out"
        with open(os.path.join(self.ensemble.workdir, filename),
                  "a",
                  newline=None) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results to a json file
        self.write_json()
