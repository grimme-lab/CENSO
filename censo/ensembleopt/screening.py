"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv implicitly and include the RRHO contribution.
"""
import os
from math import isclose, exp
from statistics import stdev
from functools import reduce

from censo.ensembleopt.prescreening import Prescreening
from censo.part import CensoPart
from censo.utilities import print, timeit, format_data
from censo.parallel import ProcessHandler, execute
from censo.core import CensoCore
from censo.datastructure import MoleculeData
from censo.params import (
    AU2KCAL,
    SOLV_MODS,
    GSOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
)
from censo.utilities import DfaHelper


class Screening(Prescreening):
    alt_name = "part1"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    # __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "threshold": {
            "default": 3.5,
            "range": [
                0.75,
                7.5
            ]
        },
        "func": {
            "default": "r2scan-3c",
            "options": DfaHelper.find_func("screening")
        },
        "basis": {
            "default": "def2-TZVP",
            "options": BASIS_SETS
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
        "grid": {
            "default": "low+",
            "options": GRIDOPTIONS
        },
        "run": {
            "default": True
        },
        "gcp": {
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

    @timeit
    # @CensoPart._create_dir - not required here because super().run() already does this
    def run(self) -> None:
        """
        Advanced screening of the ensemble by doing single-point calculations on the input geometries,
        but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

        Basically consists of two parts:
            - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
            - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
        """
        super().run()

        # NOTE: the following is only needed if 'evaluate_rrho' is enabled, since 'screening' runs the same procedure as prescreening before
        # therefore the sorting and filtering only needs to be redone if the rrho contributions are going to be included
        if self._instructions["evaluate_rrho"]:
            # PART (2)
            threshold = self._instructions["threshold"] / AU2KCAL

            self._instructions["jobtype"] = ["xtb_rrho"]

            # append results to previous results
            results = execute(self.core.conformers, self._instructions, self.dir)
            for conf in self.core.conformers:
                # update results for each conformer
                conf.results[self._name].update(results[id(conf)])

                # calculate new gtot including RRHO contribution
                conf.results[self._name]["gtot"] = self.grrho(conf)

            # sort conformers list
            self.core.conformers.sort(key=lambda conf: conf.results[self._name]["gtot"])

            # pick the free enthalpy of the lowest conformer
            limit = min([conf.results[self._name]["gtot"] for conf in self.core.conformers])

            # calculate fuzzyness of threshold (adds 1 kcal/mol at max to the threshold)
            fuzzy = (1 / AU2KCAL) * (1 - exp(-5 * stdev(
                [conf.results[self._name]["xtb_rrho"]["energy"] for conf in
                 self.core.conformers]) * AU2KCAL))
            threshold += fuzzy
            print(f"Updated fuzzy threshold: {threshold * AU2KCAL:.2f} kcal/mol.")

            # filter out all conformers above threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.grrho(x) - limit > threshold,
                    self.core.conformers
                )
            ]

            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(filtered)

            # calculate boltzmann weights from gtot values calculated here
            # trying to get temperature from instructions, set it to room temperature if that fails for some reason
            self.core.calc_boltzmannweights(
                self._instructions.get("temperature", 298.15),
                self._name
            )

            # if no conformers are filtered basically nothing happens

            # second 'write_results' for the updated sorting with RRHO contributions
            self.write_results2()

        # dump ensemble
        self.core.dump_ensemble(self._name)

        # DONE

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
            "[kcal/mol]",
            "[Eh]",
            "[kcal/mol]",
        ]

        # variables for printmap
        # minimal xtb single-point energy (taken from prescreening)
        # TODO - where do prescreening and screening get xtb single-point from? - if prescreening is done grab from there, otherwise xtb_sp should be run
        xtbmin = min(
            conf.results["prescreening"]['xtb_gsolv']['energy_xtb_gas']
            for conf in self.core.conformers
        )

        # minimal total free enthalpy (single-point and potentially gsolv)
        gtotmin = min(self.gtot(conf) for conf in self.core.conformers)

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: f"{conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas']:.6f}",  # TODO
            "ΔE (xTB)": lambda conf: f"{(conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}",  # TODO
            "E (DFT)": lambda conf: f"{conf.results[self._name]['sp']['energy']:.6f}",
            "ΔGsolv": lambda conf:
            f"{self.gtot(conf) - conf.results[self._name]['sp']['energy']:.6f}"
            if "xtb_gsolv" in conf.results[self._name].keys() or "gsolv" in conf.results[
                self._name].keys()
            else "---",
            "Gtot": lambda conf: f"{self.gtot(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.gtot(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # write everything to a file
        with open(os.path.join(self.core.workdir, f"{self._name}.out"), "w",
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
        # TODO - what if there was no prescreening?
        gxtbmin = min(
            conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] +
            conf.results[self._name]['xtb_rrho']['gibbs'][
                self._instructions["temperature"]]  # TODO?
            if self._instructions["evaluate_rrho"] else conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas']
            for conf in self.core.conformers
        )

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            self.grrho(conf)
            for conf in self.core.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "G (xTB)": lambda conf: f"{conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self._name]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}",
            # TODO
            "ΔG (xTB)": lambda conf: f"{(conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self._name]['xtb_rrho']['gibbs'][self._instructions['temperature']] - gxtbmin) * AU2KCAL:.2f}",
            # TODO
            "E (DFT)": lambda conf: f"{conf.results[self._name]['sp']['energy']:.6f}",
            "ΔGsolv": lambda conf:
            f"{self.gtot(conf) - conf.results[self._name]['sp']['energy']:.6f}"
            if not isclose(self.gtot(conf), conf.results[self._name]['sp']['energy'])
            else "---",
            "GmRRHO": lambda conf:
            f"{conf.results[self._name]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}"
            if self._instructions["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # append lines to already existing file
        with open(os.path.join(self.core.workdir, f"{self._name}.out"), "a",
                  newline=None) as outfile:
            outfile.writelines(lines)
