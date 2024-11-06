"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv implicitly and include the RRHO contribution.
"""

import os
from math import exp
from statistics import stdev

from ..datastructure import MoleculeData
from ..logging import setup_logger
from ..parallel import execute
from ..params import AU2KCAL, PLENGTH, Config
from ..utilities import format_data, h1, print, DfaHelper, Factory
from .prescreening import Prescreening

logger = setup_logger(__name__)


class Screening(Prescreening):
    """
    Advanced screening of the ensemble by doing single-point calculations on the input geometries,
    but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

    Basically consists of two parts:
        - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
        - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
    """

    _grid = "low+"

    __solv_mods = {prog: Config.SOLV_MODS[prog] for prog in Config.PROGS}
    # __gsolv_mods = reduce(lambda x, y: x + y, GConfig.SOLV_MODS.values())

    _options = {
        "threshold": {"default": 3.5},
        "func": {
            "default": "r2scan-3c",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in Config.PROGS},
        },
        "basis": {"default": "def2-TZVP"},
        "prog": {"default": "tm", "options": Config.PROGS},
        "sm": {"default": "cosmors", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": Config.GFNOPTIONS},
        "run": {"default": True},
        "implicit": {"default": False},
        "template": {"default": False},
    }

    _settings = {}

    def _optimize(self, cut: bool = True) -> None:
        """
        Screening ensemble optimization logic. Basically run prescreening without cutting down the ensemble
        first, then calculate mrrho contributions if necessary.
        """
        super()._optimize(cut=False)

        # NOTE: the following is only needed if 'evaluate_rrho' is enabled, since 'screening' runs the same procedure as prescreening before
        # therefore the sorting and filtering only needs to be redone if the rrho contributions are going to be included
        if self.get_general_settings()["evaluate_rrho"]:
            # PART (2)
            threshold = self.get_settings()["threshold"] / AU2KCAL

            jobtype = ["xtb_rrho"]
            prepinfo = self._setup_prepinfo(jobtype)

            # append results to previous results
            results, failed = execute(
                self._ensemble.conformers,
                self._dir,
                self.get_settings()["prog"],
                prepinfo,
                jobtype,
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                retry_failed=self.get_general_settings()["retry_failed"],
            )

            # Remove failed conformers
            self._ensemble.remove_conformers(failed)

            # Update results
            self._update_results(results)

            for conf in self._ensemble.conformers:
                # calculate new gtot including RRHO contribution
                self.data["results"][conf.name]["gtot"] = self._grrho(conf)

            if cut and len(self._ensemble.conformers) > 1:
                # calculate fuzzyness of threshold (adds 1 kcal/mol at max to the threshold)
                fuzzy = (1 / AU2KCAL) * (
                    1
                    - exp(
                        -AU2KCAL
                        * stdev(
                            [
                                self.data["results"][conf.name]["xtb_rrho"]["energy"]
                                for conf in self._ensemble.conformers
                            ]
                        )
                    )
                )
                threshold += fuzzy
                print(f"Updated fuzzy threshold: {threshold * AU2KCAL:.2f} kcal/mol.")

                limit = min(self._grrho(conf) for conf in self._ensemble.conformers)
                filtered = list(
                    filter(
                        lambda conf: self._grrho(conf) - limit > threshold,
                        self._ensemble.conformers,
                    )
                )

                # update the conformer list in ensemble (remove confs if below threshold)
                self._ensemble.remove_conformers([conf.name for conf in filtered])
                for confname in filtered:
                    print(f"No longer considering {confname}.")

            # calculate boltzmann weights from gtot values calculated here
            # trying to get temperature from instructions, set it to room temperature if that fails for some reason
            self._update_results(self._calc_boltzmannweights())

    def _gsolv(self, conf: MoleculeData) -> float:
        """
        Override of the function from Prescreening.
        """
        # If solvation contributions should be included and the solvation free enthalpy
        # should not be included in the single-point energy the 'gsolv' job should've been run
        if "gsolv" in self.data["results"][conf.name]:
            return self.data["results"][conf.name]["gsolv"]["energy_solv"]
        # Otherwise, return just the single-point energy
        else:
            return self.data["results"][conf.name]["sp"]["energy"]

    def _grrho(self, conf: MoleculeData) -> float:
        """
        Calculate the total Gibbs free energy (Gtot) of a given molecule using DFT single-point and gsolv (if included) and RRHO contributions.

        Parameters:
            conf (MoleculeData): The MoleculeData object containing the necessary information for the calculation.

        Returns:
            float: The total Gibbs free energy (Gtot) of the molecule.
        """
        # Gtot = E(DFT) + Gsolv + Grrho
        # NOTE: grrho should only be called if evaluate_rrho is True
        return self._gsolv(conf) + self.data["results"][conf.name]["xtb_rrho"]["energy"]

    # Unused atm
    # def _write_results(self) -> None:
    #    """
    #    Similar to the _write_results from prescreening.
    #    Write the results to a file in formatted way.
    #    writes (1):
    #        E (xtb),
    #        δE (xtb),
    #        E (DFT),
    #        δGsolv (DFT),
    #        Gtot,
    #        δGtot

    #    Generates NO csv file. All info is included in the file written in write_results2.
    #    """
    #    print(h1(f"{self.name.upper()} SINGLE-POINT RESULTS"))
    #    # PART (1) of writing
    #    # column headers
    #    headers = [
    #        "CONF#",
    #        "E (xTB)",
    #        "ΔE (xTB)",
    #        "E (DFT)",
    #        "ΔGsolv",
    #        "Gtot",
    #        "ΔGtot",
    #    ]

    #    # column units
    #    units = [
    #        "",
    #        "[Eh]",
    #        "[kcal/mol]",
    #        "[Eh]",
    #        "[Eh]",
    #        "[Eh]",
    #        "[kcal/mol]",
    #    ]

    #    # variables for printmap
    #    # minimum xtb single-point energy (taken from prescreening)
    #    xtb_energies = None
    #    xtbmin = None
    #    if (
    #        any(type(p) is Prescreening for p in self._ensemble.results)
    #        and not self.get_general_settings()["gas-phase"]
    #    ):
    #        # Get the most recent prescreening part
    #        using_part = [
    #            p for p in self._ensemble.results if type(p) is Prescreening
    #        ][-1]

    #        xtb_energies = {
    #            conf.name: using_part.data["results"][conf.name]["xtb_gsolv"][
    #                "energy_xtb_gas"
    #            ]
    #            for conf in self._ensemble.conformers
    #        }
    #        xtbmin = min(xtb_energies.values())

    #    # minimum total free enthalpy (single-point and potentially gsolv)
    #    gsolvmin = min(self._gsolv(conf) for conf in self._ensemble.conformers)

    #    # collect all dft single point energies
    #    dft_energies = (
    #        {
    #            conf.name: self.data["results"][conf.name]["sp"]["energy"]
    #            for conf in self._ensemble.conformers
    #        }
    #        if not all(
    #            "gsolv" in self.data["results"][conf.name].keys()
    #            for conf in self._ensemble.conformers
    #        )
    #        else {
    #            conf.name: self.data["results"][conf.name]["gsolv"]["energy_gas"]
    #            for conf in self._ensemble.conformers
    #        }
    #    )

    #    # determines what to print for each conformer in each column
    #    printmap = {
    #        "CONF#": lambda conf: conf.name,
    #        "E (xTB)": lambda conf: (
    #            f"{xtb_energies[conf.name]:.6f}" if xtb_energies is not None else "---"
    #        ),
    #        "ΔE (xTB)": lambda conf: (
    #            f"{(xtb_energies[conf.name] - xtbmin) * AU2KCAL:.2f}"
    #            if xtb_energies is not None
    #            else "---"
    #        ),
    #        "E (DFT)": lambda conf: f"{dft_energies[conf.name]:.6f}",
    #        "ΔGsolv": lambda conf: (
    #            f"{self._gsolv(conf) - dft_energies[conf.name]:.6f}"
    #            if "xtb_gsolv" in self.data["results"][conf.name].keys()
    #            or "gsolv" in self.data["results"][conf.name].keys()
    #            else "---"
    #        ),
    #        "Gtot": lambda conf: f"{self._gsolv(conf):.6f}",
    #        "ΔGtot": lambda conf: f"{(self._gsolv(conf) - gsolvmin) * AU2KCAL:.2f}",
    #    }

    #    rows = [
    #        [printmap[header](conf) for header in headers]
    #        for conf in self._ensemble.conformers
    #    ]

    #    lines = format_data(headers, rows, units=units)

    #    # Print everything
    #    for line in lines:
    #        print(line, flush=True, end="")

    #    print("".ljust(PLENGTH, "-"))

    #    # write everything to a file
    #    filename = f"{self._part_nos[self.name]}_{self.name.upper()}.out"
    #    logger.debug(f"Writing to {os.path.join(os.getcwd(), filename)}.")
    #    with open(os.path.join(os.getcwd(), filename), "w", newline=None) as outfile:
    #        outfile.writelines(lines)

    def _write_results(self) -> None:
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
        print(h1(f"{self.name.upper()} SINGLE-POINT (+ mRRHO) RESULTS"))

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
            "Boltzmann weight",
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
            f"% at {self.get_general_settings().get('temperature', 298.15)} K",
        ]

        # minimal xtb energy from single-point (and mRRHO)
        gxtb = None
        gxtbmin = None
        if (
            any(type(p) is Prescreening for p in self._ensemble.results)
            and not self.get_general_settings()["gas-phase"]
        ):
            # Get the most recent prescreening part
            using_part = [p for p in self._ensemble.results if type(p) is Prescreening][
                -1
            ]

            gxtb = {
                conf.name: using_part.data["results"][conf.name]["xtb_gsolv"][
                    "energy_xtb_solv"
                ]
                for conf in self._ensemble.conformers
            }
        elif all(conf.xtb_energy is not None for conf in self._ensemble.conformers):
            gxtb = {conf.name: conf.xtb_energy for conf in self._ensemble.conformers}

        if self.get_general_settings()["evaluate_rrho"] and gxtb is not None:
            for conf in self._ensemble.conformers:
                gxtb[conf.name] += self.data["results"][conf.name]["xtb_rrho"]["gibbs"][
                    self.get_general_settings()["temperature"]
                ]
            gxtbmin = min(gxtb.values())

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            self.data["results"][conf.name]["gtot"]
            for conf in self._ensemble.conformers
        )

        # collect all dft single point energies
        dft_energies = (
            {
                conf.name: self.data["results"][conf.name]["sp"]["energy"]
                for conf in self._ensemble.conformers
            }
            if not all(
                "gsolv" in self.data["results"][conf.name]
                for conf in self._ensemble.conformers
            )
            else {
                conf.name: self.data["results"][conf.name]["gsolv"]["energy_gas"]
                for conf in self._ensemble.conformers
            }
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "G (xTB)": lambda conf: (
                f"{gxtb[conf.name]:.6f}" if gxtb is not None else "---"
            ),
            "ΔG (xTB)": lambda conf: (
                f"{(gxtb[conf.name] - gxtbmin) * AU2KCAL:.2f}"
                if gxtb is not None and gxtbmin is not None
                else "---"
            ),
            "E (DFT)": lambda conf: f"{dft_energies[conf.name]:.6f}",
            "ΔGsolv": lambda conf: (
                f"{self._gsolv(conf) - dft_energies[conf.name]:.6f}"
                if "gsolv" in self.data["results"][conf.name]
                else "---"
            ),
            "GmRRHO": lambda conf: (
                f"{self.data['results'][conf.name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
                if self.get_general_settings()["evaluate_rrho"]
                else "---"
            ),
            "Gtot": lambda conf: f"{self.data['results'][conf.name]['gtot']:.6f}",
            "ΔGtot": lambda conf: f"{(self.data['results'][conf.name]['gtot'] - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{self.data['results'][conf.name]['bmw'] * 100:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self._ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # list the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):\n"
        )
        lines.append(
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n"
        )

        # calculate averaged free enthalpy
        avG = sum(
            self.data["results"][conf.name]["bmw"]
            * self.data["results"][conf.name]["gtot"]
            for conf in self._ensemble.conformers
        )

        # calculate averaged free energy
        avE = (
            sum(
                self.data["results"][conf.name]["bmw"]
                * self.data["results"][conf.name]["sp"]["energy"]
                for conf in self._ensemble.conformers
            )
            if all(
                "sp" in self.data["results"][conf.name]
                for conf in self._ensemble.conformers
            )
            else sum(
                self.data["results"][conf.name]["bmw"]
                * self.data["results"][conf.name]["gsolv"]["energy_gas"]
                for conf in self._ensemble.conformers
            )
        )

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part1==\n"
        )
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # append lines to already existing file
        filename = f"{self._part_nos[self.name]}_{self.name.upper()}.out"
        with open(os.path.join(os.getcwd(), filename), "a", newline=None) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results to a json file
        self._write_json()


Factory.register_builder("screening", Screening)
