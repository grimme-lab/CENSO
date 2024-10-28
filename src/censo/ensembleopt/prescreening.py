import os

from ..datastructure import MoleculeData
from ..ensembledata import EnsembleData
from ..logging import setup_logger
from ..parallel import execute
from ..params import Params
from ..utilities import format_data, h1, print, DfaHelper
from .optimizer import EnsembleOptimizer

logger = setup_logger(__name__)


class Prescreening(EnsembleOptimizer):
    _grid = "low"

    _options = {
        "threshold": {"default": 4.0},
        "func": {
            "default": "pbe-d4",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in Params.PROGS},
        },
        "basis": {"default": "def2-SV(P)"},
        "prog": {"default": "tm", "options": Params.PROGS},
        "gfnv": {"default": "gfn2", "options": Params.GFNOPTIONS},
        "run": {"default": True},
        "template": {"default": False},
    }

    _settings = {}

    def _optimize(self, cut: bool = True) -> None:
        """
        first screening of the ensemble by doing single-point calculation on the input geometries,
        using a (cheap) DFT method. if the ensemble ensembleopt is not taking place in the gas-phase,
        the gsolv contribution is calculated using xtb.

        The list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).
        """
        # set jobtype to pass to handler
        # TODO - it is not very nice to partially handle 'Screening' settings here
        if self.get_general_settings()["gas-phase"]:
            jobtype = ["sp"]
        elif self.get_settings().get("implicit", False):
            if self.get_settings().get("sm", None) in [
                "cosmors",
                "cosmors-fine",
            ]:
                # If cosmors is used as solvent model the gsolv calculation needs to be done explicitely
                logger.warning(
                    "COSMORS detected as solvation model, this requires explicit calculation of ΔGsolv."
                )
                jobtype = ["gsolv"]
            else:
                # 'implicit' is a special option of Screening that makes CENSO skip the explicit computation of Gsolv
                # Gsolv will still be included in the DFT energy though
                jobtype = ["sp"]
        elif not self.get_settings().get("implicit", False):
            # Only for prescreening the solvation should be calculated with xtb
            if self._name == "prescreening":
                jobtype = ["xtb_gsolv"]

                # Compile all information required for the preparation of input files in parallel execution step
                prepinfo = self._setup_prepinfo(jobtype)

                # compute results
                # for structure of results from handler.execute look there
                results, failed = execute(
                    self.ensemble.conformers,
                    self.dir,
                    self.get_settings()["prog"],
                    prepinfo,
                    jobtype,
                    copy_mo=self.get_general_settings()["copy_mo"],
                    balance=self.get_general_settings()["balance"],
                    omp=self.get_general_settings()["omp"],
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Remove failed conformers
                self.ensemble.remove_conformers(failed)

                # Update results
                self.results.update(results)

                jobtype = ["sp"]
            else:
                jobtype = ["gsolv"]

        # Compile all information required for the preparation of input files in parallel execution step
        prepinfo = self._setup_prepinfo(jobtype)

        # compute results
        # for structure of results from handler.execute look there
        results, failed = execute(
            self.ensemble.conformers,
            self.dir,
            self.get_settings()["prog"],
            prepinfo,
            jobtype,
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            omp=self.get_general_settings()["omp"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

        # Update results
        self.results.update(results)

        # update results for each conformer
        for conf in self.ensemble.conformers:
            # calculate free enthalpy
            self.results[conf.name]["gtot"] = self._gsolv(conf)

        # sort conformers list with prescreening key (gtot)
        self.ensemble.conformers.sort(
            key=lambda conf: self.results[conf.name]["gtot"],
        )

        # calculate boltzmann weights from gtot values calculated here
        self.results.update(self._calc_boltzmannweights())

        self._write_results()

        if cut:
            print("\n")
            # update conformers with threshold
            threshold = self.get_settings()["threshold"] / Params.Params.AU2KCAL
            limit = min(
                self.results[conf.name]["gtot"] for conf in self.ensemble.conformers
            )
            filtered = list(
                filter(
                    lambda conf: self.results[conf.name]["gtot"] - limit > threshold,
                    self.ensemble.conformers,
                )
            )

            # update the conformer list in ensemble (remove confs if below threshold)
            self.ensemble.remove_conformers([conf.name for conf in filtered])
            for conf in filtered:
                print(f"No longer considering {conf.name}.")

    def _gsolv(self, conf: MoleculeData) -> float:
        """
        Prescreening key for conformer sorting
        Calculates Gtot = E (DFT) + Gsolv (xtb) for a given conformer
        """

        # Gtot = E (DFT) + Gsolv (xtb)
        if not self.get_general_settings()["gas-phase"]:
            gtot = (
                self.results[conf.name]["sp"]["energy"]
                + self.results[conf.name]["xtb_gsolv"]["gsolv"]
            )
        else:
            gtot = self.results[conf.name]["sp"]["energy"]

        return gtot

    def _write_results(self) -> None:
        """
        writes:
            E (xtb),
            δE (xtb),
            G_solv (xtb),
            δG_solv,

            E(DFT),
            δE(DFT),

            E(DFT) + G_solv,
            δ(E(DFT) + G_solv)

        also writes data in easily digestible format
        """
        print(h1(f"{self._name.upper()} SINGLE-POINT RESULTS"))

        # column headers
        headers = [
            "CONF#",
            "E (xTB)",
            "ΔE (xTB)",
            "E (DFT)",
            "ΔE (DFT)",
            "ΔGsolv (xTB)",
            # "δΔGsolv",
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
            "[kcal/mol]",
            "[kcal/mol]",
            # "[kcal/mol]",
            "[Eh]",
            "[kcal/mol]",
            f"% at {self.get_general_settings().get('temperature', 298.15)} K",
        ]

        # variables for printmap
        # minimal xtb single-point energy
        xtbmin = None
        if all(
            "xtb_gsolv" in self.results[conf.name].keys()
            for conf in self.ensemble.conformers
        ):
            xtbmin = min(
                self.results[conf.name]["xtb_gsolv"]["energy_xtb_gas"]
                for conf in self.ensemble.conformers
            )

        # minimal dft single-point energy
        dft_energies = (
            {
                id(conf): self.results[conf.name]["sp"]["energy"]
                for conf in self.ensemble.conformers
            }
            if not all(
                "gsolv" in self.results[conf.name].keys()
                for conf in self.ensemble.conformers
            )
            else {
                id(conf): self.results[conf.name]["gsolv"]["energy_gas"]
                for conf in self.ensemble.conformers
            }
        )

        dftmin = min(dft_energies.values())

        # minimal solvation free enthalpy
        if self.get_general_settings()["gas-phase"]:
            gsolvmin = 0.0
        else:
            # NOTE: there might still be an error if a (xtb_)gsolv calculation failed for a conformer, therefore this should be handled before this step
            if all(
                "xtb_gsolv" in self.results[conf.name].keys()
                for conf in self.ensemble.conformers
            ):
                gsolvmin = min(
                    self.results[conf.name]["xtb_gsolv"]["gsolv"]
                    for conf in self.ensemble.conformers
                )
            elif all(
                "gsolv" in self.results[conf.name].keys()
                for conf in self.ensemble.conformers
            ):
                gsolvmin = min(
                    self.results[conf.name]["gsolv"]["gsolv"]
                    for conf in self.ensemble.conformers
                )
            else:
                raise RuntimeError(
                    "The calculations should have used implicit or additive solvation for all conformers, "
                    "but it is missing for at least some conformers."
                )

        # minimal total free enthalpy
        gtotmin = min(self._gsolv(conf) for conf in self.ensemble.conformers)

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: (
                f"{self.results[conf.name]['xtb_gsolv']['energy_xtb_gas']:.6f}"
                if "xtb_gsolv" in self.results[conf.name].keys()
                else "---"
            ),
            "ΔE (xTB)": lambda conf: (
                f"{(self.results[conf.name]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * Params.AU2KCAL:.2f}"
                if "xtb_gsolv" in self.results[conf.name].keys()
                else "---"
            ),
            "E (DFT)": lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔE (DFT)": lambda conf: f"{(dft_energies[id(conf)] - dftmin) * Params.AU2KCAL:.2f}",
            "ΔGsolv (xTB)": lambda conf: (
                f"{self.results[conf.name]['xtb_gsolv']['gsolv'] * Params.AU2KCAL:.6f}"
                if "xtb_gsolv" in self.results[conf.name].keys()
                else "---"
            ),
            "Gtot": lambda conf: f"{self._gsolv(conf):.6f}",
            # "δΔGsolv": lambda conf: f"{(self.results[conf.name]['xtb_gsolv']['gsolv'] - gsolvmin) * Params.AU2KCAL:.2f}"
            # if "xtb_gsolv" in self.results[conf.name].keys()
            # else "---",
            "ΔGtot": lambda conf: f"{(self._gsolv(conf) - gtotmin) * Params.AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{self.results[conf.name]['bmw'] * 100:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
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
            [
                self.results[conf.name]["bmw"] * self.results[conf.name]["gtot"]
                for conf in self.ensemble.conformers
            ]
        )

        # calculate averaged free energy
        avE = sum(
            [
                self.results[conf.name]["bmw"] * self.results[conf.name]["sp"]["energy"]
                for conf in self.ensemble.conformers
            ]
        )

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0==\n"
        )
        lines.append("".ljust(int(Params.PLENGTH), "-"))

        # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # write everything to a file
        filename = f"{self._part_nos[self._name]}_{self._name.upper()}.out"
        logger.debug(f"Writing to {os.path.join(self.ensemble.workdir, filename)}.")
        with open(
            os.path.join(self.ensemble.workdir, filename), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write results in json format
        self._write_json()
