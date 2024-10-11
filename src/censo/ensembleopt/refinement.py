import os
from functools import reduce

from ..logging import setup_logger
from ..parallel import execute
from ..params import AU2KCAL, GFNOPTIONS, PLENGTH, PROGS, SOLV_MODS
from ..utilities import format_data, h1, print, DfaHelper
from .prescreening import Prescreening
from .screening import Screening

logger = setup_logger(__name__)


class Refinement(Screening):
    _part_no = "3"

    _grid = "high+"

    __solv_mods = {prog: SOLV_MODS[prog] for prog in PROGS}
    # __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "threshold": {"default": 0.95},
        "func": {
            "default": "wb97x-v",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in PROGS},
        },
        "basis": {"default": "def2-TZVP"},
        "prog": {"default": "tm", "options": PROGS},
        "sm": {"default": "cosmors", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
        "run": {"default": True},
        "implicit": {"default": False},
        "template": {"default": False},
    }

    _settings = {}

    def optimize(self, ncores: int, cut: bool = True) -> None:
        """
        Similar to Screening.optimize, however here we use a Boltzmann population cutoff instead of kcal cutoff.
        """
        Prescreening.optimize(self, ncores, cut=False)

        if self.get_general_settings()["evaluate_rrho"]:
            # Check if evaluate_rrho, then check if optimization was run and use that value, otherwise do xtb_rrho
            if not all(
                "optimization" in conf.results.keys()
                for conf in self.ensemble.conformers
            ):
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
                    maxcores=ncores,
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Remove failed conformers
                self.ensemble.remove_conformers(failed)

                for conf in self.ensemble.conformers:
                    # calculate new gtot including RRHO contribution
                    conf.results[self._name]["gtot"] = self.grrho(conf)
            else:
                # Use values from optimization rrho
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["xtb_rrho"] = conf.results["optimization"][
                        "xtb_rrho"
                    ]
                    conf.results[self._name]["gtot"] = self.grrho(conf)

        # sort conformers list
        self.ensemble.conformers.sort(key=lambda conf: conf.results[self._name]["gtot"])

        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        if cut:
            # Get Boltzmann population threshold from settings
            threshold = self.get_settings()["threshold"]

            # Update ensemble using Boltzman population threshold
            for confname in self.ensemble.update_conformers(
                lambda conf: conf.results[self._name]["bmw"], threshold, boltzmann=True
            ):
                print(f"No longer considering {confname}.")

        # Recalculate boltzmann weights after cutting down the ensemble
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        # second 'write_results' for the updated sorting with RRHO contributions
        self.write_results2()

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
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
            f"% at {self.get_general_settings().get('temperature', 298.15)} K",
        ]

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            conf.results[self._name]["gtot"] for conf in self.ensemble.conformers
        )

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
            "ΔGsolv": lambda conf: (
                f"{self.gtot(conf) - dft_energies[id(conf)]:.6f}"
                if not self.get_settings().get("implicit", False)
                else "---"
            ),
            "GmRRHO": lambda conf: (
                f"{conf.results[self._name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
                if self.get_general_settings()["evaluate_rrho"]
                else "---"
            ),
            "Gtot": lambda conf: f"{conf.results[self._name]['gtot']:.6f}",
            "ΔGtot": lambda conf: f"{(conf.results[self._name]['gtot'] - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # list the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble (high level single-points):\n"
        )
        lines.append(
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n"
        )

        # calculate averaged free enthalpy
        avG = sum(
            [
                conf.results[self._name]["bmw"] * conf.results[self._name]["gtot"]
                for conf in self.ensemble.conformers
            ]
        )

        # calculate averaged free energy
        avE = sum(
            [
                conf.results[self._name]["bmw"]
                * conf.results[self._name]["sp"]["energy"]
                for conf in self.ensemble.conformers
            ]
        )

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part3==\n"
        )
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # append lines to already existing file
        filename = f"{self._part_no}_{self._name.upper()}.out"
        logger.debug(f"Writing to {os.path.join(self.ensemble.workdir, filename)}.")
        with open(
            os.path.join(self.ensemble.workdir, filename), "a", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results to a json file
        self.write_json()
