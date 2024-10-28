import os
from functools import reduce

from ..logging import setup_logger
from ..parallel import execute
from ..params import AU2KCAL, GFNOPTIONS, PLENGTH, PROGS, SOLV_MODS
from ..utilities import format_data, h1, print, DfaHelper
from .prescreening import Prescreening
from .screening import Screening
from .optimization import Optimization

logger = setup_logger(__name__)


class Refinement(Screening):
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

    def _optimize(self, cut: bool = True) -> None:
        """
        Similar to Screening.optimize, however here we use a Boltzmann population cutoff instead of kcal cutoff.
        """
        Prescreening._optimize(self, cut=False)

        if self.get_general_settings()["evaluate_rrho"]:
            # Check if evaluate_rrho, then check if optimization was run and use that value, otherwise do xtb_rrho
            if not any(
                isinstance(part, Optimization) for part in self.ensemble.results
            ):
                jobtype = ["xtb_rrho"]
                prepinfo = self._setup_prepinfo(jobtype)

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
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Remove failed conformers
                self.ensemble.remove_conformers(failed)

                # Update results
                self.results.update(results)

                for conf in self.ensemble.conformers:
                    # calculate new gtot including RRHO contribution
                    self.results[conf.name]["gtot"] = self._grrho(conf)
            else:
                # Use values from (the most recent) optimization rrho
                part = next(
                    p
                    for p in self.ensemble.results[::-1]
                    if isinstance(p, Optimization)
                )
                for conf in self.ensemble.conformers:
                    self.results[conf.name]["xtb_rrho"] = part.results[conf.name][
                        "xtb_rrho"
                    ]
                    self.results[conf.name]["gtot"] = self._grrho(conf)

        # sort conformers list
        self.ensemble.conformers.sort(key=lambda conf: self.results[conf.name]["gtot"])

        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.results.update(self._calc_boltzmannweights())

        if cut:
            # Get Boltzmann population threshold from settings
            threshold = self.get_settings()["threshold"]

            # Update ensemble using Boltzman population threshold
            filtered = []
            bmws = 0
            for conf in self.ensemble.conformers:
                if bmws < threshold:
                    bmws += self.results[conf.name]["bmw"]
                else:
                    # Start filtering out conformers as soon as bmw sum crosses over threshold
                    filtered.append(conf.name)

            # Remove conformers
            self.ensemble.remove_conformers(filtered)
            for confname in filtered:
                print(f"No longer considering {confname}.")

            # Recalculate boltzmann weights after cutting down the ensemble
            self.results.update(self._calc_boltzmannweights())

        # second 'write_results' for the updated sorting with RRHO contributions
        self._write_results2()

    def _write_results2(self) -> None:
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
            self.results[conf.name]["gtot"] for conf in self.ensemble.conformers
        )

        # collect all dft single point energies
        dft_energies = (
            {
                id(conf): self.results[conf.name]["sp"]["energy"]
                for conf in self.ensemble.conformers
            }
            if not all(
                "gsolv" in self.results[conf.name] for conf in self.ensemble.conformers
            )
            else {
                id(conf): self.results[conf.name]["gsolv"]["energy_gas"]
                for conf in self.ensemble.conformers
            }
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT)": lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv": lambda conf: (
                f"{self._gsolv(conf) - dft_energies[id(conf)]:.6f}"
                if "gsolv" in self.results[conf.name]
                else "---"
            ),
            "GmRRHO": lambda conf: (
                f"{self.results[conf.name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
                if self.get_general_settings()["evaluate_rrho"]
                else "---"
            ),
            "Gtot": lambda conf: f"{self.results[conf.name]['gtot']:.6f}",
            "ΔGtot": lambda conf: f"{(self.results[conf.name]['gtot'] - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{self.results[conf.name]['bmw'] * 100:.2f}",
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
                self.results[conf.name]["bmw"] * self.results[conf.name]["gtot"]
                for conf in self.ensemble.conformers
            ]
        )

        # calculate averaged free energy
        avE = (
            sum(
                self.results[conf.name]["bmw"] * self.results[conf.name]["sp"]["energy"]
                for conf in self.ensemble.conformers
            )
            if all("sp" in self.results[conf.name] for conf in self.ensemble.conformers)
            else sum(
                self.results[conf.name]["bmw"]
                * self.results[conf.name]["gsolv"]["energy_gas"]
                for conf in self.ensemble.conformers
            )
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
        filename = f"{self._part_nos[self._name]}_{self._name.upper()}.out"
        logger.debug(f"Writing to {os.path.join(self.ensemble.workdir, filename)}.")
        with open(
            os.path.join(self.ensemble.workdir, filename), "a", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results to a json file
        self._write_json()
