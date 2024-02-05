import os

from .optimizer import EnsembleOptimizer
from ..ensembledata import EnsembleData
from ..datastructure import MoleculeData
from ..parallel import execute
from ..params import PLENGTH, AU2KCAL
from ..params import (
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
    SOLVENTS_DB,
)
from ..part import CensoPart
from ..utilities import (
    timeit,
    format_data,
    DfaHelper,
    setup_logger,
)

logger = setup_logger(__name__)


class Prescreening(EnsembleOptimizer):
    alt_name = "part0"

    _options = {
        "threshold": {"default": 4.0, "range": [1.0, 10.0]},
        "func": {"default": "pbe-d4", "options": []},
        "basis": {"default": "def2-SV(P)", "options": []},
        "prog": {"default": "orca", "options": PROGS},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
        "grid": {"default": "low", "options": GRIDOPTIONS},
        "run": {"default": True},
        "gcp": {"default": True},
        "template": {"default": False},
    }

    _settings = {}

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    def optimize(self, cut: bool = True) -> None:
        """
        first screening of the ensemble by doing single-point calculation on the input geometries,
        using a (cheap) DFT method. if the ensemble ensembleopt is not taking place in the gas-phase,
        the gsolv contribution is calculated using xtb.

        The list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).
        """
        # set jobtype to pass to handler
        # TODO - it is not very nice to partially handle 'Screening' settings here
        if self.get_general_settings()["gas-phase"] or self.get_settings().get("implicit", False):
            # 'implicit' is a special option of Screening that makes CENSO skip the explicit computation of Gsolv
            # Gsolv will still be included in the DFT energy though
            jobtype = ["sp"]
        elif not self.get_settings().get("implicit", False):
            jobtype = ["xtb_gsolv", "sp"]
        else:
            jobtype = ["gsolv"]

        # Compile all information required for the preparation of input files in parallel execution step
        prepinfo = self.setup_prepinfo(jobtype)

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
            maxcores=self.get_general_settings()["maxcores"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

        # update results for each conformer
        for conf in self.ensemble.conformers:
            # store results of single jobs for each conformer
            conf.results.setdefault(self._name, {}).update(results[id(conf)])

            # calculate free enthalpy values for every conformer
            conf.results[self._name]["gtot"] = self.gsolv(conf)

        # sort conformers list with prescreening key (gtot)
        self.ensemble.conformers.sort(
            key=lambda conf: conf.results[self._name]["gtot"],
        )

        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        self.write_results()

        if cut:
            # update conformers with threshold
            threshold = self.get_settings()["threshold"] / AU2KCAL

            # update the conformer list in ensemble (remove confs if below threshold)
            for confname in self.ensemble.update_conformers(self.gsolv, threshold):
                print(f"No longer considering {confname}.")

    def setup_prepinfo(self, jobtype: list[str]) -> dict[str, dict]:
        prepinfo = {jt: {} for jt in jobtype}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        prepinfo["sp"] = {
            "func_name": DfaHelper.get_name(
                self.get_settings()["func"], self.get_settings()["prog"]
            ),
            "func_type": DfaHelper.get_type(
                self.get_settings()["func"]),
            "disp": DfaHelper.get_disp(
                self.get_settings()["func"]),
            "basis": self.get_settings()["basis"],
            "grid": self.get_settings()["grid"],
            "template": self.get_settings()["template"],
            "gcp": self.get_settings()["gcp"],
        }

        # Add the solvent key if a solvent model exists in the part settings (this method is also used for Screening)
        # NOTE: 'sm' in key catches also cases like NMR (sm_s and sm_j)
        if any("sm" in key for key in self.get_settings().keys()):
            prepinfo["sp"]["sm"] = self.get_settings()["sm"]
            prepinfo["sp"]["solvent_key_prog"] = SOLVENTS_DB.get(
                self.get_general_settings()["solvent"])[self.get_settings()["sm"]][1]

        # TODO - this doesn't look very nice
        if "xtb_gsolv" in jobtype:
            # NOTE: [1] auto-selects replacement solvent (TODO - print warning!)
            prepinfo["xtb_sp"] = {
                "gfnv": self.get_settings()["gfnv"],
                "solvent_key_xtb": SOLVENTS_DB.get(self.get_general_settings()["solvent"])["xtb"][1],
            }

        if "xtb_rrho" in jobtype:
            prepinfo["xtb_rrho"] = {
                "gfnv": self.get_settings()["gfnv"],
                "solvent_key_xtb": SOLVENTS_DB.get(self.get_general_settings()["solvent"])["xtb"][1],
            }

        return prepinfo

    def gsolv(self, conf: MoleculeData) -> float:
        """
        Prescreening key for conformer sorting
        Calculates Gtot = E (DFT) + Gsolv (xtb) for a given conformer
        """

        # Gtot = E (DFT) + Gsolv (xtb)
        if (
            not self.get_general_settings()["gas-phase"]
        ):
            gtot = (
                conf.results[self._name]["sp"]["energy"]
                + conf.results[self._name]["xtb_gsolv"]["gsolv"]
            )
        else:
            gtot = conf.results[self._name]["sp"]["energy"]

        return gtot

    def write_results(self) -> None:
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

        # column headers
        headers = [
            "CONF#",
            "E (xTB)",
            "ΔE (xTB)",
            "E (DFT)",
            "ΔGsolv (xTB)",
            "Gtot",
            "ΔE (DFT)",
            "δΔGsolv",
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
            "[kcal/mol]",
            "[kcal/mol]",
            "[kcal/mol]",
            f"% at {self.get_general_settings().get('temperature', 298.15)} K",
        ]

        # variables for printmap
        # minimal xtb single-point energy
        if all(
            "xtb_gsolv" in conf.results[self._name].keys()
            for conf in self.ensemble.conformers
        ):
            xtbmin = min(
                conf.results[self._name]["xtb_gsolv"]["energy_xtb_gas"]
                for conf in self.ensemble.conformers
            )

        # minimal dft single-point energy
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

        dftmin = min(dft_energies.values())

        # minimal solvation free enthalpy
        if self.get_general_settings()["gas-phase"]:
            gsolvmin = 0.0
        else:
            # NOTE: there might still be an error if a (xtb_)gsolv calculation failed for a conformer, therefore this should be handled before this step
            if all(
                "xtb_gsolv" in conf.results[self._name].keys()
                for conf in self.ensemble.conformers
            ):
                gsolvmin = min(
                    conf.results[self._name]["xtb_gsolv"]["gsolv"]
                    for conf in self.ensemble.conformers
                )
            elif all(
                "gsolv" in conf.results[self._name].keys()
                for conf in self.ensemble.conformers
            ):
                gsolvmin = min(
                    conf.results[self._name]["gsolv"]["gsolv"]
                    for conf in self.ensemble.conformers
                )
            else:
                raise RuntimeError(
                    "The calculations should have used implicit or additive solvation for all conformers, "
                    "but it is missing for at least some conformers."
                )

        # minimal total free enthalpy
        gsolvmin = min(self.gsolv(conf) for conf in self.ensemble.conformers)

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: f"{conf.results[self._name]['xtb_gsolv']['energy_xtb_gas']:.6f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            else "---",
            "ΔE (xTB)": lambda conf: f"{(conf.results[self._name]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            else "---",
            "E (DFT)": lambda conf: f"{dft_energies[id(conf)]:.6f}",
            "ΔGsolv (xTB)": lambda conf: f"{conf.results[self._name]['xtb_gsolv']['gsolv']:.6f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            else "---",
            "Gtot": lambda conf: f"{self.gsolv(conf):.6f}",
            "ΔE (DFT)": lambda conf: f"{(dft_energies[id(conf)] - dftmin) * AU2KCAL:.2f}",
            "δΔGsolv": lambda conf: f"{(conf.results[self._name]['xtb_gsolv']['gsolv'] - gsolvmin) * AU2KCAL:.2f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            else "---",
            "ΔGtot": lambda conf: f"{(self.gsolv(conf) - gsolvmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
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
        print("".ljust(int(PLENGTH), "-") + "\n")

        # calculate averaged free enthalpy
        avG = sum(
            [
                conf.results[self._name]["bmw"] *
                conf.results[self._name]["gtot"]
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
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0==\n"
        )
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")

        # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # write everything to a file
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, f'{self._name}.out')}."
        )
        with open(
            os.path.join(self.ensemble.workdir, f"{self._name}.out"), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write results in json format
        self.write_json()
