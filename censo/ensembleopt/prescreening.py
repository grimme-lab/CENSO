import os
from typing import List
import csv
from functools import reduce

from censo.params import PLENGTH, AU2KCAL
from censo.utilities import (
    print,
    timeit,
    format_data,
    DfaHelper,
)
from censo.part import CensoPart
from censo.core import CensoCore
from censo.parallel import ProcessHandler
from censo.datastructure import MoleculeData
from censo.params import (
    SOLV_MODS,
    GSOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
)

"""

"""


class Prescreening(CensoPart):
    alt_name = "part0"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "threshold": {
            "default": 4.0,
            "range": [
                1.0,
                10.0
            ]
        },
        "func": {
            "default": "pbe-d4",
            "options": DfaHelper.find_func("prescreening")
        },
        "basis": {
            "default": "def2-SV(P)",
            "options": BASIS_SETS
        },
        "prog": {
            "default": "orca",
            "options": PROGS
        },
        "gfnv": {
            "default": "gfn2",
            "options": GFNOPTIONS
        },
        "grid": {
            "default": "low",
            "options": GRIDOPTIONS
        },
        "run": {
            "default": True
        },
        "gcp": {
            "default": True
        }
    }

    _settings = {}

    def __init__(self, core: CensoCore):
        super().__init__(core)

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        first screening of the ensemble by doing single-point calculation on the input geometries,
        using a (cheap) DFT method. if the ensemble ensembleopt is not taking place in the gas-phase,
        the gsolv contribution is calculated using xtb.

        the list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).
        """
        # print instructions
        self.print_info()

        # initialize process handler with conformer geometries
        handler = ProcessHandler(self._instructions, [conf.geom for conf in self.core.conformers])

        # set jobtype to pass to handler
        if not self._instructions["gas-phase"]:
            if self._instructions.get("implicit", False):
                jobtype = ["sp", "gsolv"]
            else:
                jobtype = ["sp", "xtb_gsolv"]
        else:
            jobtype = ["sp"]

        # compute results
        # for structure of results from handler.execute look there
        results = handler.execute(jobtype, self.dir)

        # update results for each conformer
        for conf in self.core.conformers:
            # if 'copy_mo' is enabled, get the mo_path from the results and store it in the respective GeometryData object
            if self._instructions["copy_mo"]:
                conf.geom.mo_path = results[id(conf)]["sp"]["mo_path"]

            # store results of single jobs for each conformer
            conf.results.setdefault(self._name, {}).update(results[id(conf)])

            # calculate free enthalpy values for every conformer
            conf.results[self._name]["gtot"] = self.gtot(conf)

        # sort conformers list with prescreening key (gtot)
        self.core.conformers.sort(
            key=lambda conf: conf.results[self._name]["gtot"],
        )

        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.core.calc_boltzmannweights(
            self._instructions.get("temperature", 298.15),
            self._name
        )

        # write results (analogous to deprecated print)
        self.write_results()

        # update conformers with threshold
        threshold = self._instructions["threshold"] / AU2KCAL

        # pick the free enthalpy of the lowest conformer
        limit = min([conf.results[self._name]["gtot"] for conf in self.core.conformers])

        # filter out all conformers above threshold
        # so that 'filtered' contains all conformers that should not be considered any further
        filtered = [
            conf for conf in filter(
                lambda x: self.gtot(x) - limit > threshold,
                self.core.conformers
            )
        ]

        # update the conformer list in core (remove conf if below threshold)
        self.core.update_conformers(filtered)

        # TODO - print out which conformers are no longer considered

        # dump ensemble
        self.core.dump_ensemble(self._name)

        # DONE

    def gtot(self, conf: MoleculeData) -> float:
        """
        prescreening key for conformer sorting
        calculates Gtot = E (DFT) + Gsolv (xtb) or Gsolv (DFT) for a given conformer
        """

        # Gtot = E (DFT) + Gsolv (xtb) or Gsolv (DFT)
        gtot: float = conf.results[self._name]["sp"]["energy"]

        # Gsolv is only calculated if prescreening is not in the gas-phase
        if "xtb_gsolv" in conf.results[self._name].keys():
            gtot += conf.results[self._name]["xtb_gsolv"]["gsolv"]
        elif "gsolv" in conf.results[self._name].keys():
            gtot += conf.results[self._name]["gsolv"]["gsolv"]

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
            f"\% at {self._instructions.get('temperature', 298.15)} K",
        ]

        # variables for printmap
        # minimal xtb single-point energy
        # TODO - where do prescreening and screening get xtb single-point from?
        # FIXME FIXME FIXME
        xtbmin = min(
            conf.results[self._name]['xtb_gsolv']['energy_xtb_gas']
            for conf in self.core.conformers
        )

        # minimal dft single-point energy
        dftmin = min(
            conf.results[self._name]['sp']['energy']
            for conf in self.core.conformers
        )

        # minimal solvation free enthalpy
        if self._instructions["gas-phase"]:
            gsolvmin = 0.0
        else:
            # NOTE: there might still be an error if a (xtb_)gsolv calculation failed for a conformer, therefore this should be handled before this step
            if all("xtb_gsolv" in conf.results[self._name].keys() for conf in
                   self.core.conformers):
                gsolvmin = min(
                    conf.results[self._name]['xtb_gsolv']['gsolv']
                    for conf in self.core.conformers
                )
            elif all("gsolv" in conf.results[self._name].keys() for conf in self.core.conformers):
                gsolvmin = min(
                    conf.results[self._name]['gsolv']['gsolv']
                    for conf in self.core.conformers
                )
            else:
                raise RuntimeError(
                    "The calculations should have used implicit or additive solvation for all conformers, but it is missing for at least some conformers.")

        # minimal total free enthalpy
        gtotmin = min(self.gtot(conf) for conf in self.core.conformers)

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda
                conf: f"{conf.results[self._name]['xtb_gsolv']['energy_xtb_gas']:.6f}",
            "ΔE (xTB)": lambda
                conf: f"{(conf.results[self._name]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}",
            "E (DFT)": lambda conf: f"{conf.results[self._name]['sp']['energy']:.6f}",
            "ΔGsolv (xTB)": lambda conf:
            f"{conf.results[self._name]['xtb_gsolv']['gsolv']:.6f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            else "---",
            "Gtot": lambda conf: f"{self.gtot(conf):.6f}",
            "ΔE (DFT)": lambda
                conf: f"{(conf.results[self._name]['sp']['energy'] - dftmin) * AU2KCAL:.2f}",
            "δΔGsolv": lambda conf:
            f"{(conf.results[self._name]['xtb_gsolv']['gsolv'] - gsolvmin) * AU2KCAL:.2f}"
            if "xtb_gsolv" in conf.results[self._name].keys()
            else "---",
            "ΔGtot": lambda conf: f"{(self.gtot(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # list the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):\n"
        )
        lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n")
        print("".ljust(int(PLENGTH), "-") + "\n")

        # calculate averaged free enthalpy
        avG = sum([
            conf.results[self._name]["bmw"]
            * conf.results[self._name]["gtot"]
            for conf in self.core.conformers
        ])

        # calculate averaged free energy
        avE = sum([
            conf.results[self._name]["bmw"]
            * conf.results[self._name]["sp"]["energy"]
            for conf in self.core.conformers
        ])

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self._instructions.get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0==\n")
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")

        # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")

        # write everything to a file
        with open(os.path.join(self.core.workdir, f"{self._name}.out"), "w",
                  newline=None) as outfile:
            outfile.writelines(lines)

        # additionally, write data in csv format
        with open(os.path.join(self.core.workdir, f"{self._name}.csv"), "w",
                  newline=None) as outfile:
            writer = csv.DictWriter(outfile, headers, delimiter=" ")
            writer.writeheader()
            rows = [{header: printmap[header](conf) for header in headers} for conf in self.core.conformers]
            writer.writerows(rows)
