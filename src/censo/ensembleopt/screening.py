"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv using higher level solvation models and include the mRRHO contributions.
"""

from math import exp
from statistics import stdev
from tabulate import tabulate
from pathlib import Path
from typing import Any
import json
from collections.abc import Callable

from ..molecules import MoleculeData
from ..logging import setup_logger
from ..parallel import execute
from ..params import AU2KCAL, PLENGTH, NCORES, OMPMIN, GridLevel
from ..utilities import h1, printf, Factory, timeit, DataDump
from ..config import PartsConfig
from ..config.parts import ScreeningConfig
from ..config.job_config import RRHOJobConfig, SPJobConfig
from ..ensembledata import EnsembleData
from ..qm import QmProc

logger = setup_logger(__name__)


@timeit
def screening(
    ensemble: EnsembleData,
    config: PartsConfig,
    ncores: int = NCORES or OMPMIN,
    omp: int = OMPMIN,
    cut: bool = True,
):
    """
    Advanced screening of the ensemble by doing single-point calculations on the input geometries,
    but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

    Basically consists of two parts:
        - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
        - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
    """

    # Setup processor and target
    proc: QmProc = Factory[QmProc].create(config.screening.prog, "1_SCREENING")

    if not config.general.gas_phase and not config.screening.gsolv_included:
        # Calculate Gsolv using qm
        job_config = SPJobConfig(
            copy_mo=config.general.copy_mo,
            func=config.screening.func,
            basis=config.screening.basis,
            grid=GridLevel.MEDIUM,
            template=config.screening.template,
            gas_phase=False,
            solvent=config.general.solvent,
            sm=config.screening.sm,
        )
        results, _ = execute(
            ensemble.conformers,
            proc.gsolv,
            job_config,
            config.screening.prog,
            ncores,
            omp,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        for conf in ensemble:
            conf.gsolv = results[conf.name].gsolv
            conf.energy = results[conf.name].energy_gas
    else:
        # Run single-point calculation with solvation
        job_config = SPJobConfig(
            copy_mo=config.general.copy_mo,
            func=config.screening.func,
            basis=config.screening.basis,
            grid=GridLevel.MEDIUM,
            template=config.screening.template,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            sm=config.screening.sm,
        )
        results, _ = execute(
            ensemble.conformers,
            proc.sp,
            job_config,
            config.screening.prog,
            ncores,
            omp,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        for conf in ensemble:
            conf.energy = results[conf.name].energy

    if config.general.evaluate_rrho:
        # Run mRRHO calculation
        job_config = RRHOJobConfig(
            gfnv=config.screening.gfnv,
            **config.general.model_dump(),  # This will just let the constructor pick the key/value pairs it needs
        )
        results, _ = execute(
            ensemble.conformers,
            proc.xtb_rrho,
            job_config,
            "xtb",
            ncores,
            omp,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        for conf in ensemble:
            conf.grrho = results[conf.name].energy

    if cut:
        # Threshold in kcal/mol
        threshold = config.screening.threshold
        if len(ensemble.conformers) > 1:
            # calculate fuzzyness of threshold (adds 1 kcal/mol at max to the threshold)
            fuzzy = 1 - exp(
                -5 * (AU2KCAL * stdev([conf.grrho for conf in ensemble])) ** 2
            )
            threshold += fuzzy
            printf(f"Updated fuzzy threshold: {threshold:.2f} kcal/mol.")

        threshold = min(conf.gtot for conf in ensemble) + threshold * AU2KCAL

        ensemble.remove_conformers(cond=lambda conf: conf.gtot > threshold)

    ensemble.set_populations(config.general.temperature)

    # Print/write out results
    _write_results(ensemble, config)

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


def _write_results(ensemble: EnsembleData, config: PartsConfig) -> None:
    """ """
    printf(h1(f"SCREENING SINGLE-POINT (+ mRRHO) RESULTS"))

    # column headers
    headers = [
        "CONF#",
        # "G (xTB)",
        # "ΔG (xTB)",
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
        # "[Eh]",
        # "[kcal/mol]",
        "[Eh]",
        "[kcal/mol]",
        "[Eh]",
        "[Eh]",
        "[kcal/mol]",
        f"% at {config.general.temperature} K",
    ]

    # variables for printmap
    # minimal xtb single-point energy
    # if all(
    #     "xtb_gsolv" in self.data["results"][conf.name]
    #     for conf in self._ensemble.conformers
    # ):
    #     xtbmin = min(
    #         self.data["results"][conf.name]["xtb_gsolv"]["energy_xtb_gas"]
    #         for conf in self._ensemble.conformers
    #     )

    gtotmin = min(conf.gtot for conf in ensemble)

    printmap = {
        "CONF#": lambda conf: conf.name,
        # "G (xTB)": lambda conf: (
        #     f"{gxtb[conf.name]:.6f}" if gxtb is not None else "---"
        # ),
        # "ΔG (xTB)": lambda conf: (
        #     f"{(gxtb[conf.name] - gxtbmin) * AU2KCAL:.2f}"
        #     if gxtb is not None and gxtbmin is not None
        #     else "---"
        # ),
        "E (DFT)": lambda conf: f"{conf.energy:.6f}",
        "ΔGsolv": lambda conf: (f"{conf.gsolv * AU2KCAL:.6f}"),
        "GmRRHO": lambda conf: (f"{conf.grrho:.6f}"),
        "Gtot": lambda conf: f"{conf.gtot:.6f}",
        "ΔGtot": lambda conf: f"{(conf.gtot - gtotmin) * AU2KCAL:.2f}",
        "Boltzmann weight": lambda conf: f"{conf.bmw * 100:.2f}",
    }

    rows = [[printmap[header](conf) for header in headers] for conf in ensemble]

    for i in range(len(headers)):
        headers[i] += "\n" + units[i]

    table = tabulate(
        rows,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )
    print(table, flush=True)

    # list the averaged free enthalpy of the ensemble
    lines = []
    lines.append(
        "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):"
    )
    lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}")

    # calculate averaged free enthalpy
    avG = sum([conf.bmw * conf.gtot for conf in ensemble])

    # calculate averaged free energy (?)
    avE = sum([conf.bmw * conf.energy for conf in ensemble])

    # append the lines for the free energy/enthalpy
    lines.append(
        f"{config.general.temperature:^15} {avE:>14.7f}  {avG:>14.7f}     <<==part1=="
    )
    lines.append("".ljust(int(PLENGTH), "-"))

    # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")

    # Print everything
    for line in lines:
        printf(line)

    # write everything to a file
    filepath = Path("1_SCREENING.out")
    filepath.write_text(table + "\n".join(lines))

    # Additionally, write results in json format
    Path("1_SCREENING.json").write_text(
        json.dumps(jsonify(ensemble, config.screening), indent=4)
    )

    ensemble.dump_xyz(Path("1_SCREENING.xyz"))


# TODO: generalize this
def jsonify(
    ensemble: EnsembleData,
    config: ScreeningConfig,
    fields: Callable[[MoleculeData], dict[str, Any]] | None = None,
):
    per_conf: Callable[[MoleculeData], dict[str, dict[str, float]]] = fields or (
        lambda conf: {
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
                "gtot": conf.gtot,
            }
        }
    )

    dump = DataDump(part_name="screening")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump()
