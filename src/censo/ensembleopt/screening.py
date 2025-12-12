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
from dask.distributed import Client

from ..molecules import MoleculeData, Contributions
from ..logging import setup_logger
from ..parallel import execute
from ..params import AU2KCAL, PLENGTH, GridLevel, Prog
from ..utilities import h1, h2, printf, Factory, DataDump
from ..config import PartsConfig
from ..config.parts import ScreeningConfig
from ..config.job_config import RRHOJobConfig, SPJobConfig
from ..ensemble import EnsembleData
from ..processing import QmProc, XtbProc

logger = setup_logger(__name__)


def screening(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    cut: bool = True,
) -> dict[str, Any]:
    """
    Advanced screening of the ensemble by doing single-point calculations on the input geometries,
    but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

    Basically consists of two parts:
        - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
        - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time

    :param ensemble: EnsembleData object containing the conformers.
    :type ensemble: EnsembleData
    :param config: PartsConfig object with configuration settings.
    :type config: PartsConfig
    :param client: dask.distributed.Client for parallel execution.
    :type client: Client
    :param cut: Whether to apply cutting conditions.
    :type cut: bool
    :return: JSON-serializable dictionary of screening results.
    :rtype: dict
    """
    printf(h2("SCREENING"))

    config = PartsConfig.model_validate(config, context={"check": ["screening"]})

    # Setup processor and target
    proc = Factory[QmProc].create(config.screening.prog, "1_SCREENING")

    contributions_dict = {conf.name: Contributions() for conf in ensemble}
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
            # multitemp=config.general.multitemp,
            # trange=config.general.trange,
            temperature=config.general.temperature,
            paths=config.paths,
        )
        results = execute(
            ensemble.conformers,
            proc.gsolv,
            job_config,
            config.screening.prog,
            "screening",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )
        if config.general.ignore_failed:
            ensemble.remove_conformers(lambda conf: conf.name not in results)

        for conf in ensemble:
            contributions_dict[conf.name].gsolv = results[conf.name].gsolv
            contributions_dict[conf.name].energy = results[conf.name].energy_gas
    else:
        # Run sp calculation
        sp_job_config = SPJobConfig(
            func=config.screening.func,
            basis=config.screening.basis,
            grid=GridLevel.MEDIUM,
            template=config.screening.template,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            sm=config.screening.sm,
            paths=config.paths,
            copy_mo=config.general.copy_mo,
        )
        sp_results = execute(
            ensemble.conformers,
            proc.sp,
            sp_job_config,
            config.screening.prog,
            "screening",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )
        if config.general.ignore_failed:
            ensemble.remove_conformers(lambda conf: conf.name not in sp_results)

        for conf in ensemble:
            contributions_dict[conf.name].energy = sp_results[conf.name].energy

    if config.general.evaluate_rrho:
        # Run mRRHO calculation
        proc_xtb = Factory[XtbProc].create(Prog.XTB, "1_SCREENING")
        rrho_job_config = RRHOJobConfig(
            gfnv=config.screening.gfnv,
            paths=config.paths,
            **config.general.model_dump(),  # This will just let the constructor pick the key/value pairs it needs
        )
        rrho_results = execute(
            ensemble.conformers,
            proc_xtb.xtb_rrho,
            rrho_job_config,
            "xtb",
            "screening",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )
        if config.general.ignore_failed:
            ensemble.remove_conformers(lambda conf: conf.name not in rrho_results)

        for conf in ensemble:
            contributions_dict[conf.name].grrho = rrho_results[conf.name].energy

    # Update molecules
    ensemble.update_contributions(contributions_dict)

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

        threshold = min(conf.gtot for conf in ensemble) + threshold / AU2KCAL

        ensemble.remove_conformers(cond=lambda conf: conf.gtot > threshold)

    # Print/write out results
    results_dict = jsonify(ensemble, config.screening)
    _write_results(ensemble, config, results_dict)

    return results_dict


def _write_results(
    ensemble: EnsembleData, config: PartsConfig, results_dict: dict[str, Any]
) -> None:
    """ """
    printf(h1("SCREENING SINGLE-POINT (+ mRRHO) RESULTS"))

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

    boltzmann_populations = ensemble.get_populations(config.general.temperature)

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
        "Boltzmann weight": lambda conf: f"{boltzmann_populations[conf.name] * 100:.2f}",
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
    avG = sum([boltzmann_populations[conf.name] * conf.gtot for conf in ensemble])

    # calculate averaged free energy (?)
    avE = sum([boltzmann_populations[conf.name] * conf.energy for conf in ensemble])

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
    Path("1_SCREENING.json").write_text(json.dumps(results_dict, indent=4))

    ensemble.dump_xyz(Path("1_SCREENING.xyz"))


# TODO: generalize this
def jsonify(
    ensemble: EnsembleData,
    config: ScreeningConfig,
    fields: Callable[[MoleculeData], dict[str, Any]] | None = None,
):
    """
    Convert ensemble data to JSON format for screening results.

    :param ensemble: EnsembleData object.
    :type ensemble: EnsembleData
    :param config: ScreeningConfig object.
    :type config: ScreeningConfig
    :param fields: Optional callable to customize fields.
    :type fields: Callable[[MoleculeData], dict[str, Any]] | None
    :return: JSON-serializable dictionary.
    :rtype: dict[str, Any]
    """
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

    return dump.model_dump(mode="json")
