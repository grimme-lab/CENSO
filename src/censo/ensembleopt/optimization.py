"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""

from copy import deepcopy
import os
import json
from pathlib import Path
from collections.abc import Callable
from typing import Any
from tabulate import tabulate

from ..ensembledata import EnsembleData
from ..molecules import MoleculeData
from ..qm import QmProc
from ..parallel import OptResult, execute
from ..params import AU2KCAL, PLENGTH, NCORES, OMPMIN, GridLevel
from ..config import PartsConfig
from ..config.parts import OptimizationConfig
from ..config.job_config import RRHOJobConfig, OptJobConfig, XTBOptJobConfig
from ..utilities import (
    printf,
    h1,
    Factory,
    timeit,
    DataDump,
)
from ..logging import setup_logger

logger = setup_logger(__name__)


@timeit
def optimization(
    ensemble: EnsembleData,
    config: PartsConfig,
    ncores: int = NCORES or OMPMIN,
    omp: int = OMPMIN,
    cut: bool = True,
):
    """
    Geometry optimization of the ensemble at DFT level (possibly with implicit solvation)

    Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
    by calculating the thermodynamics of the optimized geometries

    Alternatively just run the complete geometry optimization for every conformer with xtb as driver (decide with 'macrocycles')
    """
    # Setup processor
    proc: QmProc = Factory[QmProc].create(config.optimization.prog, "2_OPTIMIZATION")

    if config.optimization.macrocycles:
        _macrocycle_opt(proc, ensemble, config, ncores, omp, cut)
    else:
        _full_opt(proc, ensemble, config, ncores, omp)

    if config.general.evaluate_rrho:
        # Run mRRHO calculation
        job_config = RRHOJobConfig(
            gfnv=config.optimization.gfnv,
            **config.general.model_dump(),
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
        threshold = (
            min(conf.gtot for conf in ensemble)
            + config.optimization.threshold * AU2KCAL
        )
        ensemble.remove_conformers(lambda conf: conf.gtot > threshold)

    ensemble.set_populations(config.general.temperature)

    # Print/write out results
    _write_results(ensemble, config)


def _macrocycle_opt(
    proc: QmProc,
    ensemble: EnsembleData,
    config: PartsConfig,
    ncores: int,
    omp: int,
    cut: bool,
):
    """
    Geometry optimization using macrocycles, whereafter every macrocycle cutting conditions are checked (if cut == True).
    """
    # Set up a list of unconverged conformers (perhaps a second ensemble?)
    unconverged_ensemble = deepcopy(ensemble)

    # Set up target
    if config.optimization.xtb_opt:
        job_config = XTBOptJobConfig(
            grid=GridLevel.MEDIUM,
            copy_mo=config.general.copy_mo,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            **config.optimization.model_dump(),
        )
        target = proc.xtb_opt
    else:
        job_config = OptJobConfig(
            grid=GridLevel.MEDIUM,
            copy_mo=config.general.copy_mo,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            **config.optimization.model_dump(),
        )
        target = proc.opt

    # Iterate and cut only the unconverged conformers
    ncyc = 0
    rrho_done = False
    while (
        len(unconverged_ensemble.conformers) > 0 and ncyc < config.optimization.maxcyc
    ):
        results, _ = execute(
            unconverged_ensemble.conformers,
            target,
            job_config,
            config.optimization.prog,
            ncores,
            omp,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        for conf in unconverged_ensemble:
            # Update energies
            conf.energy = results[conf.name].energy

            # Update geometries
            conf.geom.xyz = results[conf.name].geom

        if (
            config.general.evaluate_rrho
            and ncyc + config.optimization.optcycles > 6
            and not rrho_done
        ):
            # Run mRRHO calculation
            job_config_rrho = RRHOJobConfig(
                gfnv=config.optimization.gfnv,
                **config.general.model_dump(),
            )
            results_rrho, _ = execute(
                unconverged_ensemble.conformers,
                proc.xtb_rrho,
                job_config_rrho,
                "xtb",
                ncores,
                omp,
                balance=config.general.balance,
                copy_mo=config.general.copy_mo,
            )

            for conf in ensemble:
                conf.grrho = results_rrho[conf.name].energy

            rrho_done = True

        # TODO: molbar topology check

        # Remove converged conformers
        unconverged_ensemble.remove_conformers(
            lambda conf: results[conf.name].converged
        )

        if len(unconverged_ensemble.conformers) != 0 and cut:
            gradthr = config.optimization.gradthr

            # Take the minimum Gtot for all the conformers (converged AND unconverged) (converged are considered removed from the unconverged_ensemble)
            threshold = (
                min(
                    [conf.gtot for conf in unconverged_ensemble]
                    + [conf.gtot for conf in unconverged_ensemble.rem]
                )
                + config.optimization.threshold
            )
            unconverged_ensemble.remove_conformers(
                lambda conf: conf.gtot > threshold
                and results[conf.name].grad_norm < gradthr
            )
        # if all confs are converged, don't cut again, instead let the cutting be done in the optimization function,
        # so we can use the updated mRRHO contributions
        ncyc += config.optimization.optcycles

        _print_update(unconverged_ensemble, results)

    # Finally sync ensembles
    for conf in unconverged_ensemble.rem:
        conf2 = next(_conf2 for _conf2 in ensemble if _conf2.name == conf.name)
        conf2.energy = conf.energy
        conf2.grrho = conf.grrho

    for conf in unconverged_ensemble:
        logger.warning(
            f"{conf.name} did not converge after {config.optimization.maxcyc} cycles."
        )

    ensemble.remove_conformers(lambda conf: conf in unconverged_ensemble)


def _full_opt(
    proc: QmProc,
    ensemble: EnsembleData,
    config: PartsConfig,
    ncores: int,
    omp: int,
):
    """
    Full geometry optimization of every conformer.
    """
    # Set up target
    if config.optimization.xtb_opt:
        job_config = XTBOptJobConfig(
            copy_mo=config.general.copy_mo,
            grid=GridLevel.MEDIUM,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            **config.optimization.model_dump(),
        )
        target = proc.xtb_opt
    else:
        job_config = OptJobConfig(
            copy_mo=config.general.copy_mo,
            grid=GridLevel.MEDIUM,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            **config.optimization.model_dump(),
        )
        target = proc.opt

    results, _ = execute(
        ensemble.conformers,
        target,
        job_config,
        config.optimization.prog,
        ncores,
        omp,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )

    for conf in ensemble:
        conf.energy = results[conf.name].energy
        conf.geom.xyz = results[conf.name].geom


def _write_results(ensemble: EnsembleData, config: PartsConfig) -> None:
    """ """
    printf(h1(f"OPTIMIZATION RESULTS"))

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
    lines: list[str] = []
    lines.append(
        "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (DFT optimized):"
    )
    lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}")

    # calculate averaged free enthalpy
    avG = sum([conf.bmw * conf.gtot for conf in ensemble])

    # calculate averaged free energy (?)
    avE = sum([conf.bmw * conf.energy for conf in ensemble])

    # append the lines for the free energy/enthalpy
    lines.append(
        f"{config.general.temperature:^15} {avE:>14.7f}  {avG:>14.7f}     <<==part2=="
    )
    lines.append("".ljust(int(PLENGTH), "-"))

    # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")

    # Print everything
    for line in lines:
        printf(line)

    # write everything to a file
    filepath = Path("2_OPTIMIZATION.out")
    filepath.write_text(table + "\n".join(lines))

    # Additionally, write results in json format
    Path("2_OPTIMIZATION.json").write_text(
        json.dumps(jsonify(ensemble, config.optimization), indent=4)
    )

    ensemble.dump_xyz(Path("2_OPTIMIZATION.xyz"))


# TODO: generalize this
def jsonify(
    ensemble: EnsembleData,
    config: OptimizationConfig,
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

    dump = DataDump(part_name="optimization")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump()


def _print_update(ensemble: EnsembleData, results: dict[str, OptResult]):
    pass
