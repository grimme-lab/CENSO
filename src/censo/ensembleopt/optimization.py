"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""

from copy import copy
import json
from pathlib import Path
from collections.abc import Callable
from typing import Any
from tabulate import tabulate
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..molecules import MoleculeData, Contributions
from ..processing import QmProc, XtbProc
from ..parallel import execute
from ..params import AU2KCAL, PLENGTH, GridLevel, Prog
from ..config import PartsConfig
from ..config.parts import OptimizationConfig
from ..config.job_config import RRHOJobConfig, OptJobConfig, XTBOptJobConfig
from ..utilities import (
    printf,
    h1,
    Factory,
    DataDump,
    h2,
)
from ..logging import setup_logger

logger = setup_logger(__name__)


def optimization(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    cut: bool = True,
) -> dict[str, Any]:
    """
    Geometry optimization of the ensemble at DFT level (possibly with implicit solvation)

    Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
    by calculating the thermodynamics of the optimized geometries

    Alternatively just run the complete geometry optimization for every conformer with xtb as driver (decide with 'macrocycles')

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param client: dask.distributed.Client for parallel execution.
    :param cut: Whether to apply cutting conditions during optimization.
    :return: JSON-serializable dictionary of optimization results.
    """
    printf(h2("OPTIMIZATION"))

    config = PartsConfig.model_validate(config, context={"check": ["optimization"]})

    # Setup processor
    proc = Factory[QmProc].create(config.optimization.prog, "2_OPTIMIZATION")

    if config.optimization.macrocycles:
        contributions_dict = _macrocycle_opt(proc, ensemble, config, client, cut)
    else:
        contributions_dict = _full_opt(proc, ensemble, config, client=client)

    printf("\n")

    # Potentially reevaluate rrho contributions
    if config.general.evaluate_rrho:
        # Run mRRHO calculation
        proc_xtb = Factory[XtbProc].create(Prog.XTB, "2_OPTIMIZATION")
        job_config = RRHOJobConfig(
            gfnv=config.optimization.gfnv,
            paths=config.paths,
            **config.general.model_dump(),
        )
        results = execute(
            ensemble.conformers,
            proc_xtb.xtb_rrho,
            job_config,
            "xtb",
            "optimization",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        if config.general.ignore_failed:
            ensemble.remove_conformers(lambda conf: conf.name not in results)

        for conf in ensemble:
            contributions_dict[conf.name].grrho = results[conf.name].energy

    ensemble.update_contributions(contributions_dict)

    if cut:
        threshold = (
            min(conf.gtot for conf in ensemble)
            + config.optimization.threshold / AU2KCAL
        )
        ensemble.remove_conformers(lambda conf: conf.gtot > threshold)

    # Print/write out results
    results_dict = jsonify(ensemble, config.optimization)
    _write_results(ensemble, config, results_dict)

    return results_dict


def _macrocycle_opt(
    proc: QmProc,
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    cut: bool,
):
    """
    Geometry optimization using macrocycles, whereafter every macrocycle cutting conditions are checked (if cut == True).

    :param proc: Quantum mechanics processor
    :type proc: QmProc
    :param ensemble: Ensemble data
    :type ensemble: EnsembleData
    :param config: Parts configuration
    :type config: PartsConfig
    :param client: Dask client
    :type client: Client
    :param cut: Whether to apply cutting conditions
    :type cut: bool
    :returns: Contributions dictionary
    :rtype: dict[str, Contributions]
    """
    # Second, independent ensemble containing unconverged conformers
    unconverged_ensemble = EnsembleData()

    # Copy the references to the original conformers into the unconverged_ensemble
    # This way all modifications to the MoleculeData objects will be made on the original objects
    # But the ensembles themselves are independent
    unconverged_ensemble.conformers = copy(ensemble.conformers)

    # Set up contributions_dict
    contributions_dict = {conf.name: Contributions() for conf in ensemble}

    # Set up target
    if config.optimization.xtb_opt:
        opt_job_config = XTBOptJobConfig(
            grid=GridLevel.HIGH,
            copy_mo=config.general.copy_mo,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            paths=config.paths,
            **config.optimization.model_dump(),
        )
        if config.optimization.constrain:
            opt_job_config.constraints = ensemble.constraints

        target = proc.xtb_opt
    else:
        opt_job_config = OptJobConfig(
            grid=GridLevel.HIGH,
            copy_mo=config.general.copy_mo,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            paths=config.paths,
            **config.optimization.model_dump(),
        )  # type: ignore
        target = proc.opt

    # Iterate and cut only the unconverged conformers
    ncyc = 0
    rrho_done = False
    while (
        len(unconverged_ensemble.conformers) > 0 and ncyc < config.optimization.maxcyc
    ):
        printf(
            h1(
                f"OPTIMIZATION CYCLE {ncyc // config.optimization.optcycles} ({config.optimization.optcycles} steps)"
            )
        )
        results = execute(
            unconverged_ensemble.conformers,
            target,
            opt_job_config,
            config.optimization.prog,
            "optimization",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )
        if config.general.ignore_failed:
            ensemble.remove_conformers(
                lambda conf: conf in unconverged_ensemble.conformers
                and conf.name not in results
            )
            unconverged_ensemble.remove_conformers(
                lambda conf: conf.name not in results
            )

        for conf in unconverged_ensemble:
            # Update energies
            contributions_dict[conf.name].energy = results[conf.name].energy

            # Update geometries in both ensemble objects
            conf.geom.xyz = results[conf.name].geom

        if (
            config.general.evaluate_rrho
            and ncyc + config.optimization.optcycles > 6
            and not rrho_done
        ):
            # Run mRRHO calculation
            proc_xtb: XtbProc = Factory.create(Prog.XTB, "2_OPTIMIZATION")
            job_config_rrho = RRHOJobConfig(
                gfnv=config.optimization.gfnv,
                paths=config.paths,
                **config.general.model_dump(),
            )
            results_rrho = execute(
                unconverged_ensemble.conformers,
                proc_xtb.xtb_rrho,
                job_config_rrho,
                "xtb",
                "optimization",
                client,
                ignore_failed=config.general.ignore_failed,
                balance=config.general.balance,
                copy_mo=config.general.copy_mo,
            )
            if config.general.ignore_failed:
                ensemble.remove_conformers(
                    lambda conf: conf in unconverged_ensemble.conformers
                    and conf.name not in results_rrho
                )
                unconverged_ensemble.remove_conformers(
                    lambda conf: conf.name not in results_rrho
                )

            for conf in unconverged_ensemble:
                contributions_dict[conf.name].grrho = results_rrho[conf.name].energy

            rrho_done = True

        # TODO: molbar topology check

        # Remove converged conformers
        unconverged_ensemble.remove_conformers(
            lambda conf: results[conf.name].converged
        )

        # Update contributions
        unconverged_ensemble.update_contributions(contributions_dict)

        if len(unconverged_ensemble.conformers) != 0 and cut:
            gradthr = config.optimization.gradthr

            # Take the minimum Gtot for all the conformers (converged AND unconverged) (converged are considered removed from the unconverged_ensemble)
            threshold = (
                min(
                    [conf.gtot for conf in unconverged_ensemble]
                    + [conf.gtot for conf in unconverged_ensemble.rem]
                )
                + config.optimization.threshold / AU2KCAL
            )
            unconverged_ensemble.remove_conformers(
                lambda conf: conf.gtot > threshold
                and results[conf.name].grad_norm < gradthr
            )
        # if all confs are converged, don't cut again, instead let the cutting be done in the optimization function,
        # so we can use the updated mRRHO contributions
        ncyc += config.optimization.optcycles

        _print_update(unconverged_ensemble)

    for conf in unconverged_ensemble:
        logger.warning(
            f"{conf.name} did not converge after {config.optimization.maxcyc} cycles."
        )

    ensemble.remove_conformers(lambda conf: conf in unconverged_ensemble)

    # Return most up to date information
    return contributions_dict


def _full_opt(
    proc: QmProc,
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
):
    """
    Full geometry optimization of every conformer.

    :param proc: Quantum mechanics processor
    :type proc: QmProc
    :param ensemble: Ensemble data
    :type ensemble: EnsembleData
    :param config: Parts configuration
    :type config: PartsConfig
    :param client: Dask client
    :type client: Client
    :returns: Contributions dictionary
    :rtype: dict[str, Contributions]
    """
    # Set up target
    if config.optimization.xtb_opt:
        opt_job_config = XTBOptJobConfig(
            copy_mo=config.general.copy_mo,
            grid=GridLevel.HIGH,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            paths=config.paths,
            **config.optimization.model_dump(),
        )
        if config.optimization.constrain:
            opt_job_config.constraints = ensemble.constraints

        target = proc.xtb_opt
    else:
        opt_job_config = OptJobConfig(
            copy_mo=config.general.copy_mo,
            grid=GridLevel.HIGH,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            paths=config.paths,
            **config.optimization.model_dump(),
        )  # type: ignore
        target = proc.opt

    # Set up contributions_dict
    contributions_dict = {conf.name: Contributions() for conf in ensemble}

    results = execute(
        ensemble.conformers,
        target,
        opt_job_config,
        config.optimization.prog,
        "optimization",
        client,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )

    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in results)

    for conf in ensemble:
        contributions_dict[conf.name].energy = results[conf.name].energy
        conf.geom.xyz = results[conf.name].geom

    return contributions_dict


def _write_results(
    ensemble: EnsembleData, config: PartsConfig, results_dict: dict[str, Any]
) -> None:
    """ """
    printf(h1("OPTIMIZATION RESULTS"))

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
    printf(table)

    # list the averaged free enthalpy of the ensemble
    lines: list[str] = []
    lines.append(
        "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (DFT optimized):"
    )
    lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}")

    # calculate averaged free enthalpy
    avG = sum([boltzmann_populations[conf.name] * conf.gtot for conf in ensemble])

    # calculate averaged free energy (?)
    avE = sum([boltzmann_populations[conf.name] * conf.energy for conf in ensemble])

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
    Path("2_OPTIMIZATION.json").write_text(json.dumps(results_dict, indent=4))

    ensemble.dump_xyz(Path("2_OPTIMIZATION.xyz"))


# TODO: generalize this
def jsonify(
    ensemble: EnsembleData,
    config: OptimizationConfig,
    fields: Callable[[MoleculeData], dict[str, Any]] | None = None,
):
    """
    Convert ensemble data to JSON format for optimization results.

    :param ensemble: EnsembleData object.
    :param config: OptimizationConfig object.
    :param fields: Optional callable to customize fields.
    :return: JSON-serializable dictionary.
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

    dump = DataDump(part_name="optimization")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")


def _print_update(ensemble: EnsembleData):
    # column headers
    headers = [
        "CONF#",
        "E (DFT)",
        "GmRRHO",
        "Gtot",
        "ΔGtot",
    ]

    # column units
    units = [
        "",
        "[Eh]",
        "[Eh]",
        "[Eh]",
        "[kcal/mol]",
    ]

    # variables for printmap
    gtotmin = min(conf.gtot for conf in ensemble.conformers + ensemble.rem)

    printmap = {
        "CONF#": lambda conf: conf.name,
        "E (DFT)": lambda conf: f"{conf.energy:.6f}",
        "GmRRHO": lambda conf: (f"{conf.grrho:.6f}"),
        "Gtot": lambda conf: f"{conf.gtot:.6f}",
        "ΔGtot": lambda conf: f"{(conf.gtot - gtotmin) * AU2KCAL:.2f}",
    }

    if len(ensemble.rem) > 0:
        printf(h1("Converged or removed conformers"))

        rows = [[printmap[header](conf) for header in headers] for conf in ensemble.rem]

        for i in range(len(headers)):
            headers[i] += "\n" + units[i]

        table = tabulate(
            rows,
            headers=headers,
            colalign=["left"] + ["center" for _ in headers[1:]],
            disable_numparse=True,
            numalign="decimal",
        )
        printf(table)

        for i in range(len(headers)):
            headers[i] = headers[i].split("\n")[0]

    if len(ensemble.conformers) > 0:
        printf(h1("Unconverged conformers"))

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
        printf(table)
