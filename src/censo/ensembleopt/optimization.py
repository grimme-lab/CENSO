"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""

from copy import copy
from typing import Any
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..molecules import MoleculeData, Contributions
from ..processing import QmProc, XtbProc
from ..parallel import execute
from ..params import AU2KCAL, GridLevel, Prog
from ..config import PartsConfig
from ..config.job_config import RRHOJobConfig, OptJobConfig, XTBOptJobConfig
from ..output import Column, format_table, write_ranking
from ..utilities import (
    printf,
    h1,
    Factory,
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

    config = PartsConfig.model_validate(
        config.model_dump(), context={"check": ["optimization"]}
    )

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

    return write_ranking(
        "2_OPTIMIZATION",
        ensemble,
        config,
        "OPTIMIZATION RESULTS",
        geometries="DFT optimized",
    )


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

        # Update contributions
        unconverged_ensemble.update_contributions(contributions_dict)

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


def _print_update(ensemble: EnsembleData) -> None:
    """
    Print converged/removed and unconverged conformers after an optimization macrocycle.

    :param ensemble: Ensemble of unconverged conformers; removed ones are in ensemble.rem.
    """
    gtotmin = min(conf.gtot for conf in ensemble.conformers + ensemble.rem)
    columns: list[Column[MoleculeData]] = [
        ("CONF#", "", lambda conf: conf.name),
        ("E (DFT)", "[Eh]", lambda conf: f"{conf.energy:.6f}"),
        ("GmRRHO", "[Eh]", lambda conf: f"{conf.grrho:.6f}"),
        ("Gtot", "[Eh]", lambda conf: f"{conf.gtot:.6f}"),
        ("ΔGtot", "[kcal/mol]", lambda conf: f"{(conf.gtot - gtotmin) * AU2KCAL:.2f}"),
    ]
    if len(ensemble.rem) > 0:
        printf(h1("Converged or removed conformers"))
        printf(format_table(columns, ensemble.rem))
    if len(ensemble.conformers) > 0:
        printf(h1("Unconverged conformers"))
        printf(format_table(columns, ensemble.conformers))
