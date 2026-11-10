"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv using higher level solvation models and include the mRRHO contributions.
"""

from math import exp
from statistics import stdev
from typing import Any
from dask.distributed import Client

from ..molecules import Contributions
from ..logging import setup_logger
from ..parallel import execute
from ..params import AU2KCAL, GridLevel, Prog
from ..utilities import h2, printf, Factory
from ..config import PartsConfig
from ..config.job_config import RRHOJobConfig, SPJobConfig
from ..ensemble import EnsembleData
from ..output import write_ranking
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

    config = PartsConfig.model_validate(
        config.model_dump(), context={"check": ["screening"]}
    )

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

    return write_ranking(
        "1_SCREENING",
        ensemble,
        config,
        "SCREENING SINGLE-POINT (+ mRRHO) RESULTS",
        geometries="not DFT optimized",
    )
