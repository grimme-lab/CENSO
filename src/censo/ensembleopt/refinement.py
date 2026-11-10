from typing import Any
from dask.distributed import Client

from ..molecules import Contributions
from ..logging import setup_logger
from ..parallel import execute
from ..params import GridLevel, Prog
from ..utilities import h2, printf, Factory
from ..config import PartsConfig
from ..config.job_config import RRHOJobConfig, SPJobConfig
from ..ensemble import EnsembleData
from ..output import write_ranking
from ..processing import QmProc, XtbProc

logger = setup_logger(__name__)


# TODO: this is not very DRY considering the screening function
def refinement(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    cut: bool = True,
) -> dict[str, Any]:
    """
    Basically the same as screening, however here we use a Boltzmann population cutoff instead of kcal cutoff.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param client: dask.distributed.Client for parallel execution.
    :param cut: Whether to apply cutting conditions.
    :return: JSON-serializable dictionary of refinement results.
    """
    printf(h2("REFINEMENT"))

    config = PartsConfig.model_validate(
        config.model_dump(), context={"check": ["refinement"]}
    )

    # Setup processor and target
    proc = Factory[QmProc].create(config.refinement.prog, "3_REFINEMENT")

    contributions_dict = {conf.name: Contributions() for conf in ensemble}
    if not config.general.gas_phase and not config.refinement.gsolv_included:
        # Calculate Gsolv using qm
        job_config = SPJobConfig(
            copy_mo=config.general.copy_mo,
            func=config.refinement.func,
            basis=config.refinement.basis,
            grid=GridLevel.VERY_HIGH,
            template=config.refinement.template,
            gas_phase=False,
            solvent=config.general.solvent,
            sm=config.refinement.sm,
            # multitemp=config.general.multitemp,
            # trange=config.general.trange,
            temperature=config.general.temperature,
            paths=config.paths,
        )
        results = execute(
            ensemble.conformers,
            proc.gsolv,
            job_config,
            config.refinement.prog,
            "refinement",
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
            copy_mo=config.general.copy_mo,
            func=config.refinement.func,
            basis=config.refinement.basis,
            grid=GridLevel.VERY_HIGH,
            template=config.refinement.template,
            gas_phase=config.general.gas_phase,
            solvent=config.general.solvent,
            sm=config.refinement.sm,
            paths=config.paths,
        )
        sp_results = execute(
            ensemble.conformers,
            proc.sp,
            sp_job_config,
            config.refinement.prog,
            "refinement",
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
        proc_xtb = Factory[XtbProc].create(Prog.XTB, "3_REFINEMENT")
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
            "refinement",
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
        # Prepare Boltzmann populations
        boltzmann_populations = ensemble.get_populations(config.general.temperature)
        ensemble.conformers.sort(
            key=lambda conf: boltzmann_populations[conf.name], reverse=True
        )

        # Threshold in cumulative % pop
        threshold = config.refinement.threshold

        ensemble.remove_conformers(
            cond=lambda conf: sum(
                boltzmann_populations[c.name]
                for c in ensemble.conformers[: ensemble.conformers.index(conf)]
            )
            > threshold
        )

    return write_ranking(
        "3_REFINEMENT",
        ensemble,
        config,
        "REFINEMENT SINGLE-POINT (+ mRRHO) RESULTS",
    )
