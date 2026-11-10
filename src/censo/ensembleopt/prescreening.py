from typing import Any
from dask.distributed import Client

from ..processing.xtb_processor import XtbProc
from ..molecules import Contributions
from ..ensemble import EnsembleData
from ..utilities import Factory, h2, printf
from ..config import PartsConfig
from ..parallel import execute
from ..processing import QmProc
from ..output import write_ranking
from ..params import AU2KCAL, GridLevel, Prog
from ..config.job_config import SPJobConfig, XTBJobConfig
from ..logging import setup_logger

logger = setup_logger(__name__)


def prescreening(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    cut: bool = True,
) -> dict[str, Any]:
    """
    This implements a cheap prescreening step using low-cost DFT and possibly
    solvation contributions calculated using xtb.

    The list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).

    :param ensemble: EnsembleData object containing the conformers.
    :type ensemble: EnsembleData
    :param config: PartsConfig object with configuration settings.
    :type config: PartsConfig
    :param client: dask.distributed.Client for parallel execution.
    :type client: Client
    :param cut: Whether to apply cutting conditions.
    :type cut: bool
    :return: JSON-serializable dictionary of prescreening results.
    :rtype: dict
    """
    printf(h2("PRESCREENING"))

    config = PartsConfig.model_validate(
        config.model_dump(), context={"check": ["prescreening"]}
    )

    # Setup processor and target
    proc = Factory[QmProc].create(config.prescreening.prog, "0_PRESCREENING")

    contributions_dict = {conf.name: Contributions() for conf in ensemble}
    if not config.general.gas_phase:
        # Calculate Gsolv using xtb
        xtb_proc = Factory[XtbProc].create(Prog.XTB, "0_PRESCREENING")
        gsolv_job_config = XTBJobConfig(
            gfnv=config.prescreening.gfnv,
            solvent=config.general.solvent,
            sm_rrho=config.general.sm_rrho,
            gas_phase=False,
            paths=config.paths,
        )
        gsolv_results = execute(
            ensemble.conformers,
            xtb_proc.gsolv,
            gsolv_job_config,
            "xtb",
            "prescreening",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        if config.general.ignore_failed:
            ensemble.remove_conformers(lambda conf: conf.name not in gsolv_results)

        for conf in ensemble:
            contributions_dict[conf.name].gsolv = gsolv_results[conf.name].gsolv

    # Calculate gas-phase single-point
    sp_job_config = SPJobConfig(
        copy_mo=config.general.copy_mo,
        func=config.prescreening.func,
        basis=config.prescreening.basis,
        grid=GridLevel.LOW,
        template=config.prescreening.template,
        gas_phase=True,
        paths=config.paths,
    )
    sp_results = execute(
        ensemble.conformers,
        proc.sp,
        sp_job_config,
        config.prescreening.prog,
        "prescreening",
        client,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )

    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in sp_results)

    for conf in ensemble:
        contributions_dict[conf.name].energy = sp_results[conf.name].energy

    # Update molecules
    ensemble.update_contributions(contributions_dict)

    if cut:
        threshold = (
            min(conf.gtot for conf in ensemble)
            + config.prescreening.threshold / AU2KCAL
        )
        ensemble.remove_conformers(cond=lambda conf: conf.gtot > threshold)

    return write_ranking(
        "0_PRESCREENING",
        ensemble,
        config,
        "PRESCREENING SINGLE-POINT RESULTS",
        columns=(
            "CONF#",
            "E (DFT)",
            "ΔE (DFT)",
            "ΔGsolv (xTB)",
            "Gtot",
            "ΔGtot",
            "Boltzmann weight",
        ),
        geometries="not DFT optimized",
    )
