"""
Calculates the ensemble UV/Vis spectrum.
"""

from typing import Any
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..config import PartsConfig
from ..config.job_config import UVVisJobConfig
from ..params import GridLevel
from ..parallel import execute
from ..output import format_table, write_output
from ..utilities import printf, Factory, h1, h2
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


def uvvis(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
) -> dict[str, Any]:
    """
    Calculation of the ensemble UV/Vis spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :return: JSON-serializable dictionary of UV/Vis results.
    """
    printf(h2("UVVIS"))

    config = PartsConfig.model_validate(
        config.model_dump(), context={"check": ["uvvis"]}
    )

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc = Factory[QmProc].create(config.uvvis.prog, "6_UVVIS")

    # Run UVVis calculations
    # TODO: if some calculations fail we would need to recalculate boltzmann populations
    job_config = UVVisJobConfig(
        copy_mo=config.general.copy_mo,
        grid=GridLevel.VERY_HIGH,
        gas_phase=config.general.gas_phase,
        solvent=config.general.solvent,
        paths=config.paths,
        **config.uvvis.model_dump(),
    )

    results = execute(
        ensemble.conformers,
        proc.uvvis,
        job_config,
        config.uvvis.prog,
        "uvvis",
        client,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )
    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in results)

    populations = ensemble.get_populations(config.general.temperature)
    excitations = [
        (
            conf.name,
            excitation["wavelength"],
            excitation["osc_str"] * populations[conf.name],
        )
        for conf in ensemble
        for excitation in results[conf.name].excitations
    ]
    table = format_table(
        [
            ("λ", "[nm]", lambda exc: f"{exc[1]:.2f}"),
            ("ε_max", "", lambda exc: f"{exc[2]:.6f}"),
            ("Origin. CONF#", "", lambda exc: f"{exc[0]}"),
        ],
        excitations,
    )
    return write_output(
        "6_UVVIS",
        ensemble,
        config,
        "\n".join([h1("Averaged UVVis excitations"), table]),
        lambda conf: {"excitations": results[conf.name].excitations},
    )
