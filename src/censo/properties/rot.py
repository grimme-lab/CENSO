"""
Calculates the ensemble optical rotation spectrum.
"""

from operator import itemgetter
from typing import Any
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..config import PartsConfig
from ..config.job_config import RotJobConfig
from ..params import GridLevel
from ..parallel import execute
from ..output import Column, format_table, write_output
from ..utilities import printf, Factory, h1, h2
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


def rot(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
) -> dict[str, Any]:
    """
    Calculation of the ensemble optical rotation spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :return: JSON-serializable dictionary of optical rotation results.
    """
    printf(h2("ROT"))

    config = PartsConfig.model_validate(config.model_dump(), context={"check": ["rot"]})

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc = Factory[QmProc].create(config.rot.prog, "5_ROT")

    # Run optical rotation calculations
    job_config = RotJobConfig(
        copy_mo=config.general.copy_mo,
        grid=GridLevel.VERY_HIGH,
        gas_phase=True,
        paths=config.paths,
        func=config.rot.func,
        basis=config.rot.basis,
        template=config.rot.template,
        temperature=config.general.temperature,
        freq=config.rot.freq,
    )

    results = execute(
        ensemble.conformers,
        proc.rot,
        job_config,
        config.rot.prog,
        "rot",
        client,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )
    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in results)

    populations = ensemble.get_populations(config.general.temperature)
    blocks: list[str] = []
    # Length representation is only printed if available for all conformers
    if all(result.rotations_length for result in results.values()):
        blocks += [
            h1("Optical Rotation Values (Length Representation)"),
            _rotation_table(
                ensemble,
                {conf.name: results[conf.name].rotations_length for conf in ensemble},
                config.rot.freq,
                populations,
            ),
        ]
    blocks += [
        h1("Optical Rotation Values (Velocity Representation)"),
        _rotation_table(
            ensemble,
            {conf.name: results[conf.name].rotations_velocity for conf in ensemble},
            config.rot.freq,
            populations,
        ),
    ]
    return write_output(
        "5_ROT",
        ensemble,
        config,
        "\n".join(blocks),
        lambda conf: {
            "rotations_length": results[conf.name].rotations_length,
            "rotations_velocity": results[conf.name].rotations_velocity,
        },
    )


def _rotation_table(
    ensemble: EnsembleData,
    rotations: dict[str, list[tuple[float, float]]],
    frequencies: list[float],
    populations: dict[str, float],
) -> str:
    """
    Tabulate the optical rotation of each conformer and the Boltzmann average at each frequency.

    :param ensemble: Conformers to tabulate.
    :param rotations: (frequency, rotation) pairs of one representation per conformer name.
    :param frequencies: Frequencies to tabulate, one column each.
    :param populations: Boltzmann population per conformer name.
    :return: The table.
    """
    values = [
        (
            conf.name,
            [
                next((r for f, r in rotations[conf.name] if abs(f - freq) < 0.1), 0.0)
                for freq in frequencies
            ],
        )
        for conf in ensemble
    ]
    averages = [
        sum(
            r * populations[conf.name]
            for conf in ensemble
            for f, r in rotations[conf.name]
            if abs(f - freq) < 0.1
        )
        for freq in frequencies
    ]
    rows = [
        [name, *(f"{value:.3f}" for value in row)]
        for name, row in [*values, ("Ensemble avg.", averages)]
    ]
    columns: list[Column[list[str]]] = [("CONF#", "", itemgetter(0))]
    columns += [
        (f"{freq} nm", "deg·dm⁻¹·(g·cm³)⁻¹", itemgetter(i))
        for i, freq in enumerate(frequencies, start=1)
    ]
    return format_table(columns, rows)
