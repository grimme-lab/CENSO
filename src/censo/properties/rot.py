"""
Calculates the ensemble optical rotation spectrum.
"""

from pathlib import Path
from collections.abc import Callable
from tabulate import tabulate
import json
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.managers import SyncManager

from ..ensemble import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.job_config import RotJobConfig, RotResult
from ..config.parts import RotConfig
from ..config.parallel_config import ParallelConfig
from ..params import GridLevel, PLENGTH
from ..parallel import execute, ResourceMonitor
from ..utilities import printf, Factory, h1, h2, timeit, DataDump
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


@timeit
def rot(
    ensemble: EnsembleData,
    config: PartsConfig,
    parallel_config: ParallelConfig | None,
    executor: ProcessPoolExecutor | None = None,
    manager: SyncManager | None = None,
    resource_monitor: ResourceMonitor | None = None,
):
    """
    Calculation of the ensemble optical rotation spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param parallel_config: ParallelConfig object for parallel execution.
    :return: None
    """
    printf(h2("ROT"))

    if executor is None or manager is None or resource_monitor is None:
        raise ValueError("executor, manager, and resource_monitor must be provided")

    config.model_validate(config, context={"check": "rot"})

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc: QmProc = Factory.create(config.rot.prog, "5_ROT")

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
        parallel_config,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
        executor=executor,  # type: ignore
        manager=manager,  # type: ignore
        resource_monitor=resource_monitor,  # type: ignore
    )
    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in results)

    _write_results(ensemble, config, results)


def _write_results(
    ensemble: EnsembleData, config: PartsConfig, results: dict[str, RotResult]
):
    boltzmann_populations = ensemble.get_populations(config.general.temperature)

    # Get all frequencies from job config
    frequencies = config.rot.freq

    # Collect rotation data for each conformer
    rotation_data_length: list[list[str | float]] = []
    rotation_data_velocity: list[list[str | float]] = []
    for conf in ensemble:
        conf_data_length: list[str | float] = [conf.name]
        conf_data_velocity: list[str | float] = [conf.name]
        # Create dictionaries to map frequencies to rotation values
        freq_to_rotation_length = {
            freq: rotation for freq, rotation in results[conf.name].rotations_length
        }
        freq_to_rotation_velocity = {
            freq: rotation for freq, rotation in results[conf.name].rotations_velocity
        }
        # Add rotation values in the order of frequencies
        for freq in frequencies:
            rotation_val_length = freq_to_rotation_length.get(freq, 0.0)
            rotation_val_velocity = freq_to_rotation_velocity.get(freq, 0.0)
            conf_data_length.append(rotation_val_length)
            conf_data_velocity.append(rotation_val_velocity)
        rotation_data_length.append(conf_data_length)
        rotation_data_velocity.append(conf_data_velocity)

    filepath = Path("5_ROT.out")

    # Create table headers with frequencies as columns
    headers = ["CONF#"] + [f"{freq:.1f} nm" for freq in frequencies]
    units = [""] + ["[deg·mL·g⁻¹·dm⁻¹]" for _ in frequencies]

    # Create table rows for length representation
    rows_length = []
    for conf_data in rotation_data_length:
        row = [conf_data[0]]  # conformer name
        for i, rotation_val in enumerate(conf_data[1:]):
            row.append(f"{rotation_val:.3f}")
        rows_length.append(row)

    # Create table rows for velocity representation
    rows_velocity = []
    for conf_data in rotation_data_velocity:
        row = [conf_data[0]]  # conformer name
        for i, rotation_val in enumerate(conf_data[1:]):
            row.append(f"{rotation_val:.3f}")
        rows_velocity.append(row)

    # Add ensemble averaged values for length representation
    ensemble_avg_length = ["Ensemble avg."]
    for i, freq in enumerate(frequencies):
        avg_rotation = sum(
            rotation * boltzmann_populations[conf.name]
            for conf in ensemble
            for f, rotation in results[conf.name].rotations_length
            if abs(f - freq) < 0.1  # Match frequency within tolerance
        )
        ensemble_avg_length.append(f"{avg_rotation:.3f}")
    rows_length.append(ensemble_avg_length)

    # Add ensemble averaged values for velocity representation
    ensemble_avg_velocity = ["Ensemble avg."]
    for i, freq in enumerate(frequencies):
        avg_rotation = sum(
            rotation * boltzmann_populations[conf.name]
            for conf in ensemble
            for f, rotation in results[conf.name].rotations_velocity
            if abs(f - freq) < 0.1  # Match frequency within tolerance
        )
        ensemble_avg_velocity.append(f"{avg_rotation:.3f}")
    rows_velocity.append(ensemble_avg_velocity)

    # Add units to headers
    for i in range(len(headers)):
        headers[i] += "\n" + units[i]

    table_length = tabulate(
        rows_length,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )

    table_velocity = tabulate(
        rows_velocity,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )

    # Print everything
    printf(h1("Optical Rotation Values (Length Representation)"))
    printf(table_length)
    printf(h1("Optical Rotation Values (Velocity Representation)"))
    printf(table_velocity)
    printf("".ljust(int(PLENGTH), "-"))

    filepath.write_text(table_length)
    filepath.write_text(table_velocity)

    # Additionally, write results in json format
    Path("5_ROT.json").write_text(
        json.dumps(jsonify(ensemble, config.rot, results), indent=4)
    )


def jsonify(ensemble: EnsembleData, config: RotConfig, results: dict[str, RotResult]):
    """
    Prepare optical rotation results for JSON serialization.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: RotConfig object with optical rotation configuration settings.
    :param results: Dictionary of RotResult objects for each conformer.
    :return: Dictionary ready for JSON serialization.
    """
    per_conf: Callable[
        [MoleculeData],
        dict[
            str,
            dict[
                str,
                float | list[tuple[float, float]],
            ],
        ],
    ] = lambda conf: {
        conf.name: {
            "energy": conf.energy,
            "gsolv": conf.gsolv,
            "grrho": conf.grrho,
            "rotations_length": results[conf.name].rotations_length,
            "rotations_velocity": results[conf.name].rotations_velocity,
        }
    }

    dump = DataDump(part_name="rot")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
