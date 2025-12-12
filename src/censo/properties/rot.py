"""
Calculates the ensemble optical rotation spectrum.
"""

from pathlib import Path
from tabulate import tabulate
from typing import Any
import json
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.job_config import RotJobConfig
from ..processing.results import RotResult
from ..config.parts import RotConfig
from ..params import GridLevel, PLENGTH
from ..parallel import execute
from ..utilities import printf, Factory, h1, h2, DataDump
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

    config = PartsConfig.model_validate(config, context={"check": ["rot"]})

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

    results_dict = jsonify(ensemble, config.rot, results)
    _write_results(ensemble, config, results, results_dict)

    return results_dict


def _write_results(
    ensemble: EnsembleData,
    config: PartsConfig,
    results: dict[str, RotResult],
    results_dict: dict[str, Any],
):
    boltzmann_populations = ensemble.get_populations(config.general.temperature)

    # Get all frequencies from job config
    frequencies = config.rot.freq

    # Define headers and units
    headers = ["CONF#"] + [f"{freq} nm" for freq in frequencies]
    units = [""] + ["deg·dm⁻¹·(g·cm³)⁻¹" for _ in frequencies]
    filepath = Path("5_ROT.dat")

    # Collect rotation data for each conformer
    rows_length: list[list[str]] = []
    rows_velocity: list[list[str]] = []

    # Only print if lengths representation is available
    if all(result.rotations_length for result in results.values()):
        for conf in ensemble:
            row = [conf.name]
            for freq in frequencies:
                rotation = next(
                    (
                        rotation
                        for f, rotation in results[conf.name].rotations_length
                        if abs(f - freq) < 0.1
                    ),
                    0.0,
                )
                row.append(f"{rotation:.3f}")
            rows_length.append(row)

    # Velocity represtntation should always be available
    for conf in ensemble:
        row = [conf.name]
        for freq in frequencies:
            rotation = next(
                (
                    rotation
                    for f, rotation in results[conf.name].rotations_velocity
                    if abs(f - freq) < 0.1
                ),
                0.0,
            )
            row.append(f"{rotation:.3f}")
        rows_velocity.append(row)

    # Add ensemble averaged values for length representation
    if all(result.rotations_length for result in results.values()):
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

    # Write length representation if available
    if rows_length:
        table_length = tabulate(
            rows_length,
            headers=headers,
            colalign=["left"] + ["center" for _ in headers[1:]],
            disable_numparse=True,
            numalign="decimal",
        )
        printf(h1("Optical Rotation Values (Length Representation)"))
        printf(table_length)
        filepath.write_text(table_length)

    # Write velocity representation
    table_velocity = tabulate(
        rows_velocity,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )
    printf(h1("Optical Rotation Values (Velocity Representation)"))
    printf(table_velocity)
    filepath.write_text(table_velocity)

    printf("".ljust(int(PLENGTH), "-"))

    # Additionally, write results in json format
    Path("5_ROT.json").write_text(json.dumps(results_dict, indent=4))


def jsonify(ensemble: EnsembleData, config: RotConfig, results: dict[str, RotResult]):
    """
    Prepare optical rotation results for JSON serialization.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: RotConfig object with optical rotation configuration settings.
    :param results: Dictionary of RotResult objects for each conformer.
    :return: Dictionary ready for JSON serialization.
    """

    def per_conf(
        conf: MoleculeData,
    ) -> dict[str, dict[str, float | list[tuple[float, float]]]]:
        return {
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
