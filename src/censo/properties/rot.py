"""
Calculates the ensemble optical rotation spectrum.
"""

from pathlib import Path
from collections.abc import Callable
from tabulate import tabulate
import json

from ..ensemble import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.job_config import RotJobConfig, RotResult
from ..config.parts import RotConfig
from ..config.parallel_config import ParallelConfig
from ..params import GridLevel
from ..parallel import execute
from ..utilities import printf, Factory, h1, h2, timeit, DataDump
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


@timeit
def rot(
    ensemble: EnsembleData,
    config: PartsConfig,
    parallel_config: ParallelConfig | None,
):
    """
    Calculation of the ensemble optical rotation spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.
    """
    printf(h2("ROT"))

    config.model_validate(config, context={"check": "rot"})

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc: QmProc = Factory[QmProc].create(config.rot.prog, "5_ROT")

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
    rotation_data: list[list[str | float]] = []
    for conf in ensemble:
        conf_data: list[str | float] = [conf.name]
        # Create a dictionary to map frequencies to rotation values
        freq_to_rotation = {
            freq: rotation for freq, rotation in results[conf.name].rotations
        }
        # Add rotation values in the order of frequencies
        for freq in frequencies:
            rotation_val = freq_to_rotation.get(freq, 0.0)
            conf_data.append(rotation_val)
        rotation_data.append(conf_data)

    filepath = Path("6_ROT.out")

    # Create table headers with frequencies as columns
    headers = ["CONF#"] + [f"{freq:.1f} nm" for freq in frequencies]
    units = [""] + ["[deg·mL·g⁻¹·dm⁻¹]" for _ in frequencies]

    # Create table rows
    rows = []
    for conf_data in rotation_data:
        row = [conf_data[0]]  # conformer name
        for i, rotation_val in enumerate(conf_data[1:]):
            row.append(f"{rotation_val:.3f}")
        rows.append(row)

    # Add ensemble averaged values
    ensemble_avg = ["Ensemble avg."]
    for i, freq in enumerate(frequencies):
        avg_rotation = sum(
            rotation * boltzmann_populations[conf.name]
            for conf in ensemble
            for f, rotation in results[conf.name].rotations
            if abs(f - freq) < 0.1  # Match frequency within tolerance
        )
        ensemble_avg.append(f"{avg_rotation:.3f}")
    rows.append(ensemble_avg)

    # Add units to headers
    for i in range(len(headers)):
        headers[i] += "\n" + units[i]

    table = tabulate(
        rows,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )

    # Print everything
    printf(h1("Optical Rotation Values"))
    printf(table)

    filepath.write_text(table)

    # Additionally, write results in json format
    Path("6_ROT.json").write_text(
        json.dumps(jsonify(ensemble, config.rot, results), indent=4)
    )


def jsonify(ensemble: EnsembleData, config: RotConfig, results: dict[str, RotResult]):
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
            "rotations": results[conf.name].rotations,
        }
    }

    dump = DataDump(part_name="rot")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
