"""
Calculates the ensemble optical rotation spectrum.
"""

from pathlib import Path
from tabulate import tabulate
from typing import Any, Literal
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
    retemp: Literal["xtb", "cosmors"] | None = None,
) -> dict[str, Any]:
    """
    Calculation of the ensemble optical rotation spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param client: dask.distributed.Client for parallel execution.
    :param retemp: If set, reuse the reloaded optical rotation results (no calculations) and only redo the Boltzmann averaging.
        Results are written to 5_ROT_RETEMP.dat/.json.
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

    if retemp is not None:
        stored = ensemble.property_results.get("rot")
        if stored is None:
            raise RuntimeError(
                "--retemp for rot requires the previous results: add 5_ROT.json to --reload."
            )
        results = {
            conf.name: RotResult(
                rotations_velocity=[
                    (f, r) for f, r in stored[conf.name]["rotations_velocity"]
                ],
                rotations_length=[
                    (f, r) for f, r in stored[conf.name]["rotations_length"]
                ],
            )
            for conf in ensemble
        }
    else:
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
    _write_results(ensemble, config, results, results_dict, "_RETEMP" if retemp else "")

    return results_dict


def _rotation(rotations: list[tuple[float, float]], freq: float) -> float:
    """
    Rotation value at the given frequency (0.0 if not calculated).

    :param rotations: (frequency, rotation) pairs.
    :param freq: Frequency to look up (matched within 0.1).
    :return: Rotation value.
    """
    return next((r for f, r in rotations if abs(f - freq) < 0.1), 0.0)


def _average(
    ensemble: EnsembleData,
    results: dict[str, RotResult],
    frequencies: list[float],
    populations: dict[str, float],
    representation: Literal["rotations_length", "rotations_velocity"],
) -> list[tuple[float, float]]:
    """
    Boltzmann averaged rotation per frequency.

    :param ensemble: EnsembleData object containing the conformers.
    :param results: Optical rotation results per conformer.
    :param frequencies: Frequencies to average.
    :param populations: Boltzmann populations per conformer.
    :param representation: Which representation to average.
    :return: (frequency, averaged rotation) pairs.
    """
    return [
        (
            freq,
            sum(
                _rotation(getattr(results[conf.name], representation), freq)
                * populations[conf.name]
                for conf in ensemble
            ),
        )
        for freq in frequencies
    ]


def _tables(
    ensemble: EnsembleData,
    results: dict[str, RotResult],
    frequencies: list[float],
    populations: dict[str, float],
) -> tuple[str | None, str]:
    """
    Format per-conformer and Boltzmann averaged rotations as tables.

    :param ensemble: EnsembleData object containing the conformers.
    :param results: Optical rotation results per conformer.
    :param frequencies: Frequencies to print.
    :param populations: Boltzmann populations per conformer.
    :return: Length representation table (None if not available) and velocity representation table.
    """
    headers = ["CONF#\n"] + [f"{freq} nm\ndeg·dm⁻¹·(g·cm³)⁻¹" for freq in frequencies]

    def table(
        representation: Literal["rotations_length", "rotations_velocity"],
    ) -> str:
        rows = [
            [conf.name]
            + [
                f"{_rotation(getattr(results[conf.name], representation), freq):.3f}"
                for freq in frequencies
            ]
            for conf in ensemble
        ]
        rows.append(
            ["Ensemble avg."]
            + [
                f"{avg:.3f}"
                for _, avg in _average(
                    ensemble, results, frequencies, populations, representation
                )
            ]
        )
        return tabulate(
            rows,
            headers=headers,
            colalign=["left"] + ["center" for _ in headers[1:]],
            disable_numparse=True,
            numalign="decimal",
        )

    # Length representation is only available for some programs, velocity always
    has_length = all(result.rotations_length for result in results.values())
    return (
        table("rotations_length") if has_length else None,
        table("rotations_velocity"),
    )


def _write_results(
    ensemble: EnsembleData,
    config: PartsConfig,
    results: dict[str, RotResult],
    results_dict: dict[str, Any],
    suffix: str,
):
    """
    Print the optical rotation tables and write them to 5_ROT<suffix>.dat/.json, including per-temperature averages with multitemp.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param results: Optical rotation results per conformer.
    :param results_dict: JSON results, per-temperature averages and populations are added with multitemp.
    :param suffix: File name suffix ("_RETEMP" for --retemp runs).
    """
    frequencies = config.rot.freq

    def text(populations: dict[str, float]) -> str:
        table_length, table_velocity = _tables(
            ensemble, results, frequencies, populations
        )
        parts = []
        if table_length is not None:
            parts += [
                h1("Optical Rotation Values (Length Representation)"),
                table_length,
            ]
        parts += [
            h1("Optical Rotation Values (Velocity Representation)"),
            table_velocity,
        ]
        return "\n".join(parts)

    out = text(ensemble.get_populations(config.general.temperature))
    printf(out)
    printf("".ljust(int(PLENGTH), "-"))

    temperatures = ensemble.temperature_range(config.general, "rot")
    if temperatures:
        results_dict["populations_t"] = {}
        results_dict["averages_t"] = {}
    representations: list[Literal["rotations_length", "rotations_velocity"]] = [
        "rotations_velocity"
    ]
    if all(result.rotations_length for result in results.values()):
        representations.insert(0, "rotations_length")
    for t in temperatures:
        populations = ensemble.get_populations(t)
        out += f"\n\nT = {t} K\n" + text(populations)
        results_dict["populations_t"][str(t)] = populations
        results_dict["averages_t"][str(t)] = {
            representation: [
                [f, avg]
                for f, avg in _average(
                    ensemble, results, frequencies, populations, representation
                )
            ]
            for representation in representations
        }

    Path(f"5_ROT{suffix}.dat").write_text(out)

    # Additionally, write results in json format
    Path(f"5_ROT{suffix}.json").write_text(json.dumps(results_dict, indent=4))


def jsonify(ensemble: EnsembleData, config: RotConfig, results: dict[str, RotResult]):
    """
    Prepare optical rotation results for JSON serialization.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: RotConfig object with optical rotation configuration settings.
    :param results: Dictionary of RotResult objects for each conformer.
    :return: Dictionary ready for JSON serialization.
    """

    def per_conf(conf: MoleculeData) -> dict[str, dict[str, Any]]:
        return {
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
                **conf.temperature_data(),
                "rotations_length": results[conf.name].rotations_length,
                "rotations_velocity": results[conf.name].rotations_velocity,
            }
        }

    dump = DataDump(part_name="rot")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
