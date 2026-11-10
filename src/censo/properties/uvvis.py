"""
Calculates the ensemble UV/Vis spectrum.
"""

from pathlib import Path
from tabulate import tabulate
from typing import Any, Literal
import json
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.job_config import UVVisJobConfig
from ..config.parts import UVVisConfig
from ..params import GridLevel, PLENGTH
from ..parallel import execute
from ..processing.results import UVVisResult
from ..utilities import printf, Factory, h1, h2, DataDump
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


def uvvis(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    retemp: Literal["xtb", "cosmors"] | None = None,
) -> dict[str, Any]:
    """
    Calculation of the ensemble UV/Vis spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param client: dask.distributed.Client for parallel execution.
    :param retemp: If set, reuse the reloaded UV/Vis results (no calculations) and only redo the Boltzmann averaging.
        Results are written to 6_UVVIS_RETEMP.out/.json.
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

    if retemp is not None:
        stored = ensemble.property_results.get("uvvis")
        if stored is None:
            raise RuntimeError(
                "--retemp for uvvis requires the previous results: add 6_UVVIS.json to --reload."
            )
        results = {
            conf.name: UVVisResult(excitations=stored[conf.name]["excitations"])
            for conf in ensemble
        }
    else:
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

    results_dict = jsonify(ensemble, config.uvvis, results)
    _write_results(ensemble, config, results, results_dict, "_RETEMP" if retemp else "")

    return results_dict


def _table(
    ensemble: EnsembleData,
    results: dict[str, UVVisResult],
    populations: dict[str, float],
) -> tuple[str, list[tuple[str, float, float]]]:
    """
    Boltzmann weight all excitations and format them as a table.

    :param ensemble: EnsembleData object containing the conformers.
    :param results: UV/Vis results per conformer.
    :param populations: Boltzmann populations per conformer.
    :return: Table and (conformer, wavelength, weighted oscillator strength) per excitation.
    """
    excitations = [
        (
            conf.name,
            excitation["wavelength"],
            excitation["osc_str"] * populations[conf.name],
        )
        for conf in ensemble
        for excitation in results[conf.name].excitations
    ]

    headers = ["λ\n[nm]", "ε_max\n", "Origin. CONF#\n"]
    rows = [[f"{exc[1]:.2f}", f"{exc[2]:.6f}", f"{exc[0]}"] for exc in excitations]

    table = tabulate(
        rows,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )
    return table, excitations


def _write_results(
    ensemble: EnsembleData,
    config: PartsConfig,
    results: dict[str, UVVisResult],
    results_dict: dict[str, Any],
    suffix: str,
):
    """
    Print the averaged UV/Vis excitations and write them to 6_UVVIS<suffix>.out/.json, including per-temperature averages with multitemp.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param results: UV/Vis results per conformer.
    :param results_dict: JSON results, per-temperature averages and populations are added with multitemp.
    :param suffix: File name suffix ("_RETEMP" for --retemp runs).
    """
    table, _ = _table(
        ensemble, results, ensemble.get_populations(config.general.temperature)
    )

    # Print everything
    printf(h1("Averaged UVVis excitations"))
    printf(table)
    printf("".ljust(int(PLENGTH), "-"))

    text = table
    temperatures = ensemble.temperature_range(config.general, "uvvis")
    if temperatures:
        results_dict["populations_t"] = {}
        results_dict["averages_t"] = {}
    for t in temperatures:
        populations = ensemble.get_populations(t)
        table_t, excitations = _table(ensemble, results, populations)
        text += f"\n\nT = {t} K\n" + table_t
        results_dict["populations_t"][str(t)] = populations
        results_dict["averages_t"][str(t)] = {
            "excitations": [
                {"conf": name, "wavelength": wavelength, "osc_str": osc_str}
                for name, wavelength, osc_str in excitations
            ]
        }

    Path(f"6_UVVIS{suffix}.out").write_text(text)

    # Additionally, write results in json format
    Path(f"6_UVVIS{suffix}.json").write_text(json.dumps(results_dict, indent=4))


def jsonify(
    ensemble: EnsembleData, config: UVVisConfig, results: dict[str, UVVisResult]
):
    """
    Prepare UV/Vis results for JSON serialization.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: UVVisConfig object with UV/Vis configuration settings.
    :param results: Dictionary of UVVisResult objects for each conformer.
    :return: Dictionary ready for JSON serialization.
    """

    def per_conf(conf: MoleculeData) -> dict[str, dict[str, Any]]:
        return {
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
                **conf.temperature_data(),
                "excitations": results[conf.name].excitations,
            }
        }

    dump = DataDump(part_name="uvvis")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
