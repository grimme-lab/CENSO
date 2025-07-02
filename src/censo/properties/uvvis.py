"""
Calculates the ensemble UV/Vis spectrum.
"""

from pathlib import Path
from collections.abc import Callable
from tabulate import tabulate
import json

from ..ensembledata import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.job_config import UVVisJobConfig
from ..config.parts import UVVisConfig
from ..config.parallel_config import ParallelConfig
from ..params import GridLevel
from ..parallel import execute
from ..data import UVVisResult
from ..utilities import printf, Factory, h1, h2, timeit, DataDump
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


@timeit
def uvvis(
    ensemble: EnsembleData,
    config: PartsConfig,
    parallel_config: ParallelConfig | None,
):
    """
    Calculation of the ensemble UV/Vis spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.
    """
    printf(h2("UVVIS"))

    # Assert that all conformers have populations defined
    try:
        assert all(conf.bmw for conf in ensemble)
    except AssertionError:
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc: QmProc = Factory[QmProc].create(config.uvvis.prog, "6_UVVIS")

    # Run UVVis calculations
    # TODO: if some calculations fail we would need to recalculate boltzmann populations
    job_config = UVVisJobConfig(
        copy_mo=config.general.copy_mo,
        grid=GridLevel.VERY_HIGH,
        gas_phase=config.general.gas_phase,
        solvent=config.general.solvent,
        **config.uvvis.model_dump(),
    )

    results = execute(
        ensemble.conformers,
        proc.uvvis,
        job_config,
        config.uvvis.prog,
        "uvvis",
        parallel_config=parallel_config,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )

    _write_results(ensemble, config, results)


def _write_results(
    ensemble: EnsembleData, config: PartsConfig, results: dict[str, UVVisResult]
):
    # Average excitations
    excitations = [
        (conf.name, excitation["wavelength"], excitation["osc_str"] * conf.bmw)
        for conf in ensemble
        for excitation in results[conf.name].excitations
    ]

    filepath = Path("6_UVVIS.out")

    # Create table
    headers = ["λ", "ε_max", "Origin. CONF#"]

    units = ["[nm]", "", ""]

    printmap = {
        "λ": lambda exc: f"{exc[1]:.2f}",
        "ε_max": lambda exc: f"{exc[2]:.6f}",
        "Origin. CONF#": lambda exc: f"{exc[0]}",
    }

    rows = [[printmap[header](exc) for header in headers] for exc in excitations]

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
    printf(h1("Averaged UVVis excitations"))
    printf(table)

    filepath.write_text(table)

    # Additionally, write results in json format
    Path("6_UVVIS.json").write_text(
        json.dumps(jsonify(ensemble, config.uvvis, results), indent=4)
    )


def jsonify(
    ensemble: EnsembleData, config: UVVisConfig, results: dict[str, UVVisResult]
):
    per_conf: Callable[
        [MoleculeData],
        dict[
            str,
            dict[
                str,
                float | list[dict[str, float]],
            ],
        ],
    ] = lambda conf: {
        conf.name: {
            "energy": conf.energy,
            "gsolv": conf.gsolv,
            "grrho": conf.grrho,
            "excitations": results[conf.name].excitations,
        }
    }

    dump = DataDump(part_name="uvvis")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
