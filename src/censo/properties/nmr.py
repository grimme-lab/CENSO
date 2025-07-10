"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""

from collections import defaultdict
from collections.abc import Callable
from pathlib import Path
from tabulate import tabulate
import json

from ..ensembledata import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.parts import NMRConfig
from ..config.job_config import NMRJobConfig
from ..config.parallel_config import ParallelConfig
from ..params import GridLevel
from ..parallel import execute
from ..data import NMRResult
from ..utilities import printf, Factory, h1, h2, timeit, DataDump
from ..logging import setup_logger
from ..processing import QmProc

logger = setup_logger(__name__)


@timeit
def nmr(
    ensemble: EnsembleData, config: PartsConfig, parallel_config: ParallelConfig | None
):
    """
    Calculation of the ensemble NMR of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.
    """
    printf(h2("NMR"))

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc: QmProc = Factory[QmProc].create(config.nmr.prog, "4_NMR")

    # Run NMR calculations
    job_config = NMRJobConfig(
        copy_mo=config.general.copy_mo,
        grid=GridLevel.NMR,
        gas_phase=config.general.gas_phase,
        solvent=config.general.solvent,
        **config.nmr.model_dump(),
    )

    results = execute(
        ensemble.conformers,
        proc.nmr,
        job_config,
        config.nmr.prog,
        "nmr",
        parallel_config,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )

    _write_results(ensemble, config, results)


def _write_results(
    ensemble: EnsembleData, config: PartsConfig, results: dict[str, NMRResult]
):
    boltzmann_populations = ensemble.get_populations(config.general.temperature)

    # Average shielding and coupling values
    shieldings = [
        (i, s * boltzmann_populations[conf.name])
        for conf in ensemble
        for (i, s) in results[conf.name].shieldings
    ]
    averaged_shieldings: defaultdict[int, float] = defaultdict(float)
    for i, v in shieldings:
        averaged_shieldings[i] += v

    couplings = [
        (i, j, c * boltzmann_populations[conf.name])
        for conf in ensemble
        for ((i, j), c) in results[conf.name].couplings
    ]
    averaged_couplings: defaultdict[tuple[int, int], float] = defaultdict(float)
    for i, j, v in couplings:
        averaged_couplings[(i, j)] += v

    filepath = Path("4_NMR.out")
    text = ""
    if len(averaged_shieldings) > 0:
        headers = ["Atom #", "σ"]
        printmap = {"Atom #": lambda k, v: str(k), "σ": lambda k, v: f"{v:.4f}"}

        rows = [
            [printmap[header](k, v) for header in headers]
            for k, v in averaged_shieldings.items()
        ]

        printf(h1("Averaged NMR Shielding Constants"))
        table = tabulate(
            rows,
            headers=headers,
            colalign=["left"] + ["center" for _ in headers[1:]],
            disable_numparse=True,
            numalign="decimal",
        )
        printf(table)
        text += table
    else:
        printf("No shieldings calculated.")

    # Print couplings
    if len(averaged_couplings) > 0:
        headers = ["Atom A #", "Atom B #", "J"]
        units = ["", "", "[Hz]"]
        printmap = {
            "Atom A #": lambda k, v: str(k[0]),
            "Atom B #": lambda k, v: str(k[1]),
            "J": lambda k, v: f"{v:.2f}",
        }

        rows = [
            [printmap[header](k, v) for header in headers]
            for k, v in averaged_couplings.items()
        ]
        for i in range(len(headers)):
            headers[i] += "\n" + units[i]

        print(h1("Averaged NMR Spin-Spin Coupling Constants"))
        table = tabulate(
            rows,
            headers=headers,
            colalign=["left"] + ["center" for _ in headers[1:]],
            disable_numparse=True,
            numalign="decimal",
        )
        printf(table)
        text += table
    else:
        printf("No couplings calculated.")

    filepath.write_text(text)

    # Additionally, write results in json format
    Path("4_NMR.json").write_text(
        json.dumps(jsonify(ensemble, config.nmr, results), indent=4)
    )


def jsonify(ensemble: EnsembleData, config: NMRConfig, results: dict[str, NMRResult]):
    per_conf: Callable[
        [MoleculeData],
        dict[
            str,
            dict[
                str,
                float | list[tuple[int, float]] | list[tuple[tuple[int, int], float]],
            ],
        ],
    ] = lambda conf: {
        conf.name: {
            "energy": conf.energy,
            "gsolv": conf.gsolv,
            "grrho": conf.grrho,
            "nat": conf.geom.nat,
            "shieldings": results[conf.name].shieldings,
            "couplings": results[conf.name].couplings,
        }
    }

    dump = DataDump(part_name="nmr")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
