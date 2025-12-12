"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""

from collections import defaultdict
from pathlib import Path
from tabulate import tabulate
from itertools import product
from typing import Any
import json
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..molecules import MoleculeData
from ..config import PartsConfig
from ..config.parts import NMRConfig
from ..config.job_config import NMRJobConfig
from ..parallel import execute
from ..processing.results import NMRResult
from ..utilities import printf, Factory, h1, h2, DataDump
from ..logging import setup_logger
from ..processing import QmProc
from ..params import GridLevel, PLENGTH

logger = setup_logger(__name__)


def nmr(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
) -> dict[str, Any]:
    """
    Calculation of the ensemble NMR of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :return: JSON-serializable dictionary of NMR results.
    """
    printf(h2("NMR"))

    config = PartsConfig.model_validate(config, context={"check": ["nmr"]})

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    # Setup processor and target
    proc = Factory[QmProc].create(config.nmr.prog, "4_NMR")

    # Run NMR calculations
    job_config = NMRJobConfig(
        copy_mo=config.general.copy_mo,
        grid=GridLevel.NMR,
        gas_phase=config.general.gas_phase,
        solvent=config.general.solvent,
        paths=config.paths,
        **config.nmr.model_dump(),
    )

    results = execute(
        ensemble.conformers,
        proc.nmr,
        job_config,
        config.nmr.prog,
        "nmr",
        client,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )
    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in results)

    results_dict = jsonify(ensemble, config.nmr, results)
    _write_results(ensemble, config, results, results_dict)

    return results_dict


def read_chemeq() -> dict[int, list[int]]:
    """
    Read chemical equivalent nuclei from anmr_nucinfo.
    NOTE: crest starts counting from 1. In CENSO, we're counting from 0.

    :return: dict[int, list[int]]
    """
    lines = Path("anmr_nucinfo").read_text().split("\n")
    atoms = []
    equiv: dict[int, list[int]] = {}
    for line in lines[1::2]:
        atoms.append(int(line.split()[0]))

    for i, line in enumerate(lines[2::2]):
        equiv[atoms[i]] = [int(x) for x in line.split()]

    return equiv


def _write_results(
    ensemble: EnsembleData,
    config: PartsConfig,
    results: dict[str, NMRResult],
    results_dict: dict[str, Any],
):
    boltzmann_populations = ensemble.get_populations(config.general.temperature)

    try:
        equiv = read_chemeq()
        printf("Applying chemical equivalence as predicted by CREST.")
    except FileNotFoundError:
        printf("Could not find anmr_nucinfo file. Did you run crest with -nmr?")
        # In this case set every atom equivalent to only itself
        conf = next(iter(ensemble.conformers))
        equiv = {i: [i] for (i, _) in results[conf.name].shieldings}

    # Average shielding and coupling values
    shieldings = [
        (i, s * boltzmann_populations[conf.name])
        for conf in ensemble
        for (i, s) in results[conf.name].shieldings
    ]
    averaged_shieldings: defaultdict[int, float] = defaultdict(float)
    for i, v in shieldings:
        # Apply chemical equivalence
        for j in equiv[i]:
            averaged_shieldings[j] += v / len(equiv[i])

    couplings = [
        (i, j, c * boltzmann_populations[conf.name])
        for conf in ensemble
        for ((i, j), c) in results[conf.name].couplings
    ]
    averaged_couplings: defaultdict[tuple[int, int], float] = defaultdict(float)
    for i, j, v in couplings:
        # Apply chemical equivalence
        equiv_i = equiv[i]
        equiv_j = equiv[j]
        n_pairs = len(equiv_i) * len(equiv_j)
        for ieq, jeq in product(equiv_i, equiv_j):
            averaged_couplings[(ieq, jeq)] += v / n_pairs

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

        printf(h1("Averaged NMR Spin-Spin Coupling Constants"))
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

    printf("".ljust(int(PLENGTH), "-"))
    filepath.write_text(text)

    # Additionally, write results in json format
    Path("4_NMR.json").write_text(json.dumps(results_dict, indent=4))


def jsonify(ensemble: EnsembleData, config: NMRConfig, results: dict[str, NMRResult]):
    """
    Prepare NMR results for JSON serialization.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: NMRConfig object with NMR configuration settings.
    :param results: Dictionary of NMRResult objects for each conformer.
    :return: Dictionary ready for JSON serialization.
    """

    def per_conf(
        conf: MoleculeData,
    ) -> dict[
        str,
        dict[
            str, float | list[tuple[int, float]] | list[tuple[tuple[int, int], float]]
        ],
    ]:
        return {
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
