"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""

from collections import defaultdict
from pathlib import Path
from tabulate import tabulate
from itertools import product
from typing import Any, Literal
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
    retemp: Literal["xtb", "cosmors"] | None = None,
) -> dict[str, Any]:
    """
    Calculation of the ensemble NMR of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param client: dask.distributed.Client for parallel execution.
    :param retemp: If set, reuse the reloaded NMR results (no calculations) and only redo the Boltzmann averaging.
        Results are written to 4_NMR_RETEMP.out/.json.
    :return: JSON-serializable dictionary of NMR results.
    """
    printf(h2("NMR"))

    config = PartsConfig.model_validate(config.model_dump(), context={"check": ["nmr"]})

    # Assert that all conformers have energies defined
    # TODO: this is not optimal since == 0 does not mean that no ensembleopt has been performed before
    if not all(conf.energy != 0 for conf in ensemble):
        raise RuntimeError(
            "Before calculating an ensemble property one has to run at least one ensemble refinement step (prescreening, screening, optimization or refinement)."
        )

    if retemp is not None:
        stored = ensemble.property_results.get("nmr")
        if stored is None:
            raise RuntimeError(
                "--retemp for nmr requires the previous results: add 4_NMR.json to --reload."
            )
        results = {
            conf.name: NMRResult(
                shieldings=[(i, s) for i, s in stored[conf.name]["shieldings"]],
                couplings=[((i, j), c) for (i, j), c in stored[conf.name]["couplings"]],
            )
            for conf in ensemble
        }
    else:
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
    _write_results(ensemble, config, results, results_dict, "_RETEMP" if retemp else "")

    return results_dict


def read_chemeq() -> dict[int, list[int]]:
    """
    Read chemical equivalent nuclei from anmr_nucinfo.
    NOTE: crest starts counting from 1. In CENSO, we're counting from 0.

    :return: dict[int, list[int]]
    """
    lines = Path("anmr_nucinfo").read_text().split("\n")
    # Filter out empty lines (e.g., from trailing newlines)
    lines = [line for line in lines if line.strip()]
    atoms = []
    equiv: dict[int, list[int]] = {}
    for line in lines[1::2]:
        atoms.append(int(line.split()[0]) - 1)  # Convert from 1-indexed to 0-indexed

    for i, line in enumerate(lines[2::2]):
        equiv[atoms[i]] = [
            int(x) - 1 for x in line.split()
        ]  # Convert from 1-indexed to 0-indexed

    return equiv


def _average(
    ensemble: EnsembleData,
    results: dict[str, NMRResult],
    populations: dict[str, float],
    equiv: dict[int, list[int]],
) -> tuple[dict[int, float], dict[tuple[int, int], float]]:
    """
    Boltzmann average shieldings and couplings, applying chemical equivalence.

    :param ensemble: EnsembleData object containing the conformers.
    :param results: NMR results per conformer.
    :param populations: Boltzmann populations per conformer.
    :param equiv: Chemically equivalent atoms per atom.
    :return: Averaged shieldings per atom and couplings per atom pair.
    """
    shieldings = [
        (i, s * populations[conf.name])
        for conf in ensemble
        for (i, s) in results[conf.name].shieldings
    ]
    averaged_shieldings: defaultdict[int, float] = defaultdict(float)
    for i, v in shieldings:
        # Apply chemical equivalence
        for j in equiv[i]:
            averaged_shieldings[j] += v / len(equiv[i])

    couplings = [
        (i, j, c * populations[conf.name])
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

    return dict(averaged_shieldings), dict(averaged_couplings)


def _tables(
    shieldings: dict[int, float], couplings: dict[tuple[int, int], float]
) -> str:
    """
    Format averaged shieldings and couplings as tables.

    :param shieldings: Averaged shieldings per atom.
    :param couplings: Averaged couplings per atom pair.
    :return: Formatted text.
    """
    parts = []
    if len(shieldings) > 0:
        headers = ["Atom #", "σ"]
        rows = [[str(k), f"{v:.4f}"] for k, v in shieldings.items()]
        parts.append(h1("Averaged NMR Shielding Constants"))
        parts.append(
            tabulate(
                rows,
                headers=headers,
                colalign=["left"] + ["center" for _ in headers[1:]],
                disable_numparse=True,
                numalign="decimal",
            )
        )
    else:
        parts.append("No shieldings calculated.")

    if len(couplings) > 0:
        headers = ["Atom A #\n", "Atom B #\n", "J\n[Hz]"]
        rows = [[str(k[0]), str(k[1]), f"{v:.2f}"] for k, v in couplings.items()]
        parts.append(h1("Averaged NMR Spin-Spin Coupling Constants"))
        parts.append(
            tabulate(
                rows,
                headers=headers,
                colalign=["left"] + ["center" for _ in headers[1:]],
                disable_numparse=True,
                numalign="decimal",
            )
        )
    else:
        parts.append("No couplings calculated.")

    return "\n".join(parts)


def _write_results(
    ensemble: EnsembleData,
    config: PartsConfig,
    results: dict[str, NMRResult],
    results_dict: dict[str, Any],
    suffix: str,
):
    """
    Print the averaged NMR data and write it to 4_NMR<suffix>.out/.json, including per-temperature averages with multitemp.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: PartsConfig object with configuration settings.
    :param results: NMR results per conformer.
    :param results_dict: JSON results, per-temperature averages and populations are added with multitemp.
    :param suffix: File name suffix ("_RETEMP" for --retemp runs).
    """
    try:
        equiv = read_chemeq()
        printf("Applying chemical equivalence as predicted by CREST.")
    except FileNotFoundError:
        printf("Could not find anmr_nucinfo file. Did you run crest with -nmr?")
        # In this case set every atom equivalent to only itself
        conf = next(iter(ensemble.conformers))
        equiv = {i: [i] for (i, _) in results[conf.name].shieldings}

    text = _tables(
        *_average(
            ensemble,
            results,
            ensemble.get_populations(config.general.temperature),
            equiv,
        )
    )
    printf(text)
    printf("".ljust(int(PLENGTH), "-"))

    temperatures = ensemble.temperature_range(config.general, "nmr")
    if temperatures:
        results_dict["populations_t"] = {}
        results_dict["averages_t"] = {}
    for t in temperatures:
        populations = ensemble.get_populations(t)
        shieldings, couplings = _average(ensemble, results, populations, equiv)
        text += f"\n\nT = {t} K\n" + _tables(shieldings, couplings)
        results_dict["populations_t"][str(t)] = populations
        results_dict["averages_t"][str(t)] = {
            "shieldings": [[i, s] for i, s in shieldings.items()],
            "couplings": [[[i, j], c] for (i, j), c in couplings.items()],
        }

    Path(f"4_NMR{suffix}.out").write_text(text)

    # Additionally, write results in json format
    Path(f"4_NMR{suffix}.json").write_text(json.dumps(results_dict, indent=4))


def jsonify(ensemble: EnsembleData, config: NMRConfig, results: dict[str, NMRResult]):
    """
    Prepare NMR results for JSON serialization.

    :param ensemble: EnsembleData object containing the conformers.
    :param config: NMRConfig object with NMR configuration settings.
    :param results: Dictionary of NMRResult objects for each conformer.
    :return: Dictionary ready for JSON serialization.
    """

    def per_conf(conf: MoleculeData) -> dict[str, dict[str, Any]]:
        return {
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
                **conf.temperature_data(),
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
