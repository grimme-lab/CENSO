"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""

from collections import defaultdict
from pathlib import Path
from itertools import product
from typing import Any
from dask.distributed import Client

from ..ensemble import EnsembleData
from ..config import PartsConfig
from ..config.job_config import NMRJobConfig
from ..parallel import execute
from ..processing.results import NMRResult
from ..output import Column, format_table, write_output
from ..utilities import printf, Factory, h1, h2
from ..logging import setup_logger
from ..processing import QmProc
from ..params import GridLevel

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

    config = PartsConfig.model_validate(config.model_dump(), context={"check": ["nmr"]})

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

    return write_output(
        "4_NMR",
        ensemble,
        config,
        _report(ensemble, config, results),
        lambda conf: {
            "nat": conf.geom.nat,
            "shieldings": results[conf.name].shieldings,
            "couplings": results[conf.name].couplings,
        },
    )


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


def _report(
    ensemble: EnsembleData, config: PartsConfig, results: dict[str, NMRResult]
) -> str:
    """
    Boltzmann-average shieldings and couplings, applying chemical equivalence from
    anmr_nucinfo, and format them as the NMR result report.

    :param ensemble: Conformers to average over.
    :param config: Settings; the temperature is taken from config.general.
    :param results: NMR results per conformer name.
    :return: The report.
    """
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

    blocks: list[str] = []
    if averaged_shieldings:
        shielding_columns: list[Column[tuple[int, float]]] = [
            ("Atom #", "", lambda item: str(item[0])),
            ("σ", "", lambda item: f"{item[1]:.4f}"),
        ]
        blocks += [
            h1("Averaged NMR Shielding Constants"),
            format_table(shielding_columns, averaged_shieldings.items()),
        ]
    else:
        blocks.append("No shieldings calculated.")
    if averaged_couplings:
        coupling_columns: list[Column[tuple[tuple[int, int], float]]] = [
            ("Atom A #", "", lambda item: str(item[0][0])),
            ("Atom B #", "", lambda item: str(item[0][1])),
            ("J", "[Hz]", lambda item: f"{item[1]:.2f}"),
        ]
        blocks += [
            h1("Averaged NMR Spin-Spin Coupling Constants"),
            format_table(coupling_columns, averaged_couplings.items()),
        ]
    else:
        blocks.append("No couplings calculated.")
    return "\n".join(blocks)
