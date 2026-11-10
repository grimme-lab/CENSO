"""
Result output shared by all CENSO parts.

A part is identified by its prefix '<n>_<PART>' (e.g. '0_PRESCREENING', '4_NMR'),
which is also its working directory name. The lower-case '<part>' is the name of the
part's settings field in PartsConfig. A part's output consists of:

- the result report, printed and written verbatim to '<prefix>.out',
- '<prefix>.json', a DataDump with per-conformer energies and part-specific fields,
- '<prefix>.xyz', the remaining ensemble (ensemble optimization parts only).
"""

import json
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any

from pydantic import BaseModel
from tabulate import tabulate

from .config import PartsConfig
from .ensemble import EnsembleData
from .molecules import MoleculeData
from .params import AU2KCAL, PLENGTH
from .utilities import h1, printf

type Column[T] = tuple[str, str, Callable[[T], str]]
"""Table column: header, unit and formatter for one item."""


class DataDump(BaseModel):
    """
    Model for dumping data with part name and settings.
    """

    part_name: str
    data: dict[str, dict[str, Any]] = {}
    settings: dict[str, Any] = {}


def format_table[T](columns: Sequence[Column[T]], items: Iterable[T]) -> str:
    """
    Render one table row per item.

    :param columns: Header, unit and formatter per column. Units are printed below the
        headers unless all units are empty.
    :param items: Objects to format, one per row.
    :return: The table; first column left-aligned, all others centered.
    """
    multiline = any(unit for _, unit, _ in columns)
    return tabulate(
        [[fmt(item) for _, _, fmt in columns] for item in items],
        headers=[
            f"{header}\n{unit}" if multiline else header for header, unit, _ in columns
        ],
        colalign=["left"] + ["center"] * (len(columns) - 1),
        disable_numparse=True,
        numalign="decimal",
    )


def write_output(
    prefix: str,
    ensemble: EnsembleData,
    config: PartsConfig,
    report: str,
    fields: Callable[[MoleculeData], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """
    Print the result report of a part and write '<prefix>.out' and '<prefix>.json'.

    :param prefix: Part prefix '<n>_<PART>'.
    :param ensemble: Conformers to dump.
    :param config: Settings; the part's own settings are stored in the json dump.
    :param report: Result report; a closing separator line is appended.
    :param fields: Part-specific json fields per conformer, added to energy, gsolv, grrho and gtot.
    :return: JSON-serializable dump, identical to the content of '<prefix>.json'.
    """
    part = prefix.split("_", 1)[1].lower()
    text = f"{report}\n{'-' * PLENGTH}"
    printf(text)
    Path(f"{prefix}.out").write_text(text + "\n")

    dump = DataDump(
        part_name=part,
        data={
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
                "gtot": conf.gtot,
            }
            | (fields(conf) if fields else {})
            for conf in ensemble
        },
        settings=getattr(config, part).model_dump(),
    ).model_dump(mode="json")
    Path(f"{prefix}.json").write_text(json.dumps(dump, indent=4))
    return dump


def write_ranking(
    prefix: str,
    ensemble: EnsembleData,
    config: PartsConfig,
    title: str,
    columns: Sequence[str] = (
        "CONF#",
        "E (DFT)",
        "ΔGsolv",
        "GmRRHO",
        "Gtot",
        "ΔGtot",
        "Boltzmann weight",
    ),
    geometries: str | None = None,
) -> dict[str, Any]:
    """
    Write the output of an ensemble optimization part: ranking table with Boltzmann
    weights and averages (printed, '<prefix>.out', '<prefix>.json') and the remaining
    ensemble as '<prefix>.xyz'.

    :param prefix: Part prefix '<n>_<PART>'.
    :param ensemble: Ranked conformers.
    :param config: Settings; the temperature is taken from config.general.
    :param title: Heading of the ranking table.
    :param columns: Table columns, any of: CONF#, E (DFT), ΔE (DFT), ΔGsolv,
        ΔGsolv (xTB), GmRRHO, Gtot, ΔGtot, Boltzmann weight.
    :param geometries: Kind of geometries the averages refer to, e.g. 'DFT optimized'.
    :return: JSON-serializable dump, identical to the content of '<prefix>.json'.
    """
    temperature = config.general.temperature
    populations = ensemble.get_populations(temperature)
    emin = min(conf.energy for conf in ensemble)
    gtotmin = min(conf.gtot for conf in ensemble)
    catalog: dict[str, tuple[str, Callable[[MoleculeData], str]]] = {
        "CONF#": ("", lambda conf: conf.name),
        "E (DFT)": ("[Eh]", lambda conf: f"{conf.energy:.6f}"),
        "ΔE (DFT)": (
            "[kcal/mol]",
            lambda conf: f"{(conf.energy - emin) * AU2KCAL:.2f}",
        ),
        "ΔGsolv": ("[kcal/mol]", lambda conf: f"{conf.gsolv * AU2KCAL:.6f}"),
        "ΔGsolv (xTB)": ("[kcal/mol]", lambda conf: f"{conf.gsolv * AU2KCAL:.6f}"),
        "GmRRHO": ("[Eh]", lambda conf: f"{conf.grrho:.6f}"),
        "Gtot": ("[Eh]", lambda conf: f"{conf.gtot:.6f}"),
        "ΔGtot": (
            "[kcal/mol]",
            lambda conf: f"{(conf.gtot - gtotmin) * AU2KCAL:.2f}",
        ),
        "Boltzmann weight": (
            f"% at {temperature} K",
            lambda conf: f"{populations[conf.name] * 100:.2f}",
        ),
    }
    av_e = sum(populations[conf.name] * conf.energy for conf in ensemble)
    av_g = sum(populations[conf.name] * conf.gtot for conf in ensemble)
    on = f" on input geometries ({geometries})" if geometries else ""
    report = "\n".join(
        [
            h1(title),
            format_table([(name, *catalog[name]) for name in columns], ensemble),
            f"\nBoltzmann averaged free energy/enthalpy of ensemble{on}:",
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}",
            (
                f"{temperature:^15} {av_e:>14.7f}  {av_g:>14.7f}"
                f"     <<==part{prefix.split('_')[0]}=="
            ),
        ]
    )
    dump = write_output(prefix, ensemble, config, report)
    ensemble.dump_xyz(Path(f"{prefix}.xyz"))
    return dump
