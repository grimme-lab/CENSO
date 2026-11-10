"""
Per-temperature output of the ensemble optimization parts (multitemp).
"""

from typing import Any

from ..ensemble import EnsembleData


def temperature_block(
    ensemble: EnsembleData,
    temperatures: list[float],
    main_temperature: float,
    tag: str,
    results_dict: dict[str, Any],
) -> list[str]:
    """
    Boltzmann averaged contributions of the ensemble for every temperature.
    Also stores the averages and populations per temperature in results_dict ("averages_t", "populations_t").

    :param ensemble: The conformer ensemble.
    :type ensemble: EnsembleData
    :param temperatures: Temperatures in K (ascending).
    :type temperatures: list[float]
    :param main_temperature: Main temperature in K, its row is marked with tag.
    :type main_temperature: float
    :param tag: Part tag used to mark the main temperature row (e.g. "part1").
    :type tag: str
    :param results_dict: JSON results of the part, updated in place.
    :type results_dict: dict[str, Any]
    :return: Header and one row per temperature.
    :rtype: list[str]
    """
    lines = [
        f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avGmRRHO(T) /a.u.':>18} {'avGsolv(T) /a.u.':>17} {'avG(T) /a.u.':>14}"
    ]
    averages_t: dict[str, dict[str, float]] = {}
    populations_t: dict[str, dict[str, float]] = {}
    for t in temperatures:
        av = ensemble.boltzmann_averages(t)
        row = f"{t:^15} {av['avE']:>14.7f} {av['avGmRRHO']:>18.7f} {av['avGsolv']:>17.7f} {av['avG']:>14.7f}"
        if t == round(main_temperature, 2):
            row += f"     <<=={tag}=="
        lines.append(row)
        averages_t[str(t)] = av
        populations_t[str(t)] = ensemble.get_populations(t)

    results_dict["averages_t"] = averages_t
    results_dict["populations_t"] = populations_t
    return lines
