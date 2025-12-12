from pathlib import Path
from collections.abc import Callable
from typing import Any
import json
from tabulate import tabulate
from dask.distributed import Client

from ..processing.xtb_processor import XtbProc
from ..config.parts.prescreening import PrescreeningConfig
from ..molecules import Contributions, MoleculeData
from ..ensemble import EnsembleData
from ..utilities import Factory, h1, h2, DataDump, printf
from ..config import PartsConfig
from ..parallel import execute
from ..processing import QmProc
from ..params import GridLevel, AU2KCAL, PLENGTH, Prog
from ..config.job_config import SPJobConfig, XTBJobConfig
from ..logging import setup_logger

logger = setup_logger(__name__)


def prescreening(
    ensemble: EnsembleData,
    config: PartsConfig,
    client: Client,
    cut: bool = True,
) -> dict[str, Any]:
    """
    This implements a cheap prescreening step using low-cost DFT and possibly
    solvation contributions calculated using xtb.

    The list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).

    :param ensemble: EnsembleData object containing the conformers.
    :type ensemble: EnsembleData
    :param config: PartsConfig object with configuration settings.
    :type config: PartsConfig
    :param client: dask.distributed.Client for parallel execution.
    :type client: Client
    :param cut: Whether to apply cutting conditions.
    :type cut: bool
    :return: JSON-serializable dictionary of prescreening results.
    :rtype: dict
    """
    printf(h2("PRESCREENING"))

    config = PartsConfig.model_validate(config, context={"check": ["prescreening"]})

    # Setup processor and target
    proc = Factory[QmProc].create(config.prescreening.prog, "0_PRESCREENING")

    contributions_dict = {conf.name: Contributions() for conf in ensemble}
    if not config.general.gas_phase:
        # Calculate Gsolv using xtb
        xtb_proc = Factory[XtbProc].create(Prog.XTB, "0_PRESCREENING")
        gsolv_job_config = XTBJobConfig(
            gfnv=config.prescreening.gfnv,
            solvent=config.general.solvent,
            sm_rrho=config.general.sm_rrho,
            gas_phase=False,
            paths=config.paths,
        )
        gsolv_results = execute(
            ensemble.conformers,
            xtb_proc.gsolv,
            gsolv_job_config,
            "xtb",
            "prescreening",
            client,
            ignore_failed=config.general.ignore_failed,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        if config.general.ignore_failed:
            ensemble.remove_conformers(lambda conf: conf.name not in gsolv_results)

        for conf in ensemble:
            contributions_dict[conf.name].gsolv = gsolv_results[conf.name].gsolv

    # Calculate gas-phase single-point
    sp_job_config = SPJobConfig(
        copy_mo=config.general.copy_mo,
        func=config.prescreening.func,
        basis=config.prescreening.basis,
        grid=GridLevel.LOW,
        template=config.prescreening.template,
        gas_phase=True,
        paths=config.paths,
    )
    sp_results = execute(
        ensemble.conformers,
        proc.sp,
        sp_job_config,
        config.prescreening.prog,
        "prescreening",
        client,
        ignore_failed=config.general.ignore_failed,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )

    if config.general.ignore_failed:
        ensemble.remove_conformers(lambda conf: conf.name not in sp_results)

    for conf in ensemble:
        contributions_dict[conf.name].energy = sp_results[conf.name].energy

    # Update molecules
    ensemble.update_contributions(contributions_dict)

    if cut:
        threshold = (
            min(conf.gtot for conf in ensemble)
            + config.prescreening.threshold / AU2KCAL
        )
        ensemble.remove_conformers(cond=lambda conf: conf.gtot > threshold)

    # print results
    results_dict = jsonify(ensemble, config.prescreening)
    _write_results(ensemble, config, results_dict)

    return results_dict


def _write_results(
    ensemble: EnsembleData, config: PartsConfig, results_dict: dict[str, Any]
) -> None:
    """ """
    printf(h1("PRESCREENING SINGLE-POINT RESULTS"))

    # column headers
    headers = [
        "CONF#",
        # "E (xTB)",
        # "ΔE (xTB)",
        "E (DFT)",
        "ΔE (DFT)",
        "ΔGsolv (xTB)",
        # "δΔGsolv",
        "Gtot",
        "ΔGtot",
        "Boltzmann weight",
    ]

    # column units
    units = [
        "",
        # "[Eh]",
        # "[kcal/mol]",
        "[Eh]",
        "[kcal/mol]",
        "[kcal/mol]",
        # "[kcal/mol]",
        "[Eh]",
        "[kcal/mol]",
        f"% at {config.general.temperature} K",
    ]

    # variables for printmap
    # minimal xtb single-point energy
    # if all(
    #     "xtb_gsolv" in self.data["results"][conf.name]
    #     for conf in self._ensemble.conformers
    # ):
    #     xtbmin = min(
    #         self.data["results"][conf.name]["xtb_gsolv"]["energy_xtb_gas"]
    #         for conf in self._ensemble.conformers
    #     )

    # minimal dft single-point energy
    dftmin = min(conf.energy for conf in ensemble)

    # minimal solvation free enthalpy
    # gsolvmin = min(conf.gsolv for conf in ensemble)

    # minimal total free enthalpy
    gtotmin = min(conf.gtot for conf in ensemble)

    boltzmann_populations = ensemble.get_populations(config.general.temperature)

    # determines what to print for each conformer in each column
    printmap = {
        "CONF#": lambda conf: conf.name,
        # "E (xTB)": lambda conf: (
        #     f"{self.data['results'][conf.name]['xtb_gsolv']['energy_xtb_gas']:.6f}"
        #     if "xtb_gsolv" in self.data["results"][conf.name]
        #     else "---"
        # ),
        # "ΔE (xTB)": lambda conf: (
        #     f"{(self.data['results'][conf.name]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}"
        #     if "xtb_gsolv" in self.data["results"][conf.name]
        #     else "---"
        # ),
        "E (DFT)": lambda conf: f"{conf.energy:.6f}",
        "ΔE (DFT)": lambda conf: f"{(conf.energy - dftmin) * AU2KCAL:.2f}",
        "ΔGsolv (xTB)": lambda conf: (f"{conf.gsolv * AU2KCAL:.6f}"),
        "Gtot": lambda conf: f"{conf.gtot:.6f}",
        # "δΔGsolv": lambda conf: f"{(self.data["results"][conf.name]['xtb_gsolv']['gsolv'] - gsolvmin) * AU2KCAL:.2f}"
        # if "xtb_gsolv" in self.data["results"][conf.name].keys()
        # else "---",
        "ΔGtot": lambda conf: f"{(conf.gtot - gtotmin) * AU2KCAL:.2f}",
        "Boltzmann weight": lambda conf: f"{boltzmann_populations[conf.name] * 100:.2f}",
    }

    rows = [[printmap[header](conf) for header in headers] for conf in ensemble]

    for i in range(len(headers)):
        headers[i] += "\n" + units[i]

    table = tabulate(
        rows,
        headers=headers,
        colalign=["left"] + ["center" for _ in headers[1:]],
        disable_numparse=True,
        numalign="decimal",
    )
    print(table, flush=True)

    # list the averaged free enthalpy of the ensemble
    lines = []
    lines.append(
        "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):"
    )
    lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}")

    # calculate averaged free enthalpy
    avG = sum([boltzmann_populations[conf.name] * conf.gtot for conf in ensemble])

    # calculate averaged free energy (?)
    avE = sum([boltzmann_populations[conf.name] * conf.energy for conf in ensemble])

    # append the lines for the free energy/enthalpy
    lines.append(
        f"{config.general.temperature:^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0=="
    )
    lines.append("".ljust(int(PLENGTH), "-"))

    # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")

    # Print everything
    for line in lines:
        printf(line)

    # write everything to a file
    filepath = Path("0_PRESCREENING.out")
    filepath.write_text(table + "\n".join(lines))

    # Additionally, write results in json format
    Path("0_PRESCREENING.json").write_text(json.dumps(results_dict, indent=4))

    ensemble.dump_xyz(Path("0_PRESCREENING.xyz"))


def jsonify(
    ensemble: EnsembleData,
    config: PrescreeningConfig,
    fields: Callable[[MoleculeData], dict[str, Any]] | None = None,
):
    """
    Convert ensemble data to JSON format for prescreening results.

    :param ensemble: EnsembleData object.
    :type ensemble: EnsembleData
    :param config: PrescreeningConfig object.
    :type config: PrescreeningConfig
    :param fields: Optional callable to customize fields.
    :type fields: Callable[[MoleculeData], dict[str, Any]] | None
    :return: JSON-serializable dictionary.
    :rtype: dict[str, Any]
    """
    per_conf: Callable[[MoleculeData], dict[str, dict[str, float]]] = fields or (
        lambda conf: {
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
                "gtot": conf.gtot,
            }
        }
    )

    dump = DataDump(part_name="prescreening")

    for conf in ensemble:
        dump.data.update(per_conf(conf))

    dump.settings = config.model_dump()

    return dump.model_dump(mode="json")
