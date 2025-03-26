from pathlib import Path
from collections.abc import Callable
from typing import Any
import json
from tabulate import tabulate


from ..config.parts.prescreening import PrescreeningConfig
from ..molecules import MoleculeData
from ..ensembledata import EnsembleData
from ..utilities import Factory, timeit, h1, DataDump, printf
from ..config import PartsConfig
from ..parallel import execute
from ..qm import QmProc
from ..params import OMPMIN, NCORES, GridLevel, AU2KCAL, PLENGTH
from ..config.job_config import SPJobConfig, XTBJobConfig
from ..logging import setup_logger

logger = setup_logger(__name__)


@timeit
def prescreening(
    ensemble: EnsembleData,
    config: PartsConfig,
    ncores: int = NCORES or OMPMIN,
    omp: int = OMPMIN,
    cut: bool = True,
):
    """
    This implements a cheap prescreening step using low-cost DFT and possibly
    solvation contributions calculated using xtb.

    The list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).
    """
    # Setup processor and target
    proc: QmProc = Factory[QmProc].create(config.prescreening.prog, "0_PRESCREENING")

    if not config.general.gas_phase:
        # Calculate Gsolv using xtb
        job_config = XTBJobConfig(
            gfnv=config.prescreening.gfnv,
            solvent=config.general.solvent,
            sm_rrho=config.general.sm_rrho,
            gas_phase=False,
        )
        results, _ = execute(
            ensemble.conformers,
            proc.xtb_gsolv,
            job_config,
            config.prescreening.prog,
            ncores,
            omp,
            balance=config.general.balance,
            copy_mo=config.general.copy_mo,
        )

        for conf in ensemble:
            conf.gsolv = results[conf.name].gsolv

    # Calculate gas-phase single-point
    job_config = SPJobConfig(
        copy_mo=config.general.copy_mo,
        func=config.prescreening.func,
        basis=config.prescreening.basis,
        grid=GridLevel.LOW,
        template=config.prescreening.template,
        gas_phase=True,
    )
    results, _ = execute(
        ensemble.conformers,
        proc.sp,
        job_config,
        config.prescreening.prog,
        ncores,
        omp,
        balance=config.general.balance,
        copy_mo=config.general.copy_mo,
    )
    for conf in ensemble:
        conf.energy = results[conf.name].energy

    if cut:
        threshold = (
            min(conf.gtot for conf in ensemble)
            + config.prescreening.threshold / AU2KCAL
        )
        ensemble.remove_conformers(cond=lambda conf: conf.gtot > threshold)

    ensemble.set_populations(config.general.temperature)

    # print results
    _write_results(ensemble, config)


def _write_results(ensemble: EnsembleData, config: PartsConfig) -> None:
    """ """
    printf(h1(f"PRESCREENING SINGLE-POINT RESULTS"))

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
        "Boltzmann weight": lambda conf: f"{conf.bmw * 100:.2f}",
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
    avG = sum([conf.bmw * conf.gtot for conf in ensemble])

    # calculate averaged free energy (?)
    avE = sum([conf.bmw * conf.energy for conf in ensemble])

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
    Path("0_PRESCREENING.json").write_text(
        json.dumps(jsonify(ensemble, config.prescreening), indent=4)
    )

    ensemble.dump_xyz(Path("0_PRESCREENING.xyz"))


# TODO: generalize this
def jsonify(
    ensemble: EnsembleData,
    config: PrescreeningConfig,
    fields: Callable[[MoleculeData], dict[str, Any]] | None = None,
):
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
