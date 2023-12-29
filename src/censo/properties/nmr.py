"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""
from functools import reduce
import os

from ..ensembledata import EnsembleData
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    SOLVENTS_DB,
)
from ..part import CensoPart
from ..utilities import print, timeit, DfaHelper, setup_logger, format_data

logger = setup_logger(__name__)


class NMR(CensoPart):
    alt_name = "part2"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "resonance_frequency": {"default": 300.0, "range": [150.0, 1000.0]},
        "threshold_bmw": {"default": 0.95, "range": [0.01, 0.99]},
        "prog": {"default": "orca", "options": PROGS},  # required
        "func_j": {"default": "pbe0-d4", "options": DfaHelper.find_func("nmr_j")},
        "basis_j": {"default": "def2-TZVP", "options": BASIS_SETS},
        "sm_j": {"default": "smd", "options": __solv_mods},
        "func_s": {"default": "pbe0-d4", "options": DfaHelper.find_func("nmr_s")},
        "basis_s": {"default": "def2-TZVP", "options": BASIS_SETS},
        "sm_s": {"default": "smd", "options": __solv_mods},
        "h_ref": {"default": "TMS", "options": ["TMS"]},
        "c_ref": {"default": "TMS", "options": ["TMS"]},
        "f_ref": {"default": "CFCl3", "options": ["CFCl3"]},
        "si_ref": {"default": "TMS", "options": ["TMS"]},
        "p_ref": {"default": "TMP", "options": ["TMP", "PH3"]},
        "grid": {"default": "high+", "options": GRIDOPTIONS},  # required
        "run": {"default": False},  # required
        "template": {"default": False},  # required
        "gcp": {"default": True},  # required
        "couplings": {"default": True},
        "shieldings": {"default": True},
        "h_active": {"default": True},
        "c_active": {"default": True},
        "f_active": {"default": False},
        "si_active": {"default": False},
        "p_active": {"default": False},
    }

    _settings = {}

    def __init__(self, censo: EnsembleData):
        super().__init__(censo)

        # Set the correct name for 'func_s' and 'func_j', as well as the solvent key
        for c in ["s", "j"]:
            self._instructions[f"func_type_{c}"] = DfaHelper.get_type(
                self._instructions[f"func_{c}"]
            )
            self._instructions[f"func_name_{c}"] = DfaHelper.get_name(
                self._instructions[f"func_{c}"], self._instructions["prog"]
            )
            self._instructions[f"disp_{c}"] = DfaHelper.get_disp(
                self._instructions[f"func_{c}"]
            )
            self._instructions[f"solvent_key_prog_{c}"] = SOLVENTS_DB.get(
                self._instructions["solvent"]
            )[self._instructions[f"sm_{c}"]][1]

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        Calculation of the ensemble NMR of a (previously) optimized ensemble.
        """

        # print instructions
        self.print_info()

        self._instructions["jobtype"] = ["nmr"]
        if not self._instructions["couplings"] and not self._instructions["shieldings"]:
            # This case should basically never happen except for user input error
            raise (
                RuntimeError(
                    "No jobtype selected. "
                    "Please select at least one of the following: couplings, shieldings"
                )
            )

        # Select conformers based on Boltzmann weight threshold, index -1 indicates to always use the most recently
        # calculated Boltzmann weight
        self.ensemble.update_conformers(
            lambda conf: conf.bmws[-1],
            self._instructions["threshold_bmw"],
            boltzmann=True,
        )

        # Store the utilized Boltzmann population in order to have it in the resulting json file
        for conf in self.ensemble.conformers:
            conf.results.setdefault(self._name, {})["bmw"] = conf.bmws[-1]

        # Execute jobs in parallel
        results = execute(self.ensemble.conformers,
                          self._instructions, self.dir)

        # Put results into conformers
        for conf in self.ensemble.conformers:
            # store results
            conf.results.setdefault(self._name, {}).update(
                results[conf.geom.id])

        # Calculate Boltzmann weighted parameters and generate files for ANMR
        self.__generate_anmr()

        # Write data
        self.write_results()

    def __generate_anmr(self):
        """
        Generate all necessary files for an ANMR run.
        """
        # Write 'anmr_enso'
        headers = [
            "ONOFF",
            "NMR",
            "CONF",
            "BW",
            "Energy",
            "Gsolv",
            "mRRHO",
            "gi",
        ]

        rows = []

        lines = format_data(headers, rows)

        # Write lines to file
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, 'anmr_enso')}."
        )
        with open(
            os.path.join(self.ensemble.workdir, "anmr_enso"), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results to a json file
        self.write_json()

        # Write 'anmrrc'
        lines = []

        general_settings = self.get_general_settings()
        lines.append("7 8 XH acid atoms\n")
        lines.append(
            f"ENSO qm= {self._settigs['prog'].upper()} "
            f"mf= {self._settings['resonance_frequency']:.2f} "
            f"lw= 1.0 J= {'on' if self._settings['couplings'] else 'off'} "
            f"S= {'on' if self._settings['shieldings'] else 'off'} "
            f"T= {general_settings['temperature']:6.2f}\n"
        )

        # TODO - since the only program right now is ORCA, the reference solvent model is always SMD
        # Check, if a geometry optimization has been done before
        if all("optimization" in conf.results.keys() for conf in self.ensemble.conformers):
            from ..ensembleopt.optimization import Optimization
            func_geomopt = Optimization.get_settings()["func"]
            basis_geomopt = Optimization.get_setting()["basis"]
        # Otherwise, warn the user
        else:
            logger.warning("Geometries of conformers have not been optimized using CENSO. "
                           "This is advised for accurate results. Also, user configuration"
                           " of .anmrrc is required. Insert functional and basis used for "
                           "reference geometry.")
            func_geomopt = "GEOMOPT_FUNC"
            basis_geomopt = "GEOMOPT_BASIS"

        # '{reference molecule}[{solvent}] {func used for nmr}[{reference solvent model}]/
        # {reference basis}//{geomopt func}[{reference solvent model geomopt}]/{geomopt basis}'
        lines.append(
            f"{self._settings['h_ref']}"
            f"[{general_settings['solvent'] if not general_settings['gas-phase'] else 'gas'}] "
            f"{self._settings['func_s']}[{'SMD' if not general_settings['gas-phase'] else None}]/"
            f"{'def2-TZVP'}//{func_geomopt}[{'SMD' if not general_settings['gas-phase'] else None}]/"
            f"{basis_geomopt}\n"
        )

        lines.append(f"")

        # Write 'nmrprop.dat's
        # Determine the active elements
        active_map = {
            "H": "h_active",
            "C": "c_active",
            "F": "f_active",
            "Si": "si_active",
            "P": "p_active",
        }
        active_elements = [
            element for element in active_map.keys()
            if self._settings[active_map[element]]
        ]

        # Get the indices of the active nuclei
        active_atoms = [
            i for i in range(self.ensemble.conformers[0].geom.nat)
            if self.ensemble.conformers[0].geom.xyz[i]["element"] in active_elements
        ]

        for conf in self.ensemble.conformers:
            confdir = os.path.join(self.ensemble.workdir,
                                   self._name, conf.name)
            lines = []

            # first: atom no. | sigma(iso)
            # atom no.s according to their appearance in the xyz-file
            # NOTE: keep in mind that ANMR is written in Fortran, so the indices have to be incremented by 1
            for i, shielding in conf.results[self._name]["shieldings"]:
                lines.append(f"{i + 1:4} {shielding:.3f}\n")

            # Fill in blank lines
            for _ in range(conf.geom.nat - len(active_atoms)):
                lines.append("\n")

            # then: atom no.1 | atom no.2 | J12
            for (i, j), coupling in conf.results[self._name]["couplings"]:
                lines.append(f"{i:4} {j:4} {coupling:.3f}\n")

            logger.debug(f"Writing to {os.path.join(confdir, 'nmrprop.dat')}.")
            with open(os.path.join(confdir, "nmrprop.dat"), "w") as f:
                f.writelines(lines)

    def write_results(self, results: dict[int, any]) -> None:
        """
        Write result shieldings and/or couplings to files.
        """
        # Write results to json file
        self.write_json()
