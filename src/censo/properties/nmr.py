"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""
from functools import reduce

from ..core import CensoCore
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
)
from ..part import CensoPart
from ..utilities import print, timeit, DfaHelper, setup_logger

logger = setup_logger(__name__)


class NMR(CensoPart):
    alt_name = "part2"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "resonance_frequency": {"default": 300.0, "range": [150.0, 1000.0]},
        "threshold_bmw": {"default": 0.95, "range": [0.01, 0.99]},
        "prog": {"default": "orca", "options": PROGS},
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
        "grid": {"default": "high+", "options": GRIDOPTIONS},
        "run": {"default": False},
        "couplings": {"default": True},
        "shieldings": {"default": True},
        "h_active": {"default": True},
        "c_active": {"default": True},
        "f_active": {"default": False},
        "si_active": {"default": False},
        "p_active": {"default": False},
    }

    _settings = {}

    def __init__(self, censo: CensoCore):
        super().__init__(censo)

        # Set the correct name for 'func_s' and 'func_j'
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
        self.core.update_conformers(
            lambda conf: conf.bmws[-1],
            self._instructions["threshold_bmw"],
            boltzmann=True,
        )

        # Store the utilized Boltzmann population in order to have it in the resulting json file
        for conf in self.core.conformers:
            conf.results[self._name]["bmw"] = conf.bmws[-1]

        # Execute jobs in parallel
        results = execute(self.core.conformers, self._instructions, self.dir)

        # Put results into conformers
        for conf in self.core.conformers:
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
        pass

    def write_results(self, results: dict[int, any]) -> None:
        """
        Write result shieldings and/or couplings to files.
        """
        # Write results to json file
        self.write_json()
