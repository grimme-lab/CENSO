"""
Calculates the ensemble UV/Vis spectrum.
"""
from functools import reduce

from ..ensembledata import EnsembleData
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    GRIDOPTIONS,
)
from ..datastructure import MoleculeData
from ..part import CensoPart
from ..utilities import timeit
from ..logging import setup_logger

logger = setup_logger(__name__)


class UVVis(CensoPart):
    alt_name = "part6"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "nroots": {"default": 20, "range": [1, 100]},
        "threshold_bmw": {"default": 0.95, "range": [0.01, 0.99]},
        "prog": {"default": "orca", "options": PROGS},  # required
        "grid": {"default": "high+", "options": GRIDOPTIONS},  # required
        "run": {"default": False},  # required
        "template": {"default": False},  # required
        "gcp": {"default": True},  # required
    }

    _settings = {}

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        Calculation of the ensemble UV/Vis spectrum of a (previously) optimized ensemble.
        """

        # print instructions
        self.print_info()

        jobtype = ["uvvis"]

        # Select conformers based on Boltzmann weight threshold, index -1 indicates to always use the most recently
        # calculated Boltzmann weight
        self.ensemble.update_conformers(
            lambda conf: conf.bmws[-1],
            self.get_settings()["threshold_bmw"],
            boltzmann=True,
        )

        # Store the utilized Boltzmann population in order to have it in the resulting json file
        for conf in self.ensemble.conformers:
            conf.results.setdefault(self._name, {})["bmw"] = conf.bmws[-1]

        # Compile all information required for the preparation of input files in parallel execution step
        prepinfo = self.setup_prepinfo()

        # compute results
        # for structure of results from handler.execute look there
        results, failed = execute(
            self.ensemble.conformers,
            self.dir,
            self.get_settings()["prog"],
            prepinfo,
            jobtype,
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            omp=self.get_general_settings()["omp"],
            maxcores=self.get_general_settings()["maxcores"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

        # Put results into conformers
        for conf in self.ensemble.conformers:
            # store results
            conf.results.setdefault(self._name, {}).update(
                results[conf.geom.id])

        # Recalculate Boltzmann populations based on new single-point energy
        for conf in self.ensemble.conformers:
            conf.results[self._name]["gtot"] = self.gtot(conf)
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings()["temperature"],
            self._name
        )

        # TODO - Generate files to plot UV/Vis spectrum

        # Write data
        self.write_results()

    def setup_prepinfo(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        # TODO

        return prepinfo

    def gtot(self, conformer: MoleculeData) -> float:
        """
        Calculates the free enthalpy of the conformer. If any RRHO energy is found, use it in order of priority:
            optimization -> screening
        Otherwise just returns the single-point energy.
        """
        if self.get_general_settings()["evaluate_rrho"]:
            if "optimization" in conformer.results.keys():
                return conformer.results[self._name]["uvvis"]["energy"] + conformer.results["optimization"]["xtb_rrho"]["energy"]
            elif "screening" in conformer.results.keys():
                return conformer.results[self._name]["uvvis"]["energy"] + conformer.results["screening"]["xtb_rrho"]["energy"]
            else:
                return conformer.results[self._name]["uvvis"]["energy"]

    def write_results(self) -> None:
        """
        Write result shieldings and/or couplings to files.
        """
        # Write results to json file
        self.write_json()
