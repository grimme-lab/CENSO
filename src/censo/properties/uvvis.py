"""
Calculates the ensemble UV/Vis spectrum.
"""
from functools import reduce

from ..ensembledata import EnsembleData
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    GFNOPTIONS
)
from ..datastructure import MoleculeData
from ..part import CensoPart
from ..utilities import timeit, SolventHelper, DfaHelper
from ..logging import setup_logger

logger = setup_logger(__name__)


class UVVis(CensoPart):
    alt_name = "part6"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "prog": {"default": "orca", "options": PROGS},  # required
        "func": {"default": "wb97x-d4", "options": []},
        "basis": {"default": "def2-TZVP", "options": []},
        "sm": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
        "nroots": {"default": 20, "range": [1, 100]},
        "threshold_bmw": {"default": 0.95, "range": [0.01, 0.99]},
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

        # Preselect conformers based on Boltzmann weight threshold, index -1 indicates to always use the most recently
        # calculated Boltzmann weight
        preselection = False
        if all(len(conf.bmws) > 0 for conf in self.ensemble.conformers):
            self.ensemble.update_conformers(
                lambda conf: conf.bmws[-1],
                self.get_settings()["threshold_bmw"],
                boltzmann=True,
            )
            preselection = True

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

        # If RRHO contribution should be included and there was no previous ensemble optimization, calculate RRHO
        if not (
                any(
                    part in conf.results.keys() for conf in self.ensemble.conformers
                    for part in ["screening", "optimization", "refinement"]
                ) and self.get_general_settings()["evaluate_rrho"]
        ):
            jobtype = ["xtb_rrho"]
            prepinfo = self.setup_prepinfo_rrho()

            # Run RRHO calculation
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

            for conf in self.ensemble.conformers:
                # update results for each conformer
                conf.results[self._name].update(results[id(conf)])

        # Recalculate Boltzmann populations based on new single-point energy
        for conf in self.ensemble.conformers:
            conf.results[self._name]["gtot"] = self.gtot(conf)

        self.ensemble.calc_boltzmannweights(
            self.get_general_settings()["temperature"],
            self._name
        )

        # In case there was no ensemble optimization done before for preselection, cut down ensemble here
        if not preselection:
            self.ensemble.update_conformers(
                lambda conf: conf.bmws[-1],
                self.get_settings()["threshold_bmw"],
                boltzmann=True,
            )

            # Recalculate Boltzmann populations to be used by ANMR
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

        prepinfo["uvvis"] = {
            "func_name": DfaHelper.get_name(
                self.get_settings()["func"],
                self.get_settings()["prog"]
            ),
            "func_type": DfaHelper.get_type(self.get_settings()["func"]),
            "disp": DfaHelper.get_disp(self.get_settings()["func"]),
            "basis": self.get_settings()["basis"],
            "grid": "high+",  # hardcoded grid settings
            "template": self.get_settings()["template"],
            # while the other functional isn't
            "gcp": True,  # by default GCP should always be used if possible
            "sm": self.get_settings()["sm"],
        }
        # Only look up solvent if solvation is used
        if not self.get_general_settings()["gas-phase"]:
            prepinfo["uvvis"]["solvent_key_prog"] = SolventHelper.get_solvent(
                self.get_settings()["sm"], self.get_general_settings()["solvent"])
            assert prepinfo["uvvis"]["solvent_key_prog"] is not None

        return prepinfo

    def setup_prepinfo_rrho(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        prepinfo["xtb_rrho"] = {
            "gfnv": self.get_settings()["gfnv"],
        }
        # Only lookup solvent if solvation should be used
        if not self.get_general_settings()["gas-phase"]:
            prepinfo["xtb_rrho"]["solvent_key_prog"] = SolventHelper.get_solvent(
                self.get_general_settings()["sm_rrho"], self.get_general_settings()["solvent"])
            assert prepinfo["xtb_rrho"]["solvent_key_prog"] is not None

        return prepinfo

    def gtot(self, conformer: MoleculeData) -> float:
        """
        Calculates the free enthalpy of the conformer. If any previous RRHO energy is found, use it in order of priority:
            optimization -> screening
        """
        if self.get_general_settings()["evaluate_rrho"]:
            if "optimization" in conformer.results.keys():
                return conformer.results[self._name]["uvvis"]["energy"] + conformer.results["optimization"]["xtb_rrho"]["energy"]
            elif "screening" in conformer.results.keys():
                return conformer.results[self._name]["uvvis"]["energy"] + conformer.results["screening"]["xtb_rrho"]["energy"]
            else:
                return conformer.results[self._name]["uvvis"]["energy"] + conformer.results["uvvis"]["xtb_rrho"]["energy"]

    def write_results(self) -> None:
        """
        Write result excitations to files.
        """
        # Write results to json file
        self.write_json()
