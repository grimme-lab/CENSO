"""
Calculates the ensemble UV/Vis spectrum.
"""
from functools import reduce
import json
import os

from ..ensembledata import EnsembleData
from ..parallel import execute
from ..params import (SOLV_MODS, PROGS, GFNOPTIONS)
from ..datastructure import MoleculeData
from ..part import CensoPart
from ..utilities import timeit, SolventHelper, DfaHelper, format_data
from ..logging import setup_logger

logger = setup_logger(__name__)


class UVVis(CensoPart):
    _part_no = "6"

    __solv_mods = tuple(t for t in reduce(lambda x, y: x + y, SOLV_MODS.values()) if t not in ("cosmors", "cosmors-fine"))

    _options = {
        "prog": {
            "default": "orca",
            "options": PROGS
        },  # required
        "func": {
            "default": "wb97x-d4"
        },
        "basis": {
            "default": "def2-TZVP"
        },
        "sm": {
            "default": "smd",
            "options": __solv_mods
        },
        "gfnv": {
            "default": "gfn2",
            "options": GFNOPTIONS
        },
        "nroots": {
            "default": 20
        },
        "run": {
            "default": False
        },  # required
        "template": {
            "default": False
        },  # required
        "gcp": {
            "default": True
        },  # required
    }

    _settings = {}

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self, ncores: int) -> None:
        """
        Calculation of the ensemble UV/Vis spectrum of a (previously) optimized ensemble.
        Note, that the ensemble will not be modified anymore.
        """

        # print instructions
        self.print_info()

        jobtype = ["uvvis"]

        # Compile all information required for the preparation of input files in parallel execution step
        prepinfo = self.setup_prepinfo()

        # compute results
        # for structure of results from handler.execute look there
        success, _, failed = execute(
            self.ensemble.conformers,
            self.dir,
            self.get_settings()["prog"],
            prepinfo,
            jobtype,
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            omp=self.get_general_settings()["omp"],
            maxcores=ncores,
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

        # Set energy values to be used later
        self.__set_energy(ncores)

        # Recalculate Boltzmann populations based on new single-point energy
        for conf in self.ensemble.conformers:
            conf.results[self._name]["gtot"] = self.gtot(conf)

        self.ensemble.calc_boltzmannweights(
            self.get_general_settings()["temperature"], self._name)

        # Ensemble averaging of excitations
        self.__excitation_averaging()

        # Write data
        self.write_results()

    def setup_prepinfo(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        prepinfo["uvvis"] = {
            "func_name":
            DfaHelper.get_name(self.get_settings()["func"],
                               self.get_settings()["prog"]),
            "func_type":
            DfaHelper.get_type(self.get_settings()["func"]),
            "disp":
            DfaHelper.get_disp(self.get_settings()["func"]),
            "basis":
            self.get_settings()["basis"],
            "grid":
            "high+",  # hardcoded grid settings
            "template":
            self.get_settings()["template"],
            # while the other functional isn't
            "gcp":
            True,  # by default GCP should always be used if possible
            "nroots":
            self.get_settings()["nroots"],
        }
        # Only look up solvent if solvation is used
        if not self.get_general_settings()["gas-phase"]:
            prepinfo["uvvis"]["sm"] = self.get_settings()["sm"]
            prepinfo["uvvis"]["solvent_key_prog"] = SolventHelper.get_solvent(
                self.get_settings()["sm"],
                self.get_general_settings()["solvent"])
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
            prepinfo["xtb_rrho"][
                "solvent_key_xtb"] = SolventHelper.get_solvent(
                    self.get_general_settings()["sm_rrho"],
                    self.get_general_settings()["solvent"])
            assert prepinfo["xtb_rrho"]["solvent_key_xtb"] is not None

        return prepinfo

    def __set_energy(self, ncores: int):
        """
        Looks through results to set energy values.
        Order of preference:
            refinement -> optimization -> screening -> prescreening

        If None of these are found, the energies of the UVVis calculations will be used.
        """
        using_part = None
        for partname in [
                "prescreening", "screening", "optimization", "refinement"
        ]:
            # This way, the most high-level partname should get stuck in using_part
            if all(partname in conf.results.keys()
                   for conf in self.ensemble.conformers):
                using_part = partname

        # If using_part stays None the UVVis energies must be used and xtb_rrho might be necessary
        # TODO - this is not nice, DRY
        if using_part is None:
            if self.get_general_settings()["evaluate_rrho"]:
                jobtype = ["xtb_rrho"]
                prepinfo = self.setup_prepinfo_rrho()

                # Run RRHO calculation
                success, _, failed = execute(
                    self.ensemble.conformers,
                    self.dir,
                    self.get_settings()["prog"],
                    prepinfo,
                    jobtype,
                    copy_mo=self.get_general_settings()["copy_mo"],
                    balance=self.get_general_settings()["balance"],
                    omp=self.get_general_settings()["omp"],
                    maxcores=self.ensemble.runinfo["maxcores"],
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Remove failed conformers
                self.ensemble.remove_conformers(failed)
            for conf in self.ensemble.conformers:
                conf.results[self._name]["energy"] = conf.results[
                    self._name]["uvvis"]["energy"]
                conf.results[self._name]["gsolv"] = 0.0
                conf.results[self._name]["grrho"] = conf.results[
                    self._name].get("xtb_rrho", {"energy": 0.0})["energy"]
        elif using_part == "optimization":
            for conf in self.ensemble.conformers:
                conf.results[self._name]["energy"] = conf.results[using_part][
                    "xtb_opt"]["energy"]
                conf.results[self._name]["gsolv"] = 0.0
                conf.results[
                    self._name]["grrho"] = conf.results[using_part].get(
                        "xtb_rrho", {"energy": 0.0})["energy"]
        elif using_part in ["screening", "refinement"]:
            if all("gsolv" in conf.results[using_part].keys()
                   for conf in self.ensemble.conformers):
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["energy"] = conf.results[
                        using_part]["gsolv"]["energy_gas"]
                    conf.results[self._name]["gsolv"] = conf.results[
                        using_part]["gsolv"]["gsolv"]
                    conf.results[
                        self._name]["grrho"] = conf.results[using_part].get(
                            "xtb_rrho", {"energy": 0.0})["energy"]
            else:
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["energy"] = conf.results[
                        using_part]["sp"]["energy"]
                    conf.results[self._name]["gsolv"] = 0.0
                    conf.results[
                        self._name]["grrho"] = conf.results[using_part].get(
                            "xtb_rrho", {"energy": 0.0})["energy"]
        elif using_part == "prescreening":
            if all("xtb_gsolv" in conf.results[using_part].keys()
                   for conf in self.ensemble.conformers):
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["energy"] = conf.results[
                        using_part]["sp"]["energy"]
                    conf.results[self._name]["gsolv"] = conf.results[
                        using_part]["xtb_gsolv"]["gsolv"]
                    conf.results[self._name]["grrho"] = 0.0

    def gtot(self, conformer: MoleculeData) -> float:
        resdict = conformer.results[self._name]
        return resdict["energy"] + resdict["gsolv"] + resdict["grrho"]

    def __excitation_averaging(self):
        """
        Calculates population weighted excitation parameters.
        """
        # Calculate epsilon_max (maximum extinctions) for each excitation, weighted by population
        # eps is a list of tuples that contain each excitation wavelength with the respective epsilon_max
        eps = []
        for conf in self.ensemble.conformers:
            for excitation in conf.results[self._name]["uvvis"]["excitations"]:
                epsilon_max = conf.bmws[-1] * excitation["osc_str"]
                eps.append((excitation["wavelength"], epsilon_max, conf.name))

        # Print table
        headers = ["λ", "ε_max", "Origin. CONF#"]

        units = ["[nm]", "", ""]

        printmap = {
            "λ": lambda exc: f"{exc[0]:.2f}",
            "ε_max": lambda exc: f"{exc[1]:.6f}",
            "Origin. CONF#": lambda exc: f"{exc[2]}",
        }

        rows = [[printmap[header](exc) for header in headers] for exc in eps]

        lines = format_data(headers, rows, units=units)

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # write lines to file
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, f'{self._part_no}_{
                                       self._name.upper()}.out')}."
        )
        with open(os.path.join(self.ensemble.workdir,
                               f"{self._part_no}_{self._name.upper()}.out"),
                  "w",
                  newline=None) as outfile:
            outfile.writelines(lines)

        # Dump data into json
        with open(os.path.join(self.ensemble.workdir, "excitations.json"),
                  "w") as f:
            json.dump(eps, f, indent=4)

    def write_results(self) -> None:
        """
        Write result excitations to files.
        """
        # Write results to json file
        self.write_json()
