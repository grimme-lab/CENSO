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
    GFNOPTIONS,
)
from ..datastructure import MoleculeData
from ..part import CensoPart
from ..utilities import print, timeit, DfaHelper, format_data, SolventHelper
from ..logging import setup_logger

logger = setup_logger(__name__)


class NMR(CensoPart):
    _part_no = "4"

    _grid = "high+"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "resonance_frequency": {"default": 300.0, "range": [150.0, 1000.0]},
        "prog": {"default": "orca", "options": PROGS},  # required
        "func_j": {"default": "pbe0-d4", "options": []},
        "basis_j": {"default": "def2-TZVP", "options": []},
        "sm_j": {"default": "smd", "options": __solv_mods},
        "func_s": {"default": "pbe0-d4", "options": []},
        "basis_s": {"default": "def2-TZVP", "options": []},
        "sm_s": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
        "run": {"default": False},  # required
        "template": {"default": False},  # required
        "couplings": {"default": True},
        "shieldings": {"default": True},
        "h_active": {"default": True},
        "c_active": {"default": True},
        "f_active": {"default": False},
        "si_active": {"default": False},
        "p_active": {"default": False},
    }

    _settings = {}

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        Calculation of the ensemble NMR of a (previously) optimized ensemble.
        Note, that the ensemble will not be modified anymore.
        """

        # print instructions
        self.print_info()

        jobtype = ["nmr"]
        if not self.get_settings()["couplings"] and not self.get_settings()["shieldings"]:
            # This case should basically never happen except for user input error
            raise (
                RuntimeError(
                    "No jobtype selected. "
                    "Please select at least one of the following: couplings, shieldings"
                )
            )

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
            maxcores=self.get_general_settings()["maxcores"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

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
            success, _, failed = execute(
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

        # Recalculate Boltzmann populations based on new single-point energy
        for conf in self.ensemble.conformers:
            conf.results[self._name]["gtot"] = self.gtot(conf)

        self.ensemble.calc_boltzmannweights(
            self.get_general_settings()["temperature"],
            self._name
        )

        # Generate files for ANMR
        self.__generate_anmr()

        # Write data
        self.write_results()

    def setup_prepinfo(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        # The first condition checks if the settings are the same for shieldings and couplings calculations
        # The second and third check which one should be executed
        conds = (
            self.get_settings()["func_s"] == self.get_settings()["func_j"]
            and self.get_settings()["basis_s"] == self.get_settings()["basis_j"]
            and self.get_settings()["sm_s"] == self.get_settings()["sm_j"],
            self.get_settings()["shieldings"],
            self.get_settings()["couplings"],
        )

        # Configure the jobtypes in prepinfo according to what can be done
        # (only one calculation if all the settings are the same and/or only one type of calculation should be done,
        # otherwise both)
        # TODO - this doesn't look very nice
        if all(conds):
            prepinfo["nmr"] = {
                "func_name": DfaHelper.get_name(
                    self.get_settings()["func_s"],
                    self.get_settings()["prog"]
                ),
                "func_type": DfaHelper.get_type(self.get_settings()["func_s"]),
                "disp": DfaHelper.get_disp(self.get_settings()["func_s"]),
                "basis": self.get_settings()["basis_s"],
                "grid": "high+",  # hardcoded grid settings
                "template": self.get_settings()["template"],
                # TODO - note that GCP will be messed up if you choose one func_s/j to be a composite
                # while the other functional isn't
                "gcp": True,  # by default GCP should always be used if possible
                "sm": self.get_settings()["sm_s"],
                "h_active": self.get_settings()["h_active"],
                "c_active": self.get_settings()["c_active"],
                "f_active": self.get_settings()["f_active"],
                "si_active": self.get_settings()["si_active"],
                "p_active": self.get_settings()["p_active"],
            }
            # Only look up solvent if solvation is used
            if not self.get_general_settings()["gas-phase"]:
                prepinfo["nmr"]["solvent_key_prog"] = SolventHelper.get_solvent(
                    self.get_settings()["sm_s"], self.get_general_settings()["solvent"])
                assert prepinfo["nmr"]["solvent_key_prog"] is not None
        else:
            todo = {
                "_s": self.get_settings()["shieldings"],
                "_j": self.get_settings()["couplings"]
            }
            endings = [ending for ending in todo.keys() if todo[ending]]
            for ending in endings:
                prepinfo[f"nmr{ending}"] = {
                    "func_name": DfaHelper.get_name(
                        self.get_settings()[f"func{ending}"],
                        self.get_settings()["prog"]
                    ),
                    "func_type": DfaHelper.get_type(self.get_settings()[f"func{ending}"]),
                    "disp": DfaHelper.get_disp(self.get_settings()[f"func{ending}"]),
                    "basis": self.get_settings()[f"basis{ending}"],
                    "grid": "high+",
                    "template": self.get_settings()["template"],
                    "gcp": True,
                    "sm": self.get_settings()[f"sm{ending}"],
                    "h_active": self.get_settings()["h_active"],
                    "c_active": self.get_settings()["c_active"],
                    "f_active": self.get_settings()["f_active"],
                    "si_active": self.get_settings()["si_active"],
                    "p_active": self.get_settings()["p_active"],
                }
            # Only look up solvent if solvation is used
            if not self.get_general_settings()["gas-phase"]:
                prepinfo[f"nmr{ending}"]["solvent_key_prog"] = SolventHelper.get_solvent(
                    self.get_settings()[f"sm{ending}"], self.get_general_settings()["solvent"])
                assert prepinfo[f"nmr{ending}"]["solvent_key_prog"] is not None

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
            refinement -> optimization -> screening
        """
        if self.get_general_settings()["evaluate_rrho"]:
            if "refinement" in conformer.results.keys():
                return conformer.results[self._name]["nmr"]["energy"] + conformer.results["refinement"]["xtb_rrho"]["energy"]
            elif "optimization" in conformer.results.keys():
                return conformer.results[self._name]["nmr"]["energy"] + conformer.results["optimization"]["xtb_rrho"]["energy"]
            elif "screening" in conformer.results.keys():
                return conformer.results[self._name]["nmr"]["energy"] + conformer.results["screening"]["xtb_rrho"]["energy"]
            else:
                return conformer.results[self._name]["nmr"]["energy"] + conformer.results["nmr"]["xtb_rrho"]["energy"]

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

        # determines what to print for each conformer in each column
        printmap = {
            "ONOFF": lambda conf: "1",
            "NMR": lambda conf: f"{conf.name[4:]}",
            "CONF": lambda conf: f"{conf.name[4:]}",
            "BW": lambda conf: f"{conf.results[self._name]['bmw']:.4f}",
            "Energy": lambda conf: f"{conf.results[self._name]['nmr']['energy']:.6f}",
            "Gsolv": lambda conf: f"{0.0:.6f}",
            "mRRHO": lambda conf: f"{conf.results['optimization']['xtb_rrho']['energy']:.6f}"
            if "optimization" in conf.results.keys() and self.get_general_settings()["evaluate_rrho"]
            else f"{0.0:.6f}",
            "gi": lambda conf: f"{conf.degen}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows)

        # Write lines to file
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, 'anmr_enso')}."
        )
        with open(
            os.path.join(self.ensemble.workdir, "anmr_enso"), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Write 'anmrrc'
        # TODO - this is not finished, also don't do this for now, it's pretty straightforward to configure .anmrrc
        # manually
        """
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
                           "reference geometry as well as reference shieldings.")
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
        """
        # Write 'nmrprop.dat's and coord files
        for conf in self.ensemble.conformers:
            confdir = os.path.join(self.dir, conf.name)
            lines = []

            # first: atom no. | sigma(iso)
            # atom no.s according to their appearance in the xyz-file
            # NOTE: keep in mind that ANMR is written in Fortran, so the indices have to be incremented by 1
            for i, shielding in conf.results[self._name]["nmr"]["shieldings"]:
                lines.append(f"{i + 1:4} {shielding:.3f}\n")

            # Fill in blank lines
            for _ in range(conf.geom.nat - len(conf.results[self._name]["nmr"]["shieldings"])):
                lines.append("\n")

            # then: atom no.1 | atom no.2 | J12
            for (i, j), coupling in conf.results[self._name]["nmr"]["couplings"]:
                lines.append(f"{i + 1:4} {j + 1:4} {coupling:.3f}\n")

            logger.debug(f"Writing to {os.path.join(confdir, 'nmrprop.dat')}.")
            with open(os.path.join(confdir, "nmrprop.dat"), "w") as f:
                f.writelines(lines)

            # Write coord files
            lines = conf.geom.tocoord()
            logger.debug(f"Writing to {os.path.join(confdir, 'coord')}.")
            with open(os.path.join(confdir, "coord"), "w") as f:
                f.writelines(lines)

        print("\nGeneration of ANMR files done. Don't forget to setup your .anmrrc file.")

    def shieldings_averaging(self):
        """
        Calculate the population weighted shielding constants for the ensemble NMR spectrum.
        """
        pass

    def write_results(self) -> None:
        """
        Write result shieldings and/or couplings to files.
        """
        # Write results to json file
        self.write_json()
