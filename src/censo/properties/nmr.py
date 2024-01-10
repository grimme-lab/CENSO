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
from ..datastructure import MoleculeData
from ..part import CensoPart
from ..utilities import print, timeit, DfaHelper, setup_logger, format_data

logger = setup_logger(__name__)


class NMR(CensoPart):
    alt_name = "part4"

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

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        Calculation of the ensemble NMR of a (previously) optimized ensemble.
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
                "grid": self.get_settings()["grid"],
                "template": self.get_settings()["template"],
                # TODO - note that GCP will be messed up if you choose one func_s/j to be a composite
                # while the other functional isn't
                "gcp": self.get_settings()["gcp"],
                "sm": self.get_settings()["sm_s"],
                "solvent_key_prog": SOLVENTS_DB.get(
                    self.get_general_settings()["solvent"]
                )[self.get_settings()["sm_s"]][1],
                "h_active": self.get_settings()["h_active"],
                "c_active": self.get_settings()["c_active"],
                "f_active": self.get_settings()["f_active"],
                "si_active": self.get_settings()["si_active"],
                "p_active": self.get_settings()["p_active"],
            }
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
                    "grid": self.get_settings()["grid"],
                    "template": self.get_settings()["template"],
                    "gcp": self.get_settings()["gcp"],
                    "sm": self.get_settings()[f"sm{ending}"],
                    "solvent_key_prog": SOLVENTS_DB.get(
                        self.get_general_settings()["solvent"]
                    )[self.get_settings()[f"sm{ending}"]][1],
                    "h_active": self.get_settings()["h_active"],
                    "c_active": self.get_settings()["c_active"],
                    "f_active": self.get_settings()["f_active"],
                    "si_active": self.get_settings()["si_active"],
                    "p_active": self.get_settings()["p_active"],
                }

        return prepinfo

    def gtot(self, conformer: MoleculeData) -> float:
        """
        Calculates the free enthalpy of the conformer. If any RRHO energy is found, use it in order of priority:
            optimization -> screening
        Otherwise just returns the single-point energy.
        """
        if self.get_general_settings()["evaluate_rrho"]:
            if "optimization" in conformer.results.keys():
                return conformer.results[self._name]["nmr"]["energy"] + conformer.results["optimization"]["xtb_rrho"]["energy"]
            elif "screening" in conformer.results.keys():
                return conformer.results[self._name]["nmr"]["energy"] + conformer.results["screening"]["xtb_rrho"]["energy"]
            else:
                return conformer.results[self._name]["nmr"]["energy"]

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
            "ONOFF": lambda conf: 1,
            "NMR": lambda conf: conf.name[4:],
            "CONF": lambda conf: conf.name[4:],
            "BW": lambda conf: f"{conf.results[self._name]['bmw']:.4f}",
            "Energy": lambda conf: f"{conf.results[self._name]['nmr']['energy']:.6f}",
            "Gsolv": lambda conf: f"{0.0:.6f}",
            "mRRHO": lambda conf: f"{conf.results['optimization']['xtb_rrho']['energy']:.6f}"
            if "optimization" in conf.results.keys() and self.get_general_settings()["evaluate_rrho"]
            else f"{0.0:.6f}",
            "gi": lambda conf: conf.degen,
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
        # Write 'nmrprop.dat's
        for conf in self.ensemble.conformers:
            confdir = os.path.join(self.ensemble.workdir,
                                   self._name, conf.name)
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
                lines.append(f"{i:4} {j:4} {coupling:.3f}\n")

            logger.debug(f"Writing to {os.path.join(confdir, 'nmrprop.dat')}.")
            with open(os.path.join(confdir, "nmrprop.dat"), "w") as f:
                f.writelines(lines)

        print("\nGeneration of ANMR files done. Don't forget to setup your .anmrrc file.")

    def write_results(self) -> None:
        """
        Write result shieldings and/or couplings to files.
        """
        # Write results to json file
        self.write_json()
