"""
Calculates the ensemble NMR spectrum for all active nuclei.
"""

import os

from ..parallel import execute
from ..params import Config
from .property_calculator import PropertyCalculator
from ..utilities import print, DfaHelper, format_data, SolventHelper, Factory
from ..logging import setup_logger
from ..part import CensoPart

logger = setup_logger(__name__)


class NMR(PropertyCalculator):
    _grid = "nmr"

    __solv_mods = {
        prog: tuple(
            t for t in Config.SOLV_MODS[prog] if t not in ("cosmors", "cosmors-fine")
        )
        for prog in Config.PROGS
    }

    _options = {
        "resonance_frequency": {"default": 300.0},
        "ss_cutoff": {"default": 8.0},
        "prog": {"default": "orca", "options": Config.PROGS},  # required
        "func_j": {
            "default": "pbe0-d4",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in Config.PROGS},
        },
        "basis_j": {"default": "def2-TZVP"},
        "sm_j": {"default": "smd", "options": __solv_mods},
        "func_s": {
            "default": "pbe0-d4",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in Config.PROGS},
        },
        "basis_s": {"default": "def2-TZVP"},
        "sm_s": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": Config.GFNOPTIONS},
        "run": {"default": False},  # required
        "template": {"default": False},  # required
        "couplings": {"default": True},
        "fc_only": {"default": True},
        "shieldings": {"default": True},
        "h_active": {"default": True},
        "c_active": {"default": True},
        "f_active": {"default": False},
        "si_active": {"default": False},
        "p_active": {"default": False},
    }

    _settings = {}

    @classmethod
    def _validate(cls, tovalidate: dict[str, any]) -> None:
        """
        Validates the type of each setting in the given dict. Also potentially validate if the setting is allowed by
        checking with cls._options.
        This is the part-specific version of the method. It will run the general validation first and then
        check part-specific logic.

        Args:
            tovalidate (dict[str, any]): The dict containing the settings to be validated.

        Returns:
            None

        Raises:
            ValueError: If the setting is not allowed or the value is not within the allowed options.
        """
        # General validation
        super()._validate(tovalidate)

        # Part-specific validation
        # NOTE: tovalidate is always complete
        for ending in ["_s", "_j"]:
            mapping = {"_s": "shieldings", "_j": "couplings"}

            # Only check settings for shieldings/couplings calculations if they are actually turned on
            if tovalidate[mapping[ending]]:
                # Check availability of func for prog
                func = tovalidate[f"func{ending}"]
                if (
                    func
                    not in cls._options[f"func{ending}"]["options"][tovalidate["prog"]]
                ):
                    raise ValueError(
                        f"Functional {func} is not available for {tovalidate['prog']}. "
                        "Check spelling w.r.t. CENSO functional naming convention (case insensitive)."
                    )

                # Check sm availability for prog
                # Remember: tovalidate is always complete so we don't need .get with default None here
                sm = tovalidate[f"sm{ending}"]
                if sm not in cls._options[f"sm{ending}"]["options"][tovalidate["prog"]]:
                    raise ValueError(
                        f"Solvent model {sm} not available for {tovalidate['prog']}."
                    )

                # Check solvent availability for sm
                if (
                    cls.get_general_settings()["solvent"]
                    not in CensoPart._options["solvent"]["options"][sm]
                ):
                    raise ValueError(
                        f"Solvent {cls.get_general_settings()['solvent']} is not available for {sm}. "
                        "Please create an issue on GitHub if you think this is incorrect."
                    )

                # dummy/template functionality not implemented yet for TM
                if tovalidate["prog"] == "tm" and (
                    func == "dummy" or tovalidate.get("template", False)
                ):
                    raise NotImplementedError(
                        "Dummy and template functionality is not implemented yet for use with TURBOMOLE."
                    )

    def _property(self) -> None:
        """
        Calculation of the ensemble NMR of a (previously) optimized ensemble.
        Note, that the ensemble will not be modified anymore.
        """
        jobtype = ["nmr"]
        if (
            not self.get_settings()["couplings"]
            and not self.get_settings()["shieldings"]
        ):
            # This case should basically never happen except for user input error
            raise (
                RuntimeError(
                    "No jobtype selected. "
                    "Please select at least one of the following: couplings, shieldings"
                )
            )

        # Compile all information required for the preparation of input files in parallel execution step
        prepinfo = self._setup_prepinfo()

        # compute results
        # for structure of results from handler.execute look there
        results, failed = execute(
            self._ensemble.conformers,
            self._dir,
            self.get_settings()["prog"],
            prepinfo,
            jobtype,
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self._ensemble.remove_conformers(failed)

        # Update results
        self._update_results(results)

        # Generate files for ANMR
        self.__generate_anmr()

    def _setup_prepinfo(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self.name
        prepinfo["charge"] = self._ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self._ensemble.runinfo.get("unpaired")
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
                    self.get_settings()["func_s"], self.get_settings()["prog"]
                ),
                "func_type": DfaHelper.get_type(self.get_settings()["func_s"]),
                "disp": DfaHelper.get_disp(self.get_settings()["func_s"]),
                "basis": self.get_settings()["basis_s"],
                "grid": self._grid,  # hardcoded grid settings
                "template": self.get_settings()["template"],
                "gcp": False,  # GCP is not necessary for spectra calculations
                "fc_only": self.get_settings()["fc_only"],
                "ss_cutoff": self.get_settings()["ss_cutoff"],
                "h_active": self.get_settings()["h_active"],
                "c_active": self.get_settings()["c_active"],
                "f_active": self.get_settings()["f_active"],
                "si_active": self.get_settings()["si_active"],
                "p_active": self.get_settings()["p_active"],
            }
            # Only look up solvent if solvation is used
            if not self.get_general_settings()["gas-phase"]:
                prepinfo["nmr"]["sm"] = self.get_settings()["sm_s"]
                prepinfo["nmr"]["solvent_key_prog"] = SolventHelper.get_solvent(
                    self.get_settings()["sm_s"], self.get_general_settings()["solvent"]
                )
        else:
            todo = {
                "_s": self.get_settings()["shieldings"],
                "_j": self.get_settings()["couplings"],
            }
            endings = [ending for ending, active in todo.items() if active]
            for ending in endings:
                prepinfo[f"nmr{ending}"] = {
                    "func_name": DfaHelper.get_name(
                        self.get_settings()[f"func{ending}"],
                        self.get_settings()["prog"],
                    ),
                    "func_type": DfaHelper.get_type(
                        self.get_settings()[f"func{ending}"]
                    ),
                    "disp": DfaHelper.get_disp(self.get_settings()[f"func{ending}"]),
                    "basis": self.get_settings()[f"basis{ending}"],
                    "grid": "high+",
                    "template": self.get_settings()["template"],
                    "gcp": False,  # GCP is not necessary for spectra calculations
                    "fc_only": self.get_settings()["fc_only"],
                    "ss_cutoff": self.get_settings()["ss_cutoff"],
                    "h_active": self.get_settings()["h_active"],
                    "c_active": self.get_settings()["c_active"],
                    "f_active": self.get_settings()["f_active"],
                    "si_active": self.get_settings()["si_active"],
                    "p_active": self.get_settings()["p_active"],
                }
            # Only look up solvent if solvation is used
            if not self.get_general_settings()["gas-phase"]:
                prepinfo[f"nmr{ending}"]["sm"] = self.get_settings()[f"sm{ending}"]
                prepinfo[f"nmr{ending}"]["solvent_key_prog"] = (
                    SolventHelper.get_solvent(
                        self.get_settings()[f"sm{ending}"],
                        self.get_general_settings()["solvent"],
                    )
                )

        return prepinfo

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
            "BW": lambda conf: f"{self.data['results'][conf.name]['bmw']:.4f}",
            "Energy": lambda conf: f"{self.data['results'][conf.name]['energy']:.6f}",
            "Gsolv": lambda conf: f"{self.data['results'][conf.name]['gsolv']:.6f}",
            "mRRHO": lambda conf: f"{self.data['results'][conf.name]['grrho']:.6f}",
            "gi": lambda conf: f"{conf.degen}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self._ensemble.conformers
        ]

        lines = format_data(headers, rows)

        # Write lines to file
        logger.debug(f"Writing to {os.path.join(os.getcwd(), 'anmr_enso')}.")
        with open(os.path.join(os.getcwd(), "anmr_enso"), "w", newline=None) as outfile:
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
        if all("optimization" in conf.results.keys() for conf in self._ensemble.conformers):
            from .._ensembleopt.optimization import Optimization
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
        for conf in self._ensemble.conformers:
            confdir = os.path.join(self._dir, conf.name)
            lines = []

            # first: atom no. | sigma(iso)
            # atom no.s according to their appearance in the xyz-file
            # NOTE: keep in mind that ANMR is written in Fortran, so the indices have to be incremented by 1
            for i, shielding in self.data["results"][conf.name]["nmr"]["shieldings"]:
                lines.append(f"{i + 1:4} {shielding:.3f}\n")

            # Fill in blank lines
            for _ in range(
                conf.geom.nat
                - len(self.data["results"][conf.name]["nmr"]["shieldings"])
            ):
                lines.append("\n")

            # then: atom no.1 | atom no.2 | J12
            for (i, j), coupling in self.data["results"][conf.name]["nmr"]["couplings"]:
                lines.append(f"{i + 1:4} {j + 1:4} {coupling:.3f}\n")

            logger.debug(f"Writing to {os.path.join(confdir, 'nmrprop.dat')}.")
            with open(os.path.join(confdir, "nmrprop.dat"), "w") as f:
                f.writelines(lines)

        print(
            "\nGeneration of ANMR files done. Don't forget to setup your .anmrrc file."
        )

    def _shieldings_averaging(self):
        """
        Calculate the population weighted shielding constants for the ensemble NMR spectrum.
        """
        # For each conf in conformers
        # add shift shielding value multiplied by bmw to a list
        pass

    def _write_results(self) -> None:
        """
        Write result shieldings and/or couplings to files.
        """
        # Write results to json file
        self._write_json()


Factory.register_builder("nmr", NMR)
