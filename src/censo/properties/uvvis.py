"""
Calculates the ensemble UV/Vis spectrum.
"""

import json
import os

from ..parallel import execute
from ..params import Config
from ..utilities import SolventHelper, DfaHelper, format_data, print, Factory
from ..logging import setup_logger
from .property_calculator import PropertyCalculator
from ..part import CensoPart

logger = setup_logger(__name__)


class UVVis(PropertyCalculator):
    """
    Calculation of the ensemble UV/Vis spectrum of a (previously) optimized ensemble.
    Note, that the ensemble will not be modified anymore.
    """

    __solv_mods = {
        prog: tuple(
            t for t in Config.SOLV_MODS[prog] if t not in ("cosmors", "cosmors-fine")
        )
        for prog in Config.PROGS
    }

    _options = {
        "prog": {"default": "orca", "options": ["orca"]},  # required
        "func": {
            "default": "wb97x-d4",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in Config.PROGS},
        },
        "basis": {"default": "def2-TZVP"},
        "sm": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": Config.GFNOPTIONS},
        "nroots": {"default": 20},
        "run": {"default": False},  # required
        "template": {"default": False},  # required
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
        # Check availability of func for prog
        func = tovalidate["func"]
        if func not in cls._options["func"]["options"][tovalidate["prog"]]:
            raise ValueError(
                f"Functional {func} is not available for {tovalidate['prog']}. "
                "Check spelling w.r.t. CENSO functional naming convention (case insensitive)."
            )

        # Check sm availability for prog
        # Remember: tovalidate is always complete so we don't need .get with default None here
        sm = tovalidate["sm"]
        if sm not in cls._options["sm"]["options"][tovalidate["prog"]]:
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
        jobtype = ["uvvis"]

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

        # Ensemble averaging of excitations
        self.__excitation_averaging()

    def _setup_prepinfo(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self.name
        prepinfo["charge"] = self._ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self._ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        prepinfo["uvvis"] = {
            "func_name": DfaHelper.get_name(
                self.get_settings()["func"], self.get_settings()["prog"]
            ),
            "func_type": DfaHelper.get_type(self.get_settings()["func"]),
            "disp": DfaHelper.get_disp(self.get_settings()["func"]),
            "basis": self.get_settings()["basis"],
            "grid": "high+",  # hardcoded grid settings
            "template": self.get_settings()["template"],
            "gcp": False,  # GCP is not necessary for spectra calculations
            "nroots": self.get_settings()["nroots"],
        }
        # Only look up solvent if solvation is used
        if not self.get_general_settings()["gas-phase"]:
            prepinfo["uvvis"]["sm"] = self.get_settings()["sm"]
            prepinfo["uvvis"]["solvent_key_prog"] = SolventHelper.get_solvent(
                self.get_settings()["sm"], self.get_general_settings()["solvent"]
            )

        return prepinfo

    def __excitation_averaging(self):
        """
        Calculates population weighted excitation parameters.
        """
        # Calculate epsilon_max (maximum extinctions) for each excitation, weighted by population
        # eps is a list of tuples that contain each excitation wavelength with the respective epsilon_max
        eps = []
        for conf in self._ensemble.conformers:
            for excitation in self.data["results"][conf.name]["excitations"]:
                epsilon_max = (
                    self.data["results"][conf.name]["bmw"] * excitation["osc_str"]
                )
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
            f"Writing to {os.path.join(os.getcwd(), f'{self._part_nos[self.name]}_{self.name.upper()}.out')}."
        )
        with open(
            os.path.join(
                os.getcwd(),
                f"{self._part_nos[self.name]}_{self.name.upper()}.out",
            ),
            "w",
            newline=None,
        ) as outfile:
            outfile.writelines(lines)

        # Dump data into json
        with open(os.path.join(os.getcwd(), "excitations.json"), "w") as f:
            json.dump(eps, f, indent=4)

    def _write_results(self) -> None:
        """
        Write result excitations to files.
        """
        # Write results to json file
        self._write_json()


Factory.register_builder("uvvis", UVVis)
