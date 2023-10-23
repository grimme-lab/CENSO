from typing import Dict, Any, Type, Union
import os
import ast
from functools import reduce
import json

from censo.utilities import timeit
from censo.core import CensoCore
from censo.settings import CensoRCParser
from censo.cfg import (
    PLENGTH,
    DIGILEN,
    OMPMIN,
    OMPMAX,
    ASSETS_PATH,
)

"""
Part class as parent class for all parts of the calculation to
implement complete OOP approach.
"""


class CensoPart:
    __cosmors_param = {
        "12-normal": "BP_TZVP_C30_1201.ctd",
        "13-normal": "BP_TZVP_C30_1301.ctd",
        "14-normal": "BP_TZVP_C30_1401.ctd",
        "15-normal": "BP_TZVP_C30_1501.ctd",
        "16-normal": "BP_TZVP_C30_1601.ctd",
        "17-normal": "BP_TZVP_C30_1701.ctd",
        "18-normal": "BP_TZVP_18.ctd",
        "19-normal": "BP_TZVP_19.ctd",
        "12-fine": "BP_TZVPD_FINE_HB2012_C30_1201.ctd",
        "13-fine": "BP_TZVPD_FINE_HB2012_C30_1301.ctd",
        "14-fine": "BP_TZVPD_FINE_C30_1401.ctd",
        "15-fine": "BP_TZVPD_FINE_C30_1501.ctd",
        "16-fine": "BP_TZVPD_FINE_C30_1601.ctd",
        "17-fine": "BP_TZVPD_FINE_C30_1701.ctd",
        "18-fine": "BP_TZVPD_FINE_18.ctd",
        "19-fine": "BP_TZVPD_FINE_19.ctd",
    }

    with open(os.path.join(ASSETS_PATH, "censo_solvents_db.json"), "r") as solv_file:
        __solvents_db = json.load(solv_file)

    _options = {
        "general": {
            "procs": {
                "default": 1,
                "range": [
                    1,
                    128
                ]
            },
            "omp": {
                "default": 4,
                "range": [
                    OMPMIN,
                    OMPMAX
                ]
            },
            "imagthr": {
                "default": -100.0,
                "range": [
                    -300.0,
                    0.0
                ]
            },
            "sthr": {
                "default": 0.0,
                "range": [
                    0.0,
                    100.0
                ]
            },
            "scale": {
                "default": 1.0,
                "range": [
                    0.0,
                    1.0
                ]
            },
            "temperature": {
                "default": 298.15,
                "range": [
                    1e-05,
                    2000.0
                ]
            },
            "solvent": {
                "default": "h2o",
                "options": [k for k in __solvents_db.keys()]
            },
            "sm_rrho": {
                "default": "alpb",
                "options": [
                    "alpb",
                    "gbsa"
                ]
            },
            "cosmorsparam": {
                "default": "12-normal",
                "options": [k for k in __cosmors_param.keys()]
            },
            "multitemp": {
                "default": True
            },
            "evaluate_rrho": {
                "default": True
            },
            "consider_sym": {
                "default": True
            },
            "bhess": {
                "default": True
            },
            "rmsdbias": {
                "default": False
            },
            "progress": {
                "default": False
            },
            "check": {
                "default": False
            },
            "balance": {
                "default": True
            },
            "vapor_pressure": {
                "default": False
            },
            "nmrmode": {
                "default": False
            },
            "gas-phase": {
                "default": False
            },
            "copy_mo": {
                "default": True
            },
            "trange": {
                "default": [
                    273.15,
                    378.15,
                    5
                ]
            }
        },
    }

    _settings = {}

    @classmethod
    def get_settings(cls):
        return cls._settings

    @classmethod
    def set_settings(cls, settings: Dict[str, Any]):
        cls._validate(settings)
        settings = cls._complete(settings)
        cls._settings = settings

    @classmethod
    def _complete(cls, tocomplete: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """
        fill in missing settings with default values
        """
        for part, settings in cls._options.items():
            for setting, definition in settings.items():
                if setting not in tocomplete[part]:
                    tocomplete[part][setting] = f"{definition['default']}"

        return tocomplete

    @classmethod
    def _validate(cls, tovalidate: Dict[str, Any]) -> None:
        """
        validate the type of each setting in the given dict
        also potentially validate if the setting is allowed by checking with cls._options

        raises exceptions if some setting's type is invalid or the settings value is not within predefined options
        """
        # go through each section and try to validate each setting's type
        for part in tovalidate.keys():
            for setting_name in tovalidate[part]:
                # try to get the settings's type from the default value in the _options dict
                try:
                    setting_type = type(cls._options[part][setting_name]["default"])
                except KeyError:
                    # KeyError means that the setting does not exist, therefore the setting is removed
                    tovalidate[part].pop(setting_name)
                    continue

                # try to cast the setting-string into the correct type
                try:
                    if setting_type == bool:
                        setting_value = {"True": True, "False": False}.get(tovalidate[part][setting_name])
                    elif setting_type == list:
                        setting_value = ast.literal_eval(tovalidate[part][setting_name])
                    else:
                        setting_value = setting_type(tovalidate[part][setting_name])
                # if that's not possible raise an exception
                # NOTE: KeyError is raised when the conversion for bools fails
                except (ValueError, KeyError):
                    raise ValueError(
                        f"Value '{tovalidate[part][setting_name]}' is not allowed for setting '{setting_name}' in part '{part}'")

                # now check if the setting is allowed
                # for strings check if string is within a list of allowed values
                if setting_type == str:
                    options = cls._options[part][setting_name]['options']
                    if setting_value not in options and len(options) > 0:
                        # Only check if there are options
                        # This is fatal so an exception is raised
                        raise ValueError(
                            f"Value '{setting_value}' is not allowed for setting '{setting_name}' in part '{part}'.")
                # for numeric values check if value is within a range
                elif setting_type in (int, float):
                    interval = cls._options[part][setting_name]['range']
                    if not interval[0] <= setting_value <= interval[1]:
                        # This is fatal so an exception is raised
                        raise ValueError(
                            f"Value '{setting_value}' is out of range ({interval[0]},{interval[1]}) for setting '{setting_name}' in part '{part}'.")
                # NOTE: there is no check for complex types yet (i.e. lists)

                # set the value in the dict tovalidate to the casted value
                tovalidate[part][setting_name] = setting_value

    def __init__(self, core: CensoCore, settings: CensoRCParser, part: str):
        self.core: CensoCore = core
        self.settings: CensoRCParser = settings

        # contains settings grabbed from CensoSettings instance, such as general settings etc.
        self._instructions: Dict[str, Any]

        # grabs the settings required for this part from the passed 'CensoSettings' instance
        paths = settings.settings_current.get("paths")
        general = settings.settings_current.get("general")
        specific = settings.settings_current.get(part)
        self._instructions = {**paths, **general, **specific}

        # add some additional settings to _instructions so that the processors don't have to do any lookups
        # NOTE: [1] auto-selects replacement solvent (TODO - print warning!)
        self._instructions["solvent_key_xtb"] = settings.solvents_db.get(self._instructions["solvent"])["xtb"][1]
        if 'sm' in self._instructions.keys():
            self._instructions[f"solvent_key_prog"] = \
                settings.solvents_db.get(self._instructions["solvent"])[self._instructions["sm"]][1]
            # TODO - doesn't work yet for parts where 'func' keyword doesn't exist or there are multiple functionals
        self._instructions["func_type"] = settings.dfa_settings.get_type(self._instructions["func"])

        # add 'charge' and 'unpaired' to instructions
        self._instructions["charge"] = core.runinfo.get("charge")
        self._instructions["unpaired"] = core.runinfo.get("unpaired")

        # set the correct name for 'func'
        self._instructions["func_name"] = settings.dfa_settings.get_name(self._instructions["func"],
                                                                         self._instructions["prog"])
        self._instructions["disp"] = settings.dfa_settings.get_disp(self._instructions["func"])

        # create/set folder to do the calculations in
        self.folder = os.path.join(self.core.workdir, self.__class__.__name__.lower())
        if os.path.isdir(self.folder):
            print(f"Folder {self.folder} already exists. Potentially overwriting files.")
        elif os.system(f"mkdir {self.folder}") != 0 and not os.path.isdir(self.folder):
            raise RuntimeError(f"Could not create directory for {self.__class__.__name__.lower()}.")

    @timeit
    def run(self) -> None:
        """
        what gets executed if the part is run
        should be implemented for every part respectively
        """
        pass

    def print_info(self) -> None:
        """
        formatted print for part instructions
        """

        # header
        lines = ["\n" + "".ljust(PLENGTH, "-") + "\n",
                 f"{self.__class__.__name__.upper()} - {self.alt_name.upper()}".center(PLENGTH, " ") + "\n",
                 "".ljust(PLENGTH, "-") + "\n", "\n"]

        for instruction, val in self._instructions.items():
            lines.append(f"{instruction}:".ljust(DIGILEN // 2, " ") + f"{val}\n")

        # print everything to console
        for line in lines:
            print(line)
