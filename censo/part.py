import functools
from typing import Dict, Any, Callable
import os
import ast

from censo.core import CensoCore
from censo.params import (
    PLENGTH,
    DIGILEN,
    OMPMIN,
    OMPMAX,
    SOLVENTS_DB,
    COSMORS_PARAM,
)
from censo.utilities import (
    DfaHelper
)

"""
Part class as parent class for all parts of the calculation to
implement complete OOP approach.
"""


class CensoPart:

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
                "options": [k for k in SOLVENTS_DB.keys()]
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
                "options": [k for k in COSMORS_PARAM.keys()]
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
        return {**cls._settings, **CensoPart._settings}

    @classmethod
    def get_part_settings(cls):
        # TODO - this is not really that useful yet
        return cls._settings

    @classmethod
    def set_settings(cls, settings: Dict[str, Any]):
        cls._validate(settings)
        settings = cls._complete(settings)
        cls._settings = settings

    @classmethod
    def get_options(cls):
        return {**cls._options, **CensoPart._options}

    @classmethod
    def get_part_options(cls):
        return cls._options

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
            remove = []
            for setting_name in tovalidate[part]:
                # try to get the settings's type from the default value in the _options dict
                try:
                    setting_type = type(cls._options[part][setting_name]["default"])
                except KeyError:
                    # KeyError means that the setting does not exist, therefore the setting is removed
                    remove.append(setting_name)
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

            # remove the invalid settings
            for setting in remove:
                tovalidate[part].pop(setting)

    @staticmethod
    def _create_dir(runner: Callable) -> Callable:
        """
        This method needs to be defined as @staticmethod to be accessible from within the class via the @_create_dir decorator.
        The wrapper function within will be able to access the instance variables of the class.
        To access this method from child classes, the decorator must be called like: @CensoPart._create_dir.
        """
        @functools.wraps(runner)
        def wrapper(self, *args, **kwargs):
            # create/set folder to do the calculations in
            self.dir = os.path.join(self.core.workdir, self._name)
            if os.path.isdir(self.dir):
                print(f"Folder {self.dir} already exists. Potentially overwriting files.")
            elif os.system(f"mkdir {self.dir}") != 0 and not os.path.isdir(self.dir):
                raise RuntimeError(f"Could not create directory for {self._name}.")

            return runner(self, *args, **kwargs)

        return wrapper

    def __init__(self, core: CensoCore, name: str = None):
        if name is None:
            name = self.__class__.__name__.lower()

        # sets the name of the part (used for printing and folder creation)
        self._name: str = name

        # every part instance depends on a core instance to manage the conformers
        self.core: CensoCore = core

        # dictionary with instructions that get passed to the processors
        # basically collapses the first level of nesting into a dict that is not divided into parts
        self._instructions: Dict[str, Any] = \
            {key: value for nested_dict in self.get_settings().values() for key, value in nested_dict.items()}

        # add some additional settings to instructions so that the processors don't have to do any lookups
        # NOTE: [1] auto-selects replacement solvent (TODO - print warning!)
        self._instructions["solvent_key_xtb"] = SOLVENTS_DB.get(self._instructions["solvent"])["xtb"][1]
        if 'sm' in self._instructions.keys():
            self._instructions["solvent_key_prog"] = \
                SOLVENTS_DB.get(self._instructions["solvent"])[self._instructions["sm"]][1]
            # TODO - doesn't work yet for parts where 'func' keyword doesn't exist or there are multiple functionals
        self._instructions["func_type"] = DfaHelper.get_type(self._instructions["func"])

        # add 'charge' and 'unpaired' to instructions
        self._instructions["charge"] = core.runinfo.get("charge")
        self._instructions["unpaired"] = core.runinfo.get("unpaired")

        # set the correct name for 'func'
        self._instructions["func_name"] = DfaHelper.get_name(self._instructions["func"], self._instructions["prog"])
        self._instructions["disp"] = DfaHelper.get_disp(self._instructions["func"])

        self.dir: str = None

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
                 f"{self.__class__.__name__.upper()} - {self._name.upper()}".center(PLENGTH, " ") + "\n",
                 "".ljust(PLENGTH, "-") + "\n", "\n"]

        for instruction, val in self.get_settings().items():
            lines.append(f"{instruction}:".ljust(DIGILEN // 2, " ") + f"{val}\n")

        # print everything to console
        for line in lines:
            print(line)
