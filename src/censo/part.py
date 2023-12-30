import functools
import json
import os
import ast
from collections.abc import Callable

from .ensembledata import EnsembleData
from .params import (
    PLENGTH,
    DIGILEN,
    OMPMIN,
    OMPMAX,
    SOLVENTS_DB,
)
from .datastructure import ParallelJob
from .utilities import setup_logger, DfaHelper

logger = setup_logger(__name__)


class CensoPart:
    """
    Part class as parent class for all parts of the calculation.

    Settings that are not part of the general settings and strictly necessary for ...
    ... DFT calculations:
        - prog
        - func
        - basis
        - grid
        - template
        - gcp
    ... xtb calculations:
        - gfnv
    Further requirements are defined in dictionaries belonging to processor classes.
    """

    _options = {
        "maxcores": {"default": 4, "range": [OMPMIN, 256]},
        "omp": {"default": 4, "range": [OMPMIN, OMPMAX]},
        "imagthr": {"default": -100.0, "range": [-300.0, 0.0]},
        "sthr": {"default": 0.0, "range": [0.0, 100.0]},
        "scale": {"default": 1.0, "range": [0.0, 1.0]},
        "temperature": {"default": 298.15, "range": [1e-05, 2000.0]},
        "solvent": {"default": "h2o", "options": [k for k in SOLVENTS_DB.keys()]},
        "sm_rrho": {"default": "alpb", "options": ["alpb", "gbsa"]},
        "multitemp": {"default": True},
        "evaluate_rrho": {"default": True},
        "consider_sym": {"default": True},
        "bhess": {"default": True},
        "rmsdbias": {"default": False},
        "balance": {"default": True},
        "gas-phase": {"default": False},
        "copy_mo": {"default": True},
        "retry_failed": {"default": True},
        "trange": {"default": [273.15, 373.15, 5]},
    }

    _settings = {}

    @staticmethod
    def set_general_settings(settings: dict[str, any]) -> None:
        CensoPart._validate(settings)
        settings = CensoPart._complete(settings)
        CensoPart._settings = settings

    @staticmethod
    def set_general_setting(setting, value):
        assert type(value) is type(CensoPart._settings[setting])
        CensoPart._settings[setting] = value

    @staticmethod
    def get_general_settings():
        return CensoPart._settings

    @classmethod
    def get_settings(cls):
        return cls._settings

    @classmethod
    def set_settings(cls, settings: dict[str, any]):
        cls._validate(settings)
        settings = cls._complete(settings)
        cls._settings = settings

    @classmethod
    def set_setting(cls, setting_name: str, setting_value: any):
        assert type(setting_value) is type(cls._settings[setting_name])
        cls._settings[setting_name] = setting_value

    @classmethod
    def get_options(cls):
        return cls._options

    @classmethod
    def _complete(cls, tocomplete: dict[str, any]) -> dict[str, any]:
        """
        fill in missing settings with default values
        """
        for setting, definition in cls._options.items():
            if setting not in tocomplete:
                tocomplete[setting] = definition["default"]

        return tocomplete

    @classmethod
    def _validate(cls, tovalidate: dict[str, any]) -> None:
        """
        Validates the type of each setting in the given dict. Also potentially validate if the setting is allowed by 
        checking with cls._options.

        Args:
            tovalidate (dict[str, any]): The dict containing the settings to be validated.

        Returns:
            None

        Raises:
            ValueError: If the setting is not allowed or the value is not within the allowed range.
        """
        # go through each section and try to validate each setting's type
        remove = []
        for setting_name in tovalidate:
            # try to get the settings's type from the default value in the _options dict
            try:
                setting_type = type(cls._options[setting_name]["default"])
            except KeyError:
                # KeyError means that the setting does not exist, therefore the setting key should be removed
                remove.append(setting_name)
                continue

            # try to cast the setting-string into the correct type
            try:
                if setting_type == bool:
                    setting_value = {"True": True, "False": False}.get(
                        tovalidate[setting_name]
                    )
                elif setting_type == list:
                    setting_value = ast.literal_eval(tovalidate[setting_name])
                else:
                    setting_value = setting_type(tovalidate[setting_name])
            # if that's not possible raise an exception
            # NOTE: KeyError is raised when the conversion for bools fails
            except (ValueError, KeyError):
                raise ValueError(
                    f"Value '{tovalidate[setting_name]}' is not allowed for setting '{setting_name}' in part of type '{cls.__name__}'"
                )

            # now check if the setting is allowed
            # for strings check if string is within a list of allowed values
            if setting_type == str:
                options = cls._options[setting_name]["options"]
                if setting_value not in options and len(options) > 0:
                    # Only check if there are options
                    # This is fatal so an exception is raised
                    raise ValueError(
                        f"Value '{setting_value}' is not allowed for setting '{setting_name}' in part of type '{cls.__name__}'."
                    )
            # for numeric values check if value is within a range
            elif setting_type in (int, float):
                interval = cls._options[setting_name]["range"]
                if not interval[0] <= setting_value <= interval[1]:
                    # This is fatal so an exception is raised
                    raise ValueError(
                        f"Value '{setting_value}' is out of range "
                        f"({interval[0]},{interval[1]}) for setting '{setting_name}' in part of type '{cls.__name__}'."
                    )
            # NOTE: there is no check for complex types yet (i.e. lists)

            # set the value in the dict tovalidate to the casted value
            tovalidate[setting_name] = setting_value

        # remove the invalid settings
        for setting in remove:
            tovalidate.pop(setting)

    @staticmethod
    def _create_dir(runner: Callable) -> Callable:
        """
        This method needs to be defined as @staticmethod to be accessible from within the class via the @_create_dir
        decorator. The wrapper function within will be able to access the instance variables of the class.
        To access this method from child classes, the decorator must be called like: @CensoPart._create_dir.

        Args:
            runner: The function to be decorated.

        Returns:
            Callable: The decorated function.
        """

        @functools.wraps(runner)
        def wrapper(self, *args, **kwargs):
            # create/set folder to do the calculations in
            self.dir = os.path.join(self.ensemble.workdir, self._name)
            if os.path.isdir(self.dir):
                global logger
                logger.warning(
                    f"Folder {self.dir} already exists. Potentially overwriting files."
                )
            elif os.system(f"mkdir {self.dir}") != 0 and not os.path.isdir(self.dir):
                raise RuntimeError(
                    f"Could not create directory for {self._name}.")

            return runner(self, *args, **kwargs)

        return wrapper

    def __init__(self, ensemble: EnsembleData):
        """
        Initializes a part instance.

        Args:
            ensemble: The ensemble instance that manages the conformers.

        Returns:
            None
        """
        # sets the name of the part (used for printing and folder creation)
        self._name: str = self.__class__.__name__.lower()

        # every part instance depends on a ensemble instance to manage the conformers
        self.ensemble: EnsembleData = ensemble

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

        # Print header with part name
        lines = [
            "\n" + "".ljust(PLENGTH, "-") + "\n",
            f"{self.__class__.__name__.upper()} - {self._name.upper()}".center(
                PLENGTH, " "
            )
            + "\n",
            "".ljust(PLENGTH, "-") + "\n",
            "\n",
        ]

        # Print all settings with name and value
        for setting, val in self._settings.items():
            lines.append(f"{setting}:".ljust(DIGILEN // 2, " ") + f"{val}\n")

        # print everything to console
        for line in lines:
            print(line)

    def write_json(self) -> None:
        """
        Writes the part's results to a json file.

        Returns:
            None
        """
        results = {conf.name: conf.results[self._name]
                   for conf in self.ensemble.conformers}
        with open(
            os.path.join(self.ensemble.workdir, f"{self._name}.json"), "w"
        ) as outfile:
            json.dump(results, outfile, indent=4)
