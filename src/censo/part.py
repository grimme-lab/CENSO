import functools
import json
import os
import ast
from collections.abc import Callable

from .ensembledata import EnsembleData
from .params import PLENGTH, DIGILEN, OMPMIN, OMPMAX, SOLV_MODS
from .logging import setup_logger
from .utilities import print, h1, h2, SolventHelper

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
        "omp": {"default": 4},
        "imagthr": {"default": -100.0},
        "sthr": {"default": 0.0},
        "scale": {"default": 1.0},
        "temperature": {"default": 298.15},
        "solvent": {
            "default": "h2o",
            "options": {
                sm: list(SolventHelper.get_solvents_dict(sm).keys())
                for sm in functools.reduce(lambda x, y: x + y, SOLV_MODS.values())
            },
        },
        "sm_rrho": {"default": "alpb", "options": SOLV_MODS["xtb"]},
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

    # this should contain the part number as a string
    # in this case it's just a placeholder, for e.g. prescreening it would be "0"
    _part_no = "NaN"

    @staticmethod
    def set_general_settings(settings: dict[str, any], complete: bool = True) -> None:
        """
        Set all general settings according to a settings dictionary. Will validate the dictionary and complete it
        if complete = True.

        Args:
            settings (dict[str, any]): The settings to be set.
            complete (bool): If True, the settings will be completed with default values if they are missing.

        Returns:
            None
        """
        if complete:
            settings = CensoPart._complete(settings)
        CensoPart._settings.update(settings)
        CensoPart._validate(CensoPart._settings)

    @staticmethod
    def set_general_setting(setting: str, value: any):
        """
        Set a general setting to a specific value. Will check the type of the setting.

        Args:
            setting (str): The setting to be set.
            value (any): The value to be set.

        Returns:
            None
        """
        assert type(value) is type(CensoPart._settings[setting])
        CensoPart._settings[setting] = value

    @staticmethod
    def get_general_settings():
        return CensoPart._settings

    @classmethod
    def get_settings(cls):
        return cls._settings

    @classmethod
    def set_settings(cls, settings: dict[str, any], complete: bool = True):
        """
        Set all part specific settings according to a settings dictionary. Will validate the dictionary and complete
        it if complete = True.

        Args:
            settings (dict[str, any]): The settings to be set.
            complete (bool): If True, the settings will be completed with default values if they are missing.

        Returns:
            None
        """
        if complete:
            settings = cls._complete(settings)
        cls._settings.update(settings)
        cls._validate(cls._settings)

    @classmethod
    def set_setting(cls, setting_name: str, setting_value: any):
        """
        Set a part specific setting to a specific value. Will check the value and type of the setting.

        Args:
            setting_name (str): The setting to be set.
            setting_value (any): The value to be set.

        Returns:
            None
        """
        # Case insensitive check
        setting_name = setting_name.lower()
        if type(setting_value) is str:
            setting_value = setting_value.lower()

        assert type(setting_value) is type(cls._settings[setting_name])
        cls._settings[setting_name] = setting_value
        cls._validate(cls._settings)

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
            ValueError: If the setting is not allowed or the value is not within the allowed options.
        """

        def extract_options(d: dict) -> set:
            """
            Helper function to extract all unique options from the lowest nesting level.
            """
            options = set()

            def extract(dd: dict) -> None:
                for value in dd.values():
                    if isinstance(value, dict):
                        extract(value)
                    else:
                        options.add(value)

            extract(d)
            return options

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

            # If necessary, try to cast the setting-string into the correct type
            if not isinstance(tovalidate[setting_name], setting_type):
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
                except (ValueError, KeyError) as e:
                    raise ValueError(
                        f"Value '{tovalidate[setting_name]}' is not allowed "
                        + f"for setting '{setting_name}' in part of type '{cls.__name__}'"
                    ) from e
            else:
                setting_value = tovalidate[setting_name]

            # now check if the setting is allowed
            # for any type check if setting is within a list of allowed values if it exists
            if setting_type == str:
                # Case insensitive
                setting_value = setting_value.lower()

                # Only check if there are options
                if "options" in cls._options[setting_name].keys():
                    options = extract_options(cls._options[setting_name]["options"])
                    if setting_value not in options and len(options) > 0:
                        # This is fatal so an exception is raised
                        raise ValueError(
                            f"Value '{setting_value}' is not allowed for setting "
                            + f"'{setting_name}' in part '{cls.__name__}'."
                        )

            # set the value in the dict tovalidate to the cast value
            tovalidate[setting_name] = setting_value

        # remove the invalid settings
        # FIXME - is this correct?
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
            self.dir = os.path.join(
                self.ensemble.workdir, f"{self._part_no}_{self._name.upper()}"
            )
            if os.path.isdir(self.dir):
                global logger
                # logger.warning(
                #    f"Folder {self.dir} already exists. Potentially overwriting files."
                # )
            elif os.system(f"mkdir {self.dir}") != 0 and not os.path.isdir(self.dir):
                raise RuntimeError(f"Could not create directory for {self._name}.")

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

        # Directory where the part executes it's calculations
        # It is set using the _create_dir method (intended to be used as wrapper for run)
        self.dir: str = None

    def run(self, ncores: int) -> None:
        """
        what gets executed if the part is run
        should be implemented for every part respectively
        """
        raise NotImplementedError

    def print_info(self) -> None:
        """
        formatted print for part instructions
        """

        # Print header with part name
        if self.__class__ == CensoPart:
            print(h2("GENERAL SETTINGS"))
        else:
            print(h2(f"{self._name.upper()}"))

        # Print all settings with name and value
        for setting, val in self._settings.items():
            print(f"{setting}:".ljust(DIGILEN // 2, " ") + f"{val}")

        print("\n")

    def write_json(self) -> None:
        """
        Writes the part's results to a json file.

        Returns:
            None
        """
        results = {
            conf.name: conf.results[self._name] for conf in self.ensemble.conformers
        }
        filename = f"{self._part_no}_{self._name.upper()}.json"
        with open(os.path.join(self.ensemble.workdir, filename), "w") as outfile:
            json.dump(results, outfile, indent=4)
