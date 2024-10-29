import functools
import json
import os
import ast
from math import exp
from collections.abc import Callable

from .ensembledata import EnsembleData
from .params import Params
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
                for sm in functools.reduce(
                    lambda x, y: x + y, Params.SOLV_MODS.values()
                )
            },
        },
        "sm_rrho": {"default": "alpb", "options": Params.SOLV_MODS["xtb"]},
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

    # Maps part names onto numbers
    _part_nos = {
        "prescreening": 0,
        "screening": 1,
        "optimization": 2,
        "refinement": 3,
        "nmr": 4,
        "optrot": 5,
        "uvvis": 6,
    }

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
    def set_settings(cls, settings: dict[str, any], complete: bool = False):
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

        def extract_options(d) -> set:
            """
            Helper function to extract all unique options from the lowest nesting level.
            """
            options = set()

            def extract(dd) -> None:
                if isinstance(dd, dict):
                    for value in dd.values():
                        extract(value)
                else:
                    options.update(dd)

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
                os.getcwd(),
                f"{self._part_nos[self._name]}_{self._name.upper()}",
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
        Initializes a part instance. Also calls the new instance in a single step.
        Yields the instance after execution of itself.

        Args:
            ensemble: The ensemble instance that manages the conformers.

        Returns:
            None
        """
        # sets the name of the part (used for printing and folder creation)
        self._name: str = (
            self.__class__.__name__.lower()
            if self.__class__ is not CensoPart
            else "general settings"
        )

        # every part instance depends on a ensemble instance to manage the conformers
        self.ensemble: EnsembleData = ensemble

        # Directory where the part executes its calculations
        # It is set using the _create_dir method (intended to be used as wrapper for __call__)
        self.dir: str = None

        # Store the part's results
        self.results = {}
        # should be structured like the following:
        # 'CONFXX': <results from part jobs/in-part-calculations>
        # => e.g. self.results["CONF2"]["gtot"]
        #    would return the free enthalpy of the conformer CONF2
        #    (not calculated with an external program)
        #
        #    self.results["CONF3"]["sp"]
        #    returns the 'result' of a DFT single point
        #    (calculated by external program)
        #    to get the single-point energy: self.results["CONF3"]["sp"]["energy"]
        #    (refer to the results for each jobtype)

        self.runtime = self()

        # Attach this part's results to the ensemble results
        self.ensemble.results[self._name].append(self.results)

    def __call__(self) -> None:
        """
        what gets executed if the part is run
        should be implemented for every part respectively
        """
        pass

    def _calc_boltzmannweights(self) -> None:
        """
        Calculate populations for boltzmann distribution of ensemble at given
        temperature given values for free enthalpy.
        """
        temp = self.get_general_settings()["temperature"]
        # find lowest gtot value
        if all(
            [
                "gtot" in self.results[conf.name].keys()
                for conf in self.ensemble.conformers
            ]
        ):
            minfree: float = min(
                [self.results[conf.name]["gtot"] for conf in self.ensemble.conformers]
            )

            # calculate boltzmann factors
            bmfactors = {
                conf.name: conf.degen
                * exp(
                    -(self.results[conf.name]["gtot"] - minfree)
                    * Params.AU2J
                    / (Params.KB * temp)
                )
                for conf in self.ensemble.conformers
            }
        else:
            # NOTE: if anything went wrong in the single-point calculation ("success": False),
            # this should be handled before coming to this step
            # since then the energy might be 'None'
            gtot_replacement = False
            for jt in ["xtb_opt", "sp"]:
                if all(
                    jt in self.results[conf.name].keys()
                    for conf in self.ensemble.conformers
                ):
                    minfree: float = min(
                        self.results[conf.name][jt]["energy"]
                        for conf in self.ensemble.conformers
                    )

                    # calculate boltzmann factors
                    bmfactors = {
                        conf.name: conf.degen
                        * exp(
                            -(self.results[conf.name][jt]["energy"] - minfree)
                            * Params.AU2J
                            / (Params.KB * temp)
                        )
                        for conf in self.ensemble.conformers
                    }
                    gtot_replacement = True
                    break

            if not gtot_replacement:
                raise RuntimeError(
                    f"Could not determine Boltzmann factors for {self._name}."
                )

        # calculate partition function from boltzmann factors
        bsum: float = sum(bmfactors.values())

        # Return Boltzmann populations
        for conf in self.ensemble.conformers:
            self.results[conf.name]["bmw"] = bmfactors[conf.name] / bsum

    def _print_info(self) -> None:
        """
        formatted print for part instructions
        """

        # Print header with part name
        print(h2(f"{self._name.upper()}"))

        # Print all settings with name and value
        for setting, val in self._settings.items():
            print(f"{setting}:".ljust(Params.DIGILEN // 2, " ") + f"{val}")

        print("\n")

    def _write_json(self) -> None:
        """
        Writes the part's results to a json file.

        Returns:
            None
        """
        results = {
            self._name: {
                conf.name: self.results[conf.name] for conf in self.ensemble.conformers
            }
        }
        filename = f"{self._part_nos.get(self._name, 'NaN')}_{self._name.upper()}.json"
        with open(os.path.join(os.getcwd(), filename), "w") as outfile:
            json.dump(results, outfile, indent=4)
