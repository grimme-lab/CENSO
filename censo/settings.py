import os
import shutil
import configparser
from typing import Any, Dict, Union

from censo.cfg import (
    PROGS,
    GRIDOPTIONS,
    GFNOPTIONS,
    CENSORCNAME,
)
from censo.optimization import *
from censo.part import CensoPart

# map the part names to their respective classes
parts = {
    "general": CensoPart,
    "prescreening": prescreening.Prescreening,
    "screening": screening.Screening,
    "optimization": optimization.Optimization,
}


class CensoRCParser:
    # Available settings are defined here
    __settings_options = {
        "screening": {
            "threshold": {
                "default": 3.5,
                "range": [
                    0.75,
                    7.5
                ]
            },
            "func": {
                "default": "r2scan-3c",
                "options": DfaHelper.find_func("screening")
            },
            "basis": {
                "default": "def2-TZVP",
                "options": BASIS_SETS
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "sm": {
                "default": "smd",
                "options": solv_mods
            },
            "smgsolv": {
                "default": "smd",
                "options": gsolv_mods
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "grid": {
                "default": "low+",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": True
            },
            "gcp": {
                "default": True
            },
            "implicit": {
                "default": True
            }
        },
        "optimization": {
            "optcycles": {
                "default": 8,
                "range": [
                    1,
                    10
                ]
            },
            # "radsize": { # ???
            #     "default": 10,
            #     "range": [
            #         1,
            #         100
            #     ]
            # },
            "maxcyc": {
                "default": 200,
                "range": [
                    10,
                    1000
                ]
            },
            "threshold": {
                "default": 1.5,
                "range": [
                    0.5,
                    5
                ]
            },
            "hlow": {
                "default": 0.01,
                "range": [
                    0.001,
                    0.1
                ]
            },
            "gradthr": {
                "default": 0.01,
                "range": [
                    0.01,
                    1.0
                ]
            },
            "boltzmannthr": {  # boltzmann sum threshold
                "default": 85.0,
                "range": [
                    1.0,
                    99.9
                ]
            },
            # "spearmanthr": {
            #     "default": 0.0,
            #     "range": [
            #         -1.0,
            #         1.0
            #     ]
            # },
            "func": {
                "default": "r2scan-3c",
                "options": dfa_settings.find_func("optimization")
            },
            "basis": {
                "default": "def2-TZVP",
                "options": basis_sets
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "sm": {
                "default": "smd",
                "options": solv_mods
            },
            "smgsolv": {
                "default": "smd",
                "options": gsolv_mods
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "optlevel": {
                "default": "normal",
                "options": [
                    "crude",
                    "sloppy",
                    "loose",
                    "lax",
                    "normal",
                    "tight",
                    "vtight",
                    "extreme",
                ]
            },
            "grid": {
                "default": "high",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": True
            },
            "gcp": {
                "default": True
            },
            "opt_spearman": {
                "default": True
            },
            "crestcheck": {
                "default": False
            }
        },
        # "refinement": {
        #     "threshold": {
        #         "default": 90.0,
        #         "range": [
        #             0.1,
        #             99.9
        #         ]
        #     },
        #     "prog": {
        #         "default": "orca",
        #         "options": PROGS
        #     },
        #     "func": {
        #         "default": "wb97x-v",
        #         "options": dfa_settings.find_func("refinement")
        #     },
        #     "basis": {
        #         "default": "def2-TZVPP",
        #         "options": basis_sets
        #     },
        #     "smgsolv": {
        #         "default": "smd",
        #         "options": gsolv_mods
        #     },
        #     "gfnv": {
        #         "default": "gfn2",
        #         "options": GFNOPTIONS
        #     },
        #     "grid": {
        #         "default": "high+",
        #         "options": GRIDOPTIONS
        #     },
        #     "run": {
        #         "default": False
        #     },
        #     "gcp": {
        #         "default": True
        #     }
        # },
        # "nmr": {
        #     "resonance_frequency": {
        #         "default": 300.0,
        #         "range": [
        #             150.0,
        #             1000.0
        #         ]
        #     },
        #     "prog4_j": {
        #         "default": "tm",
        #         "options": PROGS
        #     },
        #     "func_j": {
        #         "default": "pbe0-d4",
        #         "options": []
        #     },
        #     "basis_j": {
        #         "default": "def2-TZVP",
        #         "options": basis_sets
        #     },
        #     "sm4_j": {
        #         "default": "smd",
        #         "options": solv_mods
        #     },
        #     "prog4_s": {
        #         "default": "tm",
        #         "options": PROGS
        #     },
        #     "func_s": {
        #         "default": "pbe0-d4",
        #         "options": []
        #     },
        #     "basis_s": {
        #         "default": "def2-TZVP",
        #         "options": basis_sets
        #     },
        #     "sm4_s": {
        #         "default": "smd",
        #         "options": solv_mods
        #     },
        #     "h_ref": {
        #         "default": "TMS",
        #         "options": [
        #             "TMS"
        #         ]
        #     },
        #     "c_ref": {
        #         "default": "TMS",
        #         "options": [
        #             "TMS"
        #         ]
        #     },
        #     "f_ref": {
        #         "default": "CFCl3",
        #         "options": [
        #             "CFCl3"
        #         ]
        #     },
        #     "si_ref": {
        #         "default": "TMS",
        #         "options": [
        #             "TMS"
        #         ]
        #     },
        #     "p_ref": {
        #         "default": "TMP",
        #         "options": [
        #             "TMP",
        #             "PH3"
        #         ]
        #     },
        #     "run": {
        #         "default": False
        #     },
        #     "couplings": {
        #         "default": True
        #     },
        #     "shieldings": {
        #         "default": True
        #     },
        #     "h_active": {
        #         "default": True
        #     },
        #     "c_active": {
        #         "default": True
        #     },
        #     "f_active": {
        #         "default": False
        #     },
        #     "si_active": {
        #         "default": False
        #     },
        #     "p_active": {
        #         "default": False
        #     }
        # },
        # "optrot": {
        #     "func": {
        #         "default": "pbe-d4",
        #         "options": dfa_settings.find_func("optrot")
        #     },
        #     "func_or_scf": {
        #         "default": "r2scan-3c",
        #         "options": []
        #     },
        #     "basis": {
        #         "default": "def2-SVPD",
        #         "options": basis_sets
        #     },
        #     "prog": {
        #         "default": "orca",
        #         "options": [
        #             "orca"
        #         ]
        #     },
        #     "run": {
        #         "default": False
        #     },
        #     "freq_or": {
        #         "default": [
        #             598.0
        #         ]
        #     }
        # },
        # "uvvis": {
        #     "nroots": {
        #         "default": 20,
        #         "range": [
        #             1,
        #             100
        #         ]
        #     },
        #     "sigma": {
        #         "default": 0.1,
        #         "range": [
        #             0.1,
        #             1.0
        #         ]
        #     },
        #     "run": {
        #         "default": False
        #     }
        # },
    }

    # try to find the .censorc2 in the user's home directory
    __censorc_path = __find_rcfile()

    # if no configuration file is found, create a new one and configure parts with default settings
    if __censorc_path is None:
        __censorc_path = os.path.join(os.path.expanduser("~"), CENSORCNAME)
        __write_rcfile(__censorc_path)  # TODO
    # otherwise, try read the configuration file and configure the parts with the settings from it
    else:
        __settings = __read_rcfile(__censorc_path)
        global parts
        for section, settings in __settings.items():
            try:
                assert section in parts
                parts[section].set_settings(settings)
            except AssertionError:
                pass

    @staticmethod
    def __read_rcfile(path) -> Dict[str, Dict[str, Any]]:
        """
        Read from config data from file located at 'path'
        """
        rcdata: Dict = {}

        # read config file
        parser: configparser.ConfigParser = configparser.ConfigParser()
        with open(path, "r") as file:
            parser.read_file(file)

        return rcdata

    @staticmethod
    def __write_rcfile(path: str) -> None:
        """
        write new configuration file with default settings into file at 'path' 
        """
        # what to do if there is an existing configuration file
        external_paths = None
        if os.path.isfile(path):
            print(
                f"An existing configuration file has been found at {path}.\n",
                f"Renaming existing file to {CENSORCNAME}_OLD.\n"
            )
            # Read program paths from the existing configuration file
            print("Reading program paths from existing configuration file ...")
            external_paths = __read_program_paths(path)

            # Rename existing file
            os.rename(path, f"{path}_OLD")

        print(f"Writing new configuration file to {path} ...")
        with open(path, "w", newline=None) as rcfile:
            parser = configparser.ConfigParser()
            parser.read_dict(
                {part: {setting: settings[setting]['default'] for setting in settings} for part, settings in
                 self._settings_options.items()})

            # Try to get paths from 'which'
            if external_paths is None:
                print("Trying to determine program paths automatically ...")
                external_paths = self.__find_program_paths()

            parser["paths"] = external_paths
            parser.write(rcfile)

        print(
            f"\nA new configuration file was written into {path}.\n"
            "You should adjust the settings to your needs and set the program paths.\n"
            "Right now only the default settings are used.\n"
        )

        if CENSORCNAME not in path:
            print(
                f"Additionally make sure that the file name is '{CENSORCNAME}'.\n"
                f"Currently it is '{os.path.split(path)[-1]}'.\n"
            )

    @staticmethod
    def __find_rcfile() -> Union[str, None]:
        """
        check for existing .censorc in $home dir
        """

        rcpath = None
        # check for .censorc in $home
        if os.path.isfile(os.path.join(os.path.expanduser("~"), CENSORCNAME)):
            rcpath = os.path.join(os.path.expanduser("~"), CENSORCNAME)

        return rcpath

    @classmethod
    def __read_program_paths(self, path: str) -> Union[Dict[str, str], None]:
        """
        Read program paths from the configuration file at 'path'
        """
        with open(path, "r") as inp:
            parser = configparser.ConfigParser()
            parser.read_file(inp)

        try:
            return parser["paths"]
        except KeyError:
            print(f"WARNING: No paths found in {path}")
            return None

    @classmethod
    def __find_program_paths(self) -> Dict[str, str]:
        """
        Try to determine program paths automatically
        """
        # TODO - for now only the most important ones are implemented
        mapping = {
            "orcapath": "orca",
            "xtbpath": "xtb",
            "crestpath": "crest",
            "cosmorssetup": None,
            "dbpath": None,
            "cosmothermversion": None,
            "mpshiftpath": None,
            "escfpath": None,
        }
        paths = {}

        for pathname, program in mapping.items():
            if program is not None:
                path = shutil.which(program)
            else:
                path = None

            if path is not None:
                paths[pathname] = path
            else:
                paths[pathname] = ""

        # if orca was found try to determine orca version from the path (kinda hacky)
        if paths["orcapath"] != "":
            try:
                paths["orcaversion"] = paths["orcapath"].split(os.sep)[-2][5:10].replace("_", ".")
            except Exception:
                paths["orcaversion"] = ""

        return paths