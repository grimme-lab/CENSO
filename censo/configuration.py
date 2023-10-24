import os
import shutil
import configparser
from typing import Any, Dict, Union

from censo.params import (
    CENSORCNAME,
    load_dbs, ASSETS_PATH
)
from censo.qm_processor import QmProc
from censo.utilities import DfaHelper


def configure(rcpath: str = None):
    """
    Configures the application based on the provided configuration file path.
    If no configuration file path is provided, it searches for the default configuration file.
    If no configuration file is found, it creates a new one with default settings.
    """

    # Try to find the .censorc2 in the user's home directory if no configuration file path is provided
    if rcpath is None:
        censorc_path = find_rcfile()
    else:
        if not os.path.isfile(rcpath):
            raise FileNotFoundError(f"No configuration file found at {rcpath}.")
        else:
            censorc_path = rcpath

    # Set up the DFAHelper
    DfaHelper.set_dfa_dict(os.path.join(ASSETS_PATH, "censo_dfa_settings.json"))

    # Load the lookup tables from the assets directory
    load_dbs()

    # map the part names to their respective classes
    # NOTE: the DFAHelper and the databases should be setup before the parts are imported,
    # otherwise there will be errors in the CensoPart._options
    from censo.part import CensoPart
    from censo.ensembleopt import prescreening, screening, optimization
    parts = {
        "general": CensoPart,
        "prescreening": prescreening.Prescreening,
        "screening": screening.Screening,
        "ensembleopt": optimization.Optimization,
    }

    # If no configuration file is found, create a new one and configure parts with default settings
    if censorc_path is None:
        censorc_path = os.path.join(os.path.expanduser("~"), CENSORCNAME)
        write_rcfile(censorc_path)
    # Otherwise, read the configuration file and configure the parts with the settings from it
    else:
        settings_dict = read_rcfile(censorc_path)
        for section, settings in settings_dict.items():
            try:
                assert section in parts
                parts[section].set_settings(settings)
            except AssertionError:
                pass

    # Update the paths for the processors
    paths = read_rcfile(censorc_path)["paths"]
    QmProc.paths = paths


def read_rcfile(path: str) -> Dict[str, Dict[str, Any]]:
    """
    Read from config data from file located at 'path'
    """
    rcdata: Dict = {}

    # read config file
    parser: configparser.ConfigParser = configparser.ConfigParser()
    with open(path, "r") as file:
        parser.read_file(file)

    return rcdata


def write_rcfile(path: str) -> None:
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
        external_paths = read_program_paths(path)

        # Rename existing file
        os.rename(path, f"{path}_OLD")

    with open(path, "w", newline=None) as rcfile:
        parser = configparser.ConfigParser()

        # collect all default settings from parts and feed them into the parser
        global parts
        parser.read_dict({
            partname: {
                settingname: setting["default"] for settingname, setting in part.get_part_options()
            } for partname, part in parts.items()
        })

        # Try to get paths from 'which'
        if external_paths is None:
            print("Trying to determine program paths automatically ...")
            external_paths = find_program_paths()

        parser["paths"] = external_paths

        print(f"Writing new configuration file to {path} ...")
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


def read_program_paths(path: str) -> Union[Dict[str, str], None]:
    """
    Read program paths from the configuration file at 'path'
    """
    with open(path, "r") as inp:
        parser = configparser.ConfigParser()
        parser.read_file(inp)

    try:
        return dict(parser["paths"])
    except KeyError:
        print(f"WARNING: No paths found in {path}")
        return None


def find_program_paths() -> Dict[str, str]:
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


def find_rcfile() -> Union[str, None]:
    """
    check for existing .censorc2 in $home dir
    """

    rcpath = None
    # check for .censorc in $home
    if os.path.isfile(os.path.join(os.path.expanduser("~"), CENSORCNAME)):
        rcpath = os.path.join(os.path.expanduser("~"), CENSORCNAME)

    return rcpath


__settings_options = {
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
