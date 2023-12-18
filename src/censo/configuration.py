import os
import shutil
import configparser

from .params import CENSORCNAME, load_dbs, ASSETS_PATH, USER_ASSETS_PATH
from .qm_processor import QmProc
from .utilities import DfaHelper

parts = {}


def configure(rcpath: str = None, create_new: bool = False):
    """
    Configures the application based on the provided configuration file path.
    If no configuration file path is provided, it searches for the default configuration file.
    If no configuration file is found, it creates a new one with default settings.
    """
    # Try to find the .censo2rc in the user's home directory if no configuration file path is provided
    if rcpath is None:
        censorc_path = find_rcfile()
    else:
        if not os.path.isfile(rcpath) and not create_new:
            raise FileNotFoundError(f"No configuration file found at {rcpath}.")
        else:
            censorc_path = rcpath

    # TODO - Set up the logger

    # Set up the DFAHelper
    DfaHelper.set_dfa_dict(os.path.join(ASSETS_PATH, "censo_dfa_settings.json"))

    # Load the lookup tables from the assets directory
    load_dbs()

    # map the part names to their respective classes
    # NOTE: the DFAHelper and the databases should be setup before the parts are imported,
    # otherwise there will be errors in the CensoPart._options
    from .part import CensoPart
    from .ensembleopt import Prescreening, Screening, Optimization
    from .properties import EnsembleNMR

    global parts
    parts = {
        "prescreening": Prescreening,
        "screening": Screening,
        "optimization": Optimization,
        "nmr": EnsembleNMR,
    }

    # If no configuration file is found, create a new one and configure parts with default settings
    if censorc_path is None:
        censorc_path = os.path.join(os.path.expanduser("~"), CENSORCNAME)
        write_rcfile(censorc_path)
    # if explicitely set to create a new configuration file, do so
    elif create_new:
        censorc_path = os.path.join(rcpath, "censo2rc_NEW")
        write_rcfile(censorc_path)
    # Otherwise, read the configuration file and configure the parts with the settings from it
    else:
        settings_dict = read_rcfile(censorc_path)

        # first set general settings
        CensoPart.set_general_settings(settings_dict["general"])

        # set settings for each part
        for section, settings in settings_dict.items():
            if section in parts.keys():
                parts[section].set_settings(settings)
            # NOTE: if section is not in the parts names, it will be ignored

    # Update the paths for the processors
    paths = read_rcfile(censorc_path)["paths"]
    QmProc._paths.update(paths)

    # create user assets folder if it does not exist
    if not os.path.isdir(USER_ASSETS_PATH):
        os.mkdir(USER_ASSETS_PATH)


def read_rcfile(path: str) -> dict[str, dict[str, any]]:
    """
    Read from config data from file located at 'path'
    """
    # read config file
    parser: configparser.ConfigParser = configparser.ConfigParser()
    with open(path, "r") as file:
        parser.read_file(file)

    returndict = {section: dict(parser[section]) for section in parser.sections()}
    return returndict


def write_rcfile(path: str) -> None:
    """
    write new configuration file with default settings into file at 'path'
    also reads program paths from preexisting configuration file or tries to determine the paths automatically
    """
    # what to do if there is an existing configuration file
    external_paths = None
    if os.path.isfile(path):
        print(
            f"An existing configuration file has been found at {path}.\n",
            f"Renaming existing file to {CENSORCNAME}_OLD.\n",
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
        from .part import CensoPart

        parts["general"] = CensoPart
        parser.read_dict(
            {
                partname: {
                    settingname: setting["default"]
                    for settingname, setting in part.get_options().items()
                }
                for partname, part in parts.items()
            }
        )

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
        "Right now the settings are at their default values.\n"
    )

    if CENSORCNAME not in path:
        print(
            f"Additionally make sure that the file name is '{CENSORCNAME}'.\n"
            f"Currently it is '{os.path.split(path)[-1]}'.\n"
        )


def read_program_paths(path: str) -> dict[str, str] | None:
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


def find_program_paths() -> dict[str, str]:
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
            paths["orcaversion"] = (
                paths["orcapath"].split(os.sep)[-2][5:10].replace("_", ".")
            )
        except Exception:
            paths["orcaversion"] = ""

    return paths


def find_rcfile() -> str | None:
    """
    check for existing .censorc2 in $home dir
    """

    rcpath = None
    # check for .censorc in $home
    if os.path.isfile(os.path.join(os.path.expanduser("~"), CENSORCNAME)):
        rcpath = os.path.join(os.path.expanduser("~"), CENSORCNAME)

    return rcpath


# __settings_options = {
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
# }
