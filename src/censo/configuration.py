import os
import shutil
import configparser
from argparse import Namespace

from .params import Config
from .qm_processor import QmProc
from .utilities import DfaHelper, SolventHelper, print

parts = {}


def configure(rcpath: str = None, create_new: bool = False):
    """
    Configures the application based on the provided configuration file path.
    If no configuration file path is provided, it searches for the default configuration file.
    If no configuration file is found, it raises an error.

    Args:
        rcpath (str): Path to the configuration file.
        create_new (bool): If True, a new configuration file will be created at rcpath.

    Returns:
        None
    """
    # Try to find the .censo2rc in the user's home directory
    # if no configuration file path is provided
    if rcpath is None:
        censorc_path = find_rcfile()
    else:
        if not os.path.isfile(rcpath) and not create_new:
            raise FileNotFoundError(f"No configuration file found at {rcpath}.")
        censorc_path = rcpath

    # Set up the DFAHelper
    DfaHelper.set_dfa_dict(os.path.join(Config.ASSETS_PATH, "censo_dfa_settings.json"))

    # Set up the SolventHelper
    SolventHelper.set_solvent_dict(
        os.path.join(Config.ASSETS_PATH, "censo_solvents_db.json")
    )

    # map the part names to their respective classes
    # NOTE: the DFAHelper and the databases should be setup before the parts are imported,
    # otherwise there will be errors in the CensoPart._options
    from .part import CensoPart
    from .ensembleopt import Prescreening, Screening, Optimization, Refinement
    from .properties import NMR, UVVis

    global parts
    parts = {
        "prescreening": Prescreening,
        "screening": Screening,
        "optimization": Optimization,
        "refinement": Refinement,
        "nmr": NMR,
        "uvvis": UVVis,
    }

    # if explicitely told to create a new configuration file, do so
    if create_new:
        if rcpath is None:
            # If not chosen otherwise, the new rcfile is written in the home dir
            censorc_path = os.path.join(os.path.expanduser("~"), "censo2rc_NEW")
        else:
            censorc_path = os.path.join(rcpath, "censo2rc_NEW")
        write_rcfile(censorc_path)
    else:
        # Initialize default settings
        # Make sure that settings are initialized even if there is no section for this part in the rcfile
        # General settings should always be configured first
        CensoPart.set_general_settings({})
        for part in parts.values():
            part.set_settings({}, complete=True)

        # Read rcfile if it exists
        if censorc_path is not None:
            # Read the actual configuration file (located at rcpath if not None, otherwise rcfile in home dir)
            settings_dict = read_rcfile(censorc_path, silent=False)

            # first set general settings
            CensoPart.set_general_settings(settings_dict["general"])

            # Then the remaining settings for each part
            for section, settings in settings_dict.items():
                if section in parts:
                    parts[section].set_settings(settings)
                # NOTE: if section is not in the parts names, it will be ignored

            paths = read_rcfile(censorc_path)["paths"]
        else:
            # Try to automatically determine program paths (not guaranteed to succeed)
            paths = find_program_paths()

        # Update the paths for the processors
        QmProc._paths.update(paths)

    # create user assets folder if it does not exist
    if not os.path.isdir(Config.USER_ASSETS_PATH):
        os.mkdir(Config.USER_ASSETS_PATH)


def read_rcfile(path: str, silent: bool = True) -> dict[str, dict[str, any]]:
    """
    Read the configuration file at 'path' and return the settings as a dictionary.

    Args:
        path (str): Path to the configuration file.
        silent (bool): If True, no messages will be printed.

    Returns:
        dict[str, dict[str, any]]: Dictionary containing the settings read from the configuration file.
    """
    # read config file
    if not silent:
        print(f"Reading configuration file from {path}.")

    parser: configparser.ConfigParser = configparser.ConfigParser()
    with open(path, "r") as file:
        parser.read_file(file)

    returndict = {section: dict(parser[section]) for section in parser.sections()}
    return returndict


def write_rcfile(path: str) -> None:
    """
    Write new configuration file with default settings into file at 'path'.
    Also reads program paths from preexisting configuration file or tries to
    determine the paths automatically.

    Args:
        path (str): Path to the new configuration file.

    Returns:
        None
    """
    # what to do if there is an existing configuration file
    external_paths = None
    if os.path.isfile(path):
        print(
            f"An existing configuration file has been found at {path}.\n",
            f"Renaming existing file to {Config.CENSORCNAME}_OLD.\n",
        )
        # Read program paths from the existing configuration file
        print("Reading program paths from existing configuration file ...")
        external_paths = read_program_paths(path)

        # Rename existing file
        os.rename(path, f"{path}_OLD")

    with open(path, "w", newline=None) as rcfile:
        parser = configparser.ConfigParser()

        # collect all default settings from parts and feed them into the parser
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

    if Config.CENSORCNAME not in path:
        print(
            f"Additionally make sure that the file name is '{Config.CENSORCNAME}'.\n"
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
        "mpshiftpath": "mpshift",
        "escfpath": "escf",
        # "crestpath": "crest",
        # "cosmorssetup": None,
        # "dbpath": None,
        # "cosmothermversion": None,
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
    if os.path.isfile(os.path.join(os.path.expanduser("~"), Config.CENSORCNAME)):
        rcpath = os.path.join(os.path.expanduser("~"), Config.CENSORCNAME)

    return rcpath


def override_rc(args: Namespace) -> None:
    """
    Override the settings from the rcfile (or default settings) with settings from the command line.

    Args:
        args(Namespace): Namespace generated by command line parser.

    Returns:
        None
    """
    # Override general and part specific settings
    from .part import CensoPart

    for part in list(parts.values()) + [CensoPart]:
        part_settings = part.get_settings()
        for setting in part_settings:
            if getattr(args, setting, None) is not None:
                part.set_setting(setting, getattr(args, setting))
