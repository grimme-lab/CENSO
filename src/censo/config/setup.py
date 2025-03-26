import os
import shutil
import toml
from argparse import Namespace
from typing import Any
from pathlib import Path


from .parts_config import PartsConfig
from ..params import CENSORCNAME
from ..logging import setup_logger
from ..qm import QmProc

logger = setup_logger(__name__)

parts = {}


def configure(rcpath: str | None = None, args: Namespace | None = None) -> PartsConfig:
    """
    Configures the application based on the provided configuration file path.
    If no configuration file path is provided, it searches for the default configuration file.
    If no configuration file is found, it raises a FileNotFoundError.

    Args:
        rcpath (str, optional): Path to the configuration file.
        args (Namespace, optional): Parsed command line arguments. Defaults to None.

    Returns:
        PartsConfig: Configuration instance.
    """
    # Try to find the .censo2rc in the user's home directory
    # if no configuration file path is provided
    if rcpath is None:
        censorc_path = find_rcfile()
    else:
        if not os.path.isfile(rcpath):
            raise FileNotFoundError(f"No configuration file found at {rcpath}.")
        censorc_path = rcpath

    if censorc_path is not None:
        # Read the actual configuration file (located at rcpath if not None, otherwise rcfile in home dir)
        settings_dict = read_rcfile(censorc_path, silent=False)
        paths = settings_dict["paths"]

        # Create configurations
        parts_config = PartsConfig(**settings_dict)
    else:
        # Create default configurations
        parts_config = PartsConfig()

        # Try to find paths
        paths = find_program_paths()

    # Override general settings only for now
    for field in parts_config.general.model_fields:
        setting = getattr(args, field, None)
        if setting is not None:
            parts_config.general.__setattr__(field, setting)

    # Update the paths for the processors
    QmProc.paths.update(paths)

    return parts_config


def read_rcfile(path: str, silent: bool = True) -> dict[str, dict[str, Any]]:
    """
    Read the configuration file at 'path' and return the settings as a dictionary.

    Args:
        path (str): Path to the configuration file.
        silent (bool): If True, no messages will be printed.

    Returns:
        dict[str, dict[str, Any]]: Dictionary containing the settings read from the configuration file.
    """
    # read config file
    if not silent:
        print(f"Reading configuration file from {path}.")

    return toml.loads(Path(path).read_text())


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
            f"Renaming existing file to {CENSORCNAME}_OLD.\n",
        )
        # Read program paths from the existing configuration file
        print("Reading program paths from existing configuration file ...")
        external_paths = read_program_paths(path)

        # Rename existing file
        os.rename(path, f"{path}_OLD")
    else:
        # Try to get paths from 'which'
        print("Trying to determine program paths automatically ...")
        external_paths = find_program_paths()

    with open(path, "w", newline=None) as rcfile:
        # collect all default settings from parts and feed them into the parser
        parts_config = PartsConfig()
        data = parts_config.model_dump(mode="json") | {"paths": external_paths}

        if external_paths is not None:
            data["paths"] = external_paths

        print(f"Writing new configuration file to {path} ...")
        toml.dump(data, rcfile)

    print(
        f"\nA new configuration file was written into {path}.\n"
        + "You should adjust the settings to your needs and set the program paths.\n"
        + "Right now the settings are at their default values.\n"
    )

    if CENSORCNAME not in path:
        print(
            f"To load it automatically make sure that the file name is '{CENSORCNAME}' and it's located in your home directory.\n"
            f"Current name: '{os.path.split(path)[-1]}'.\n"
        )


def read_program_paths(path: str) -> dict[str, str] | None:
    """
    Read program paths from the configuration file at 'path'
    """
    data = toml.loads(Path(path).read_text())

    try:
        return dict(data["paths"])
    except KeyError:
        logger.warning(f"No paths found in {path}")
        return None


def find_program_paths() -> dict[str, str]:
    """
    Try to determine program paths automatically
    """
    programs = ["orca", "xtb", "mpshift", "escf"]
    paths = {}

    for program in programs:
        if program is not None:
            path = shutil.which(program)
        else:
            path = None

        if path is not None:
            paths[program] = path
        else:
            paths[program] = ""

    # if orca was found try to determine orca version from the path (kinda hacky)
    if paths["orca"] != "":
        try:
            paths["orcaversion"] = (
                paths["orca"].split(os.sep)[-2][5:10].replace("_", ".")
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
