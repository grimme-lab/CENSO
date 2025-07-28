import os
import shutil
from argparse import Namespace
from typing import Any
from pathlib import Path
from configparser import ConfigParser
import warnings

from censo.config.paths import PathsConfig


from .parts_config import PartsConfig
from ..params import CENSORCNAME, QmProg
from ..logging import setup_logger

logger = setup_logger(__name__)


def configure(
    rcpath: str | None = None,
    args: Namespace | None = None,
    context: dict[str, list[str]] | None = None,
) -> PartsConfig:
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
    # NOTE: order of priority for rcfile path:
    # - args.inprcfile
    # - home dir
    # - env variable
    if rcpath is None:
        censorc_path = find_rcfile()
    else:
        if not Path(rcpath).resolve().is_file():
            raise FileNotFoundError(f"No configuration file found at {rcpath}.")
        censorc_path = Path(rcpath).resolve()

    if censorc_path is not None:
        # Read the actual configuration file (located at rcpath if not None, otherwise rcfile in home dir)
        settings_dict = read_rcfile(censorc_path, silent=False)

        # Create configurations without solvlent and paths validation for now (will be validated in the end after cml args)
        with warnings.catch_warnings():
            parts_config = PartsConfig.model_validate(settings_dict)
    else:
        # Try to find paths
        paths = find_program_paths()

        # Create default configurations with auto-detected paths (no path or solvent validation)
        parts_config = PartsConfig.model_validate({"paths": paths})

    # Override general settings only for now
    for field in parts_config.general.__class__.model_fields:
        setting = getattr(args, field, None)
        if setting is not None:
            parts_config.general.__setattr__(field, setting)

    # Revalidate after command line args incl. solvents and paths with context
    parts_config = PartsConfig.model_validate(parts_config, context=context)

    return parts_config


def read_rcfile(path: Path, silent: bool = True) -> dict[str, dict[str, Any]]:
    """
    Read the configuration file at 'path' and return the settings as a dictionary.

    Args:
        path (Path): Path to the configuration file.
        silent (bool): If True, no messages will be printed.

    Returns:
        dict[str, dict[str, Any]]: Dictionary containing the settings read from the configuration file.
    """
    # read config file
    if not silent:
        print(f"Reading configuration file from {path}.")

    parser = ConfigParser()
    parser.read_string(path.read_text())

    return {section: dict(parser[section]) for section in parser.sections()}


def write_rcfile(path: Path) -> None:
    """
    Write new configuration file with default settings into file at 'path'.
    Also reads program paths from preexisting configuration file or tries to
    determine the paths automatically.

    Args:
        path (Path): Path to the new configuration file.

    Returns:
        None
    """
    # what to do if there is an existing configuration file
    paths = None
    if path.is_file():
        print(
            f"An existing configuration file has been found at {path}.\n",
            f"Renaming existing file to {path.name}_OLD.\n",
        )
        # Read program paths from the existing configuration file
        print("Reading program paths from existing configuration file ...")
        paths = read_program_paths(path)

        # Rename existing file
        path.rename(f"{path}_OLD")

    with open(path, "w", newline=None) as rcfile:
        parser = ConfigParser()

        # Try to get paths from 'which'
        print("Trying to determine program paths automatically ...")
        paths = find_program_paths()

        # collect all default settings from parts and feed them into the parser
        # Mute warnings
        with warnings.catch_warnings():
            parts_config = PartsConfig.model_validate({"paths": paths})

        parser.read_dict(parts_config.model_dump(mode="json"))

        print(f"Writing new configuration file to {path} ...")
        parser.write(rcfile)

    print(
        f"\nA new configuration file was written into {path}.\n"
        + "You should adjust the settings to your needs and set the program paths.\n"
        + "Right now the settings are at their default values.\n"
    )

    if CENSORCNAME not in path.stem:
        print(
            f"To load it automatically make sure that the file name is '{CENSORCNAME}' and it's located in your home directory.\n"
            f"Current name: '{os.path.split(path)[-1]}'.\n"
        )


def read_program_paths(path: Path) -> dict[str, str]:
    """
    Read program paths from the configuration file at 'path'
    """
    parser = ConfigParser()
    parser.read_string(path.read_text())
    paths = parser.has_section("paths")

    if not paths:
        logger.warning(f"No paths found in {path}")
        return dict()

    return dict(parser["paths"])


def find_program_paths() -> dict[str, str]:
    """
    Try to determine program paths automatically
    """
    paths: dict[str, str] = {}
    for program in PathsConfig.model_fields.keys():
        path = shutil.which(program)

        if path:
            paths[program] = path

    if QmProg.TM.value in PathsConfig.model_fields.keys():
        path = shutil.which("ridft")
        if path:
            paths[QmProg.TM.value] = str(Path(path).parent)

    # If cosmotherm is found try to set cosmorssetup automatically
    if "cosmotherm" in paths:
        # cosmotherm path parent should be BIN-LINUX, CTDATA-FILES is on the same level
        ctdata = (Path(paths["cosmotherm"]).parent / ".." / "CTDATA-FILES").resolve()
        for file in ctdata.glob("*.ctd"):
            if "BP_TZVP" in file.name and file.is_file():
                paths["cosmorssetup"] = file.name

    return paths


def find_rcfile() -> Path | None:
    """
    Check for existing .censorc2 in $home dir or rcfile path in environment variable.
    """

    rcpath = None
    homepath = Path("~").expanduser() / CENSORCNAME
    try:
        envpath = Path(os.environ.get("CENSORC_PATH"))
    except TypeError:
        envpath = None

    # check for .censorc in $home
    if homepath.is_file():
        rcpath = homepath
    # check for path in env variable
    elif envpath is not None and envpath.is_file():
        rcpath = envpath

    return rcpath
