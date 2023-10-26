"""
CENSO run code:

"""
from os import getcwd
import os
import shutil
import sys

from censo import configure
from censo.params import DESCR, ASSETS_PATH, __version__
from censo.inputhandling import cml
from censo.utilities import print
from censo.core import CensoCore
from censo.ensembleopt import *


# TODO - MAJOR - fix compatibility with old json and censorc files
# TODO - MAJOR - introduce option to return all user customizable dbs to default
# TODO - give meaningful help messages e.g. to avoid pitfalls for repeated calculations (too high threshold etc.)
# TODO - MAJOR - make censo available as package to be easily installed via pip?
# TODO - restore coverage of cml args and settings_options
# TODO - function for user input?
# TODO - add option to read out only specific conformers from input (e.g. 1, 3, 4)
def main(argv=None):
    """
    Execute the CENSO code.
    """

    # put every step for startup into a method for convenience
    # makes testing easier and may also be called for customized workflows
    # standard censo setup
    core = startup(argv)

    run = filter(
        lambda x: x.get_settings()["run"],
        [prescreening.Prescreening, screening.Screening, optimization.Optimization]
    )

    for part in run:
        tmp = part(core)
        print(f"Ran {tmp._name} in {tmp.run()} seconds!")

    print("\nCENSO all done!")
    return 0


# sets up a core object for you using the given cml arguments and censorc
def startup(argv):
    # first, check program integrity
    # TODO - proper implementation?
    if not os.path.isdir(ASSETS_PATH):
        raise FileNotFoundError(ASSETS_PATH)

    # get most important infos for current run
    args = cml(DESCR, argv)
    cwd = getcwd()

    # run actions for which no complete setup is needed
    if args.version:
        print(__version__)
        sys.exit(0)
    # elif args.tutorial:
    #     interactiv_doc()
    #     sys.exit(0)
    elif args.cleanup:
        cleanup_run(cwd)
        print("Removed files and going to exit!")
        sys.exit(0)
    elif args.cleanup_all:
        cleanup_run(cwd, complete=True)
        print("Removed files and going to exit!")
        sys.exit(0)
    # TODO - insert option to just create a new censorc

    if args.inprcpath is not None:
        configure(args.inprcpath)

    # initialize core, constructor get runinfo from args
    core = CensoCore(cwd, args=args)

    # read input and setup conformers
    # NOTE: read_input used without keyword arguments because usage of this function implies cml usage of CENSO
    core.read_input(args.inp)

    ### END of setup
    # -> core.conformers contains all conformers with their info from input (sorted by preliminary xtb energy if possible)

    return core


def cleanup_run(cwd, complete=False):
    """
        Delete all unneeded files.
        """

    # files containing these patterns are deleted
    to_delete = [
        "enso.json.",
        "enso_ensemble_part1.xyz.",
        "enso_ensemble_part2.xyz.",
        "enso_ensemble_part3.xyz.",
        "a3mat.tmp",
        "a2mat.tmp",
        "amat.tmp",
    ]

    # remove conformer_rotamer_check folder if complete cleanup
    if complete:
        print("Cleaning up the directory from ALL unneeded files!")
        to_delete[0] = "enso.json"
        if os.path.isdir(os.path.join(cwd, "conformer_rotamer_check")):
            print("Removing conformer_rotamer_check")
            shutil.rmtree(os.path.join(cwd, "conformer_rotamer_check"))
    else:
        print("Cleaning up the directory from unneeded files!")

    print(
        f"Be aware that files in {cwd} and subdirectories with names containing the following substrings will be deleted:")
    for sub in to_delete:
        print(sub)

    print("Do you wish to continue?")
    print("Please type 'yes' or 'no':")

    ui = input()
    if ui.strip().lower() not in ["yes", "y"]:
        print("Aborting cleanup!")
        sys.exit(0)

    # iterate over files in cwd and subdirs recursively and remove them if to delete
    deleted = 0
    for subdir, dirs, files in os.walk(cwd):
        for file in files:
            if any([str in file for str in to_delete]) and int(file.split(".")[2]) > 1:
                print(f"Removing: {file}")
                os.remove(os.path.join(subdir, file))
                deleted += os.path.getsize(os.path.join(subdir, file))

    print(f"Removed {deleted / (1024 * 1024): .2f} MB")
