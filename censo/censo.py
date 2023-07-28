"""
CENSO run code:

"""
from os import getcwd
import os
import shutil
from time import perf_counter
import sys
from traceback import print_exc

from censo.cfg import PLENGTH, DESCR, ASSETS_PATH, __version__
from censo.inputhandling import cml
from censo.utilities import print
from censo.tutorial import interactiv_doc
from censo.core import CensoCore
from censo.prescreening import Prescreening
from censo.settings import CensoSettings

# use generators for reduced memory usage?
# dict.setdefault()
# join dicts with merged_dict = {**d1, **d2}
# remove all mutable defaults in functions
# TODO - MAJOR - fix compatibility with old json and censorc files
# TODO - MAJOR - introduce option to return all user customizable dbs to default
# TODO - more stringent folder naming (not part1, part2, part3, nmr..., etc.)
# TODO - rename functions, change respective variables and assure backwards compatibility
# TODO - reuse gbw files for orca calculations
# TODO - test performance of orca vs tm (same results?)
# TODO - improve thread management
# TODO - give meaningful help messages e.g. to avoid pitfalls for repeated calculations (too high threshold etc.)
# TODO - ask if CENSO should do an automatic cleanup for every run?
# TODO - assign meaning to different return values of main
# TODO - MAJOR - make censo available as package to be easily installed via pip?
# TODO - define custom error types
# TODO - restore coverage of cml args and settings_options
# TODO - error handling
# TODO - function for user input?
# TODO - add option to read out only specific conformers from input (e.g. 1, 3, 4)
# TODO - introduce uniform formatting for print (utilities print redundant?)
# TODO - fix all paths
# TODO - output data in an easily processable format
# TODO - print formatting
def main(argv=None):
    """
    Execute the CENSO code.
    """

    # put every step for startup into a method for convenience
    # makes testing easier and may also be called for customized workflows
    # standard censo setup
    core, settings = startup(argv)

    ### default: part1 - 3
    # TODO - reduce copy/paste code with list of functions which is iterated over
    run = [Prescreening, ]

    for part in run:
        print(f"Ran {part.name} in {part(core, settings).run()} seconds!")
        
    print("\nCENSO all done!")
    return 0


# sets up a core and settings object for you using the given cml arguments and censorc
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
    elif args.tutorial:
        interactiv_doc()
        sys.exit(0)
    elif args.cleanup:
        cleanup_run(cwd)
        print("Removed files and going to exit!")
        sys.exit(0)
    elif args.cleanup_all:
        cleanup_run(cwd, complete=True)
        print("Removed files and going to exit!")
        sys.exit(0)
        
    # setup internal settings with default values
    settings = CensoSettings()
    
    # look for censorc
    settings.find_rcfile(cwd, args.inprcpath)
    
    # TODO - where to put this?
    # no censorc found at standard dest./given dest.
    if settings.censorc_path == "":
        print(
            f"No rcfile has been found. Do you want to create a new one?\n"
        )

        user_input = ""
        while user_input.strip().lower() not in ["yes", "y", "no", "n"]:
            print("Please type 'yes/y' or 'no/n':")
            user_input = input()
        
        if user_input.strip().lower() in ("y", "yes"):
            settings.censorc_path = settings.write_config(args, cwd)
        elif user_input.strip().lower() in ("n", "no"):
            print(
                "Configuration file needed to run CENSO!\n"
                "Going to exit!"
            )
            sys.exit(1)
    
    # read paths for external programs (definition in rcfile?) ????
    settings.read_program_paths()
    
    # update the settings with rcdata and cml args
    settings.settings_current = args
    
    # initialize core with args and cwd
    core = CensoCore(args, cwd)
    
    # read input and setup conformers
    core.read_input()

    ### END of setup
    # -> core.conformers contains all conformers with their info from input (sorted)
    # -> core.ensembledata is blank
    
    return core, settings


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

        print(f"Be aware that files in {cwd} and subdirectories with names containing the following substrings will be deleted:")
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