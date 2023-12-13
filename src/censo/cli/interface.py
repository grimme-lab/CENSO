import os
import shutil
import sys
from os import getcwd

from .cml_parser import parse
from ..configuration import configure
from ..core import CensoCore
from ..ensembleopt import Prescreening, Screening, Optimization
from ..params import DESCR, __version__
from ..utilities import print, setup_logger

logger = setup_logger(__name__)


def entry_point(argv: list[str] | None = None) -> int:
    """
    Console entry point to execute CENSO from the command line.
    """
    # put every step for startup into a method for convenience
    # makes testing easier and may also be called for customized workflows
    # standard censo setup
    try:
        args = parse(DESCR, argv)
    except SystemExit:
        return 0

    if not any(vars(args).values()):
        print("CENSO needs at least one argument!")
        return 1

    core = startup(args)
    if core is None:
        return 0

    run = filter(
        lambda x: x.get_settings()["run"],
        [Prescreening, Screening, Optimization]
    )

    for part in run:
        tmp = part(core)
        print(f"Ran {tmp._name} in {tmp.run()} seconds!")

    print("\nCENSO all done!")
    return 0


# sets up a core object for you using the given cml arguments and censorc
def startup(args) -> CensoCore | None:
    # get most important infos for current run
    cwd = getcwd()

    # run actions for which no complete setup is needed
    if args.version:
        print(__version__)
        return None
    elif args.cleanup:
        cleanup_run(cwd)
        print("Removed files and going to exit!")
        return None
    elif args.cleanup_all:
        cleanup_run(cwd, complete=True)
        print("Removed files and going to exit!")
        return None
    elif args.writeconfig:
        configure(rcpath=cwd, create_new=True)
        return None
    elif args.inprcpath is not None:
        configure(args.inprcpath)

    # initialize core, constructor get runinfo from args
    core = CensoCore(cwd, args=args)

    # read input and setup conformers
    core.read_input(args.inp)

    ### END of setup
    # -> core.conformers contains all conformers with their info from input (sorted by CREST energy if possible)

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
        "censo.log",
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

    print(f"Be aware that files in {cwd} and subdirectories with names containing the following substrings "
          f"will be deleted:")
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
            if any([s in file for s in to_delete]) and int(file.split(".")[2]) > 1:
                print(f"Removing: {file}")
                os.remove(os.path.join(subdir, file))
                deleted += os.path.getsize(os.path.join(subdir, file))

    print(f"Removed {deleted / (1024 * 1024): .2f} MB")
