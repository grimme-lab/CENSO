import os
import shutil
import sys
from os import getcwd
from argparse import ArgumentError
from datetime import timedelta

from .cml_parser import parse
from ..configuration import configure, override_rc
from ..ensembledata import EnsembleData
from ..ensembleopt import Prescreening, Screening, Optimization, Refinement
from ..part import CensoPart
from ..properties import NMR, UVVis
from ..params import __version__, Config
from ..utilities import print
from ..logging import setup_logger, set_loglevel

logger = setup_logger(__name__)


def entry_point(argv: list[str] | None = None) -> int:
    """
    Console entry point to execute CENSO from the command line.
    """
    try:
        args = parse(argv=argv)
    except ArgumentError as e:
        print(e.message)
        return 1
    except SystemExit:
        return 0

    if not any(vars(args).values()):
        print("CENSO needs at least one argument!")
        return 1

    # Print program call
    # FIXME - what happens here? argv is always none, yet it is parsed above without problems?
    print("CALL: " + " ".join(arg for arg in sys.argv))

    ensemble = startup(args)
    if ensemble is None:
        return 0

    # Print general settings once
    CensoPart(ensemble)

    run = filter(
        lambda x: x.get_settings()["run"],
        [Prescreening, Screening, Optimization, Refinement, NMR, UVVis],
    )

    time = 0.0
    for part in run:
        res, runtime = part.run(ensemble)
        print(f"Ran {res.name} in {runtime:.2f} seconds!")
        time += runtime

    time = timedelta(seconds=int(time))
    hours, r = divmod(time.seconds, 3600)
    minutes, seconds = divmod(r, 60)
    if time.days:
        hours += time.days * 24

    print(f"\nRan CENSO in {hours:02d}:{minutes:02d}:{seconds:02d}")

    print("\nCENSO all done!")
    return 0


# sets up a ensemble object for you using the given cml arguments and censorc
def startup(args) -> EnsembleData | None:
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

    if args.loglevel:
        set_loglevel(args.loglevel)

    # Override settings with command line arguments
    override_rc(args)

    # initialize ensemble, constructor get runinfo from args
    ensemble = EnsembleData()

    # read input and setup conformers
    ensemble.read_input(
        args.inp, charge=args.charge, unpaired=args.unpaired, nconf=args.nconf
    )

    # if data should be reloaded, do it here
    if args.reload:
        for filename in args.reload:
            ensemble.read_output(os.path.join(cwd, filename))

    # Set multiprocessing variables
    if args.maxcores:
        Config.NCORES = args.maxcores

    if args.omp:
        Config.OMP = args.omp

    if args.ompmin:
        Config.OMPMIN = args.ompmin

    # if data should be reloaded, do it here
    if args.reload:
        for filename in args.reload:
            ensemble.read_output(os.path.join(cwd, filename))

    # END of setup
    # -> ensemble.conformers contains all conformers with their info from input (sorted by CREST energy if possible)
    # -> output data is reloaded if wanted

    return ensemble


def cleanup_run(cwd, complete=False):
    """
    Delete all unneeded files.
    """

    # files containing these patterns are deleted
    to_delete = [
        "censo.log",
        "0_PRESCREENING",
        "1_SCREENING",
        "2_OPTIMIZATION",
        "3_REFINEMENT",
        "4_NMR",
        "6_UVVIS",
    ]

    if complete:
        print(
            "Removing ALL files generated by previous CENSO runs, including ensembles!"
        )

    print(
        f"Be aware that files in {cwd} and subdirectories with names containing the following substrings "
        f"will be deleted:"
    )
    for sub in to_delete:
        print(sub)

    print("Do you wish to continue?")
    print("Please type 'yes' or 'no':")

    ui = input()
    if ui.strip().lower() not in ["yes", "y"]:
        print("Aborting cleanup!")
        sys.exit(0)

    # iterate over files in cwd and subdirs recursively and remove them if to delete
    for subdir, dirs, files in os.walk(cwd):
        if any(s in subdir for s in to_delete):
            print(f"Removing: {subdir}")
            shutil.rmtree(subdir)
        for file in files:
            if any(s in file for s in to_delete) and (
                complete or "ensemble" not in file
            ):
                print(f"Removing: {file}")
                os.remove(os.path.join(subdir, file))
