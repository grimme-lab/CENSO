from logging import DEBUG
import os
import shutil
import sys
from argparse import ArgumentError, Namespace
from datetime import timedelta
import traceback
from typing import cast
from pathlib import Path
from tabulate import tabulate

from ..config import PartsConfig
from ..config.setup import configure, write_rcfile
from .cml_parser import parse
from ..ensembledata import EnsembleData
from ..ensembleopt import prescreening, screening, optimization, refinement
from ..properties import nmr, uvvis
from ..params import AU2KCAL, DESCR, __version__, Returncodes
from ..utilities import printf, h1, PLENGTH
from ..logging import set_filehandler, setup_logger, set_loglevel

logger = setup_logger(__name__)


def entry_point(argv: list[str] | None = None) -> Returncodes:
    """
    Console entry point to execute CENSO from the command line.
    """
    try:
        args = parse(argv)
    except ArgumentError as e:
        printf(e.message)
        return Returncodes.ARGUMENT_ERROR
    except SystemExit as e:
        return cast(Returncodes, e.code)

    if not any(vars(args).values()):
        printf("CENSO needs at least one argument!")
        return Returncodes.ARGUMENT_ERROR

    # Startup description
    print(DESCR + "\n")

    # Print program call
    printf("CALL: " + " ".join(arg for arg in sys.argv))

    try:
        ensemble, parts_config = startup(args)
    except SystemExit as e:
        return cast(Returncodes, e.code)

    # Print all active parts settings once
    if logger.level == DEBUG:
        for fieldname, config in parts_config:
            # Revalidate
            setattr(parts_config, fieldname, config.model_validate(config))

            # Print
            printf(config)
    else:
        for fieldname, config in parts_config:
            if "run" in config.model_fields:
                if config.run:
                    # Revalidate
                    setattr(parts_config, fieldname, config.model_validate(config))

                    # Print
                    printf(config)
            else:
                printf(config)

    printf("\n" + "".ljust(int(PLENGTH), "-") + "\n")

    tasks = [
        (parts_config.prescreening.run, prescreening),
        (parts_config.screening.run, screening),
        (parts_config.optimization.run, optimization),
        (parts_config.refinement.run, refinement),
        (parts_config.nmr.run, nmr),
        (parts_config.uvvis.run, uvvis),
    ]

    comparison: dict[str, dict[str, float]] = {}
    cut: bool = not args.keep_all
    time = 0.0
    try:
        for func in [task for enabled, task in tasks[:4] if enabled]:
            runtime = func(ensemble, parts_config, cut=cut)
            printf(f"Ran {func.__name__} in {runtime:.2f} seconds!")
            time += runtime

            # Collect results for comparison
            mingtot = min(conf.gtot for conf in ensemble)
            comparison[func.__name__] = {
                conf.name: (conf.gtot - mingtot) * AU2KCAL for conf in ensemble
            }

        if len(comparison) > 0:
            print("\nFinished ensemble optimization\n")

            # Print final comparison
            print_comparison(comparison)

        if any(enabled for enabled, _ in tasks[4:]):
            print("\nRunning property calculations\n")

        for func in [task for enabled, task in tasks[4:] if enabled]:
            runtime = func(ensemble, parts_config)
            printf(f"Ran {func.__name__} in {runtime:.2f} seconds!")
            time += runtime
    except:
        tb = traceback.format_exc()
        logger.debug(f"Encountered exception:\n{tb}")

        # Save as much data as possible
        printf(
            "Encountered exception. Stopping CENSO and dumping most recent ensemble."
        )
        ensemble.dump_json(Path("CRASH_DUMP.json"))
        ensemble.dump_xyz(Path("CRASH_DUMP.xyz"))
        sys.exit(Returncodes.GENERIC_ERROR)

    time = timedelta(seconds=int(time))
    hours, r = divmod(time.seconds, 3600)
    minutes, seconds = divmod(r, 60)
    if time.days:
        hours += time.days * 24

    printf(f"\nTotal CENSO runtime: {hours:02d}:{minutes:02d}:{seconds:02d}")

    printf("\nCENSO all done!")
    return Returncodes.OK


# sets up a ensemble object using the given cml arguments and censorc
def startup(args: Namespace) -> tuple[EnsembleData, PartsConfig]:
    # get most important infos for current run
    cwd = os.getcwd()

    # run actions for which no complete setup is needed
    if args.version:
        printf(__version__)
        sys.exit(Returncodes.OK)
    elif args.cleanup:
        cleanup_run(cwd)
        printf("Removed files and going to exit!")
        sys.exit(Returncodes.OK)
    elif args.cleanup_all:
        cleanup_run(cwd, complete=True)
        printf("Removed files and going to exit!")
        sys.exit(Returncodes.OK)
    elif args.writeconfig:
        write_rcfile(Path() / "censo2rc_NEW")
        sys.exit(Returncodes.OK)

    parts_config = configure(rcpath=args.inprcpath, args=args)

    # Set up logging
    set_loglevel(args.loglevel)
    logpath = Path(args.logpath or Path(args.inp).parent / "censo.log").resolve()
    set_filehandler(logpath)

    # initialize ensemble
    ensemble = EnsembleData()

    # read input and setup conformers
    try:
        ensemble.read_input(args.inp, args.charge or 0, args.unpaired or 0, args.nconf)
    except FileNotFoundError:
        printf(f"Could not find input file {Path(args.inp).resolve()}.")
        sys.exit(Returncodes.INPUT_NOT_FOUND)

    if args.constraints:
        try:
            ensemble.constraints = Path(args.constraints).read_text()
        except FileNotFoundError:
            printf(
                f"Could not find constraints file {Path(args.constraints).resolve()}."
            )
            sys.exit(Returncodes.CONSTRAINTS_NOT_FOUND)

    # if data should be reloaded, do it here
    if args.reload:
        for filename in args.reload:
            ensemble.read_output(Path(filename))

    if args.maxcores:
        NCORES = args.maxcores

    if args.ompmin:
        OMPMIN = args.ompmin

    # Set working directory to input parent
    # NOTE: this is the most convenient way to solve this for CLI version
    wd = Path(args.inp).parent.resolve()
    printf(f"Setting working directory to {wd}")
    os.chdir(wd)

    # END of setup
    # -> ensemble.conformers contains all conformers with their info from input (sorted by CREST energy if possible)
    # -> settings are updated with cml args
    # -> output data is reloaded if wanted
    # -> working directory is set

    return ensemble, parts_config


def cleanup_run(cwd: str | Path, complete: bool = False):
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
        printf(
            "Removing ALL files generated by previous CENSO runs, including ensembles and output files (*.xyz, *.out, *.json)!"
        )

    printf(
        f"Be aware that files in {cwd} and subdirectories with names containing the following substrings "
        f"will be deleted:"
    )
    for sub in to_delete:
        printf(sub)

    printf("Do you wish to continue?")
    printf("Please type 'yes'/'y' or 'no'/'n':")

    ui = input()

    while ui.strip().lower() not in ["yes", "y", "no", "n"]:
        printf("Please type 'yes'/'y' or 'no'/'n':")
        ui = input()

    if ui.strip().lower() in ["no", "n"]:
        printf("Aborting cleanup!")
        sys.exit(0)

    # iterate over files in cwd and subdirs recursively and remove them if to delete
    for subdir, dirs, files in os.walk(cwd):
        if any(s in subdir for s in to_delete):
            printf(f"Removing: {subdir}")
            shutil.rmtree(subdir)
        for file in files:
            if any(s in file for s in to_delete) and (
                complete
                or all(
                    not file.endswith(pattern) for pattern in [".xyz", ".out", ".json"]
                )
            ):
                printf(f"Removing: {file}")
                os.remove(os.path.join(subdir, file))


def print_comparison(comparison: dict[str, dict[str, float]]):
    if len(comparison) > 1:
        printf(h1(f"FINAL RANKING COMPARISON"))

        headers = ["CONF#"]

        headers.extend([f"ΔGtot {part.capitalize()}" for part in comparison])

        units = [
            "",
        ]

        units.extend(["[kcal/mol]" for _ in headers[1:]])

        confs = [conf for conf in list(comparison.values())[-1]]

        printmap = {"CONF#": lambda confname: confname}

        for header in headers[1:]:
            printmap[header] = (
                lambda confname, h=header: f"{comparison[h.split()[1].lower()][confname]:.2f}"
            )

        rows = [
            [printmap[header](confname) for header in headers] for confname in confs
        ]

        for i in range(len(headers)):
            headers[i] += "\n" + units[i]

        table = tabulate(
            rows,
            headers=headers,
            colalign=["center" for _ in headers],
            disable_numparse=True,
            numalign="decimal",
        )
        print(table, flush=True)

        printf("".ljust(int(PLENGTH), "-") + "\n")
