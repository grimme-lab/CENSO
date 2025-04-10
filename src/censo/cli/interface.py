from logging import DEBUG
import os
import shutil
import sys
from os import getcwd
from argparse import ArgumentError, Namespace
from datetime import timedelta
from typing import cast
from pathlib import Path
from tabulate import tabulate

from ..config import PartsConfig
from ..config.parts_config import GeneralConfig
from ..config.setup import configure, write_rcfile
from .cml_parser import parse
from ..ensembledata import EnsembleData
from ..ensembleopt import prescreening, screening, optimization, refinement
from ..properties import nmr, uvvis
from ..params import __version__, NCORES, OMPMIN
from ..utilities import printf, h1, PLENGTH
from ..logging import setup_logger, set_loglevel

logger = setup_logger(__name__)


def entry_point(argv: list[str] = sys.argv) -> int:
    """
    Console entry point to execute CENSO from the command line.
    """
    try:
        args = parse(argv)
    except ArgumentError as e:
        printf(e.message)
        return 1
    except SystemExit as e:
        return cast(int, e.code)

    if not any(vars(args).values()):
        printf("CENSO needs at least one argument!")
        return 1

    # Print program call
    printf("CALL: " + " ".join(arg for arg in argv))

    try:
        ensemble, parts_config = startup(args)
    except SystemExit as e:
        return cast(int, e.code)

    # Print all active parts settings once
    if logger.level == DEBUG:
        printf(parts_config)
    else:
        for _, config in parts_config:
            if "run" in config.model_fields:
                if config.run:
                    # Revalidate
                    config.model_validate(config)

                    # Print
                    printf(config)
            else:
                printf(config)

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
    for enabled, func in tasks:
        if enabled:
            runtime = func(ensemble, parts_config, cut=cut)
            printf(f"Ran {func.__name__} in {runtime:.2f} seconds!")
            time += runtime

            if func in (prescreening, screening, optimization, refinement):
                # Collect results for comparison
                mingtot = min(conf.gtot for conf in ensemble)
                comparison[func.__name__] = {
                    conf.name: conf.gtot - mingtot for conf in ensemble
                }

    # Print final comparison
    print_comparison(comparison)

    time = timedelta(seconds=int(time))
    hours, r = divmod(time.seconds, 3600)
    minutes, seconds = divmod(r, 60)
    if time.days:
        hours += time.days * 24

    printf(f"\nRan CENSO in {hours:02d}:{minutes:02d}:{seconds:02d}")

    printf("\nCENSO all done!")
    return 0


# sets up a ensemble object using the given cml arguments and censorc
def startup(args) -> tuple[EnsembleData, PartsConfig]:
    # get most important infos for current run
    cwd = getcwd()

    # run actions for which no complete setup is needed
    if args.version:
        printf(__version__)
        sys.exit()
    elif args.cleanup:
        cleanup_run(cwd)
        printf("Removed files and going to exit!")
        sys.exit()
    elif args.cleanup_all:
        cleanup_run(cwd, complete=True)
        printf("Removed files and going to exit!")
        sys.exit()
    elif args.writeconfig:
        write_rcfile(Path() / "censo2rc_NEW")
        sys.exit()

    parts_config = configure(rcpath=args.inprcpath, args=args)

    if args.loglevel:
        set_loglevel(args.loglevel)

    # initialize ensemble
    ensemble = EnsembleData()

    # read input and setup conformers
    ensemble.read_input(args.inp, args.charge or 0, args.unpaired or 0, args.nconf)

    # if data should be reloaded, do it here
    if args.reload:
        for filename in args.reload:
            ensemble.read_output(os.path.join(cwd, filename))

    if args.maxcores:
        NCORES = args.maxcores

    if args.ompmin:
        OMPMIN = args.ompmin

    # if data should be reloaded, do it here
    if args.reload:
        for filename in args.reload:
            ensemble.read_output(os.path.join(cwd, filename))

    # END of setup
    # -> ensemble.conformers contains all conformers with their info from input (sorted by CREST energy if possible)
    # -> settings are updated with cml args
    # -> output data is reloaded if wanted

    return ensemble, parts_config


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
    printf(h1(f"FINAL RANKING COMPARISON"))

    headers = ["CONF#"]

    headers.extend([f"ΔGtot {part.capitalize()}" for part in comparison])

    units = [
        "",
    ]

    units.extend(["[kcal/mol]" for _ in headers[1:]])

    printmap = {"CONF#": lambda column: [confname for confname in column]}

    for header in headers[1:]:
        printmap[header] = lambda column: [value for value in column.values()]

    rows = [
        printmap[header](column) for header, column in zip(headers, comparison.values())
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
