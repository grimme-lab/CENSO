import os
import sys
from argparse import ArgumentError, Namespace
from typing import cast, TYPE_CHECKING
from pathlib import Path
import time
from time import sleep

from pydantic import ValidationError

if TYPE_CHECKING:
    from ..ensemble import EnsembleData
    from ..config import PartsConfig

from ..params import DESCR, __version__, Returncode, PLENGTH
from ..utilities import printf
from ..logging import setup_logger, set_loglevel
from .cml_parser import parse
from ..utilities import print_validation_errors, get_time

logger = setup_logger(__name__)


def entry_point(argv: list[str] | None = None) -> Returncode:
    """
    Console entry point to execute CENSO from the command line.

    :param argv: Command line arguments
    :return: Return code
    """
    try:
        args = parse(argv)
    except ArgumentError as e:
        printf(e.message)
        return Returncode.ARGUMENT_ERROR
    except SystemExit as e:
        return cast(Returncode, e.code)

    if not any(vars(args).values()):
        printf("CENSO needs at least one argument!")
        return Returncode.ARGUMENT_ERROR

    # Startup description
    print(DESCR + "\n")

    # Print program call
    printf("CALL: " + " ".join(arg for arg in sys.argv))

    from ..ensembleopt import prescreening, screening, optimization, refinement
    from ..properties import nmr, uvvis, rot

    tasks = [
        ("prescreening", getattr(args, "prescreening", False), prescreening),
        ("screening", getattr(args, "screening", False), screening),
        ("optimization", getattr(args, "optimization", False), optimization),
        ("refinement", getattr(args, "refinement", False), refinement),
        ("nmr", getattr(args, "nmr", False), nmr),
        ("rot", getattr(args, "rot", False), rot),
        ("uvvis", getattr(args, "uvvis", False), uvvis),
    ]

    context = {"check": [name for name, enabled, _ in tasks if enabled]}

    try:
        ensemble, parts_config = startup(args, context)
    except SystemExit as e:
        return cast(Returncode, e.code)

    # Set up parallel managers once for the entire run
    from ..parallel import get_cluster

    try:
        cluster = get_cluster(maxcores=args.maxcores, ompmin=args.ompmin)
        client = cluster.get_client()
    except Exception as e:
        printf(f"Encountered error in setting up parallelization: {e}")
        return Returncode.GENERIC_ERROR

    comparison: dict[str, dict[str, float]] = {}
    cut: bool = not args.keep_all
    times: dict[str, float] = {}

    logger.debug("Setting up parallel managers for CLI run...")
    try:
        for partname, func in [
            (partname, task) for partname, enabled, task in tasks[:4] if enabled
        ]:
            start = time.perf_counter()
            _ = func(ensemble, parts_config, client, cut=cut)  # type: ignore[operator]
            end = time.perf_counter()
            runtime = end - start
            printf(f"Ran {func.__name__} in {runtime:.2f} seconds!")
            times[partname] = runtime

            # Collect results for comparison
            from ..params import AU2KCAL

            mingtot = min(conf.gtot for conf in ensemble)
            comparison[func.__name__] = {
                conf.name: (conf.gtot - mingtot) * AU2KCAL for conf in ensemble
            }

        if len(comparison) > 0:
            print("\nFinished ensemble optimization\n")

            # Print final comparison
            print_comparison(comparison)

        if any(enabled for _, enabled, _ in tasks[4:]):
            print("\nRunning property calculations\n")

        for partname, func in [
            (partname, task) for partname, enabled, task in tasks[4:] if enabled
        ]:
            start = time.perf_counter()
            _ = func(ensemble, parts_config, client)  # type: ignore[operator]
            end = time.perf_counter()
            runtime = end - start
            printf(f"Ran {func.__name__} in {runtime:.2f} seconds!")
            times[partname] = runtime
    except Exception:
        import traceback

        tb = traceback.format_exc()
        logger.debug(f"Encountered exception:\n{tb}")

        # Save as much data as possible
        printf(
            "Encountered exception. Stopping CENSO and dumping most recent ensemble."
        )
        ensemble.dump_json(Path("CRASH_DUMP.json"))
        ensemble.dump_xyz(Path("CRASH_DUMP.xyz"))
        sleep(5)  # Make dask exit more graceful
        return Returncode.GENERIC_ERROR

    total_time = sum(times.values())

    printf("\nTimings:")
    for part in times:
        seconds, minutes, hours = get_time(times[part])
        printf(
            f"{part.capitalize():>19}: {hours:02d}:{minutes:02d}:{seconds:02d}  ({times[part] / total_time * 100:5.1f} %)"
        )

    printf("=" * 40)
    seconds, minutes, hours = get_time(total_time)
    printf(f"\n{'Total CENSO runtime':>19}: {hours:02d}:{minutes:02d}:{seconds:02d}")

    printf("\nCENSO all done!")
    sleep(5)  # Make dask exit more graceful
    return Returncode.OK


# sets up a ensemble object using the given cml arguments and censorc
def startup(
    args: Namespace, context: dict[str, list[str]]
) -> tuple["EnsembleData", "PartsConfig"]:
    """
    Set up the ensemble and configuration using command line arguments and context.

    :param args: Parsed command line arguments
    :param context: Context dictionary with enabled tasks
    :return: Tuple of ensemble data and parts configuration
    """
    from logging import DEBUG

    from ..config.setup import configure, write_rcfile
    from ..ensemble import EnsembleData
    from ..logging import set_filehandler

    # Load up the ensemble and configure settings
    cwd = os.getcwd()

    # run actions for which no complete setup is needed
    if args.version:
        printf(__version__)
        sys.exit(Returncode.OK)
    elif args.cleanup:
        cleanup_run(cwd)
        printf("Removed files and going to exit!")
        sys.exit(Returncode.OK)
    elif args.cleanup_all:
        cleanup_run(cwd, complete=True)
        printf("Removed files and going to exit!")
        sys.exit(Returncode.OK)
    elif args.writeconfig:
        write_rcfile(Path() / "censo2rc_NEW")
        sys.exit(Returncode.OK)

    if len(context.get("check", [])) == 0:
        printf("No tasks enabled! Stopping CENSO.")
        sys.exit(Returncode.OK)

    try:
        parts_config = configure(rcpath=args.inprcpath, args=args, context=context)
    except ValidationError as e:
        print_validation_errors(e)
        sys.exit(Returncode.CONFIG_ERROR)

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
        sys.exit(Returncode.INPUT_NOT_FOUND)

    if args.constraints:
        try:
            ensemble.constraints = Path(args.constraints).read_text()
        except FileNotFoundError:
            printf(
                f"Could not find constraints file {Path(args.constraints).resolve()}."
            )
            sys.exit(Returncode.CONSTRAINTS_NOT_FOUND)

    # if data should be reloaded, do it here
    if args.reload:
        for filename in args.reload:
            ensemble.read_output(Path(filename))

    # Set working directory to input parent
    # NOTE: this is the most convenient way to solve this for CLI version
    wd = Path(args.inp).parent.resolve()
    printf(f"Setting working directory to {wd}")
    os.chdir(wd)

    # Print all active parts settings once
    if logger.level == DEBUG:
        for _, config in parts_config:
            # Print
            printf(config)
    else:
        for fieldname, config in parts_config:
            if (
                getattr(args, fieldname, False)
                or fieldname == "general"
                or fieldname == "paths"
            ):
                # Print
                printf(config)

    printf("\n" + "".ljust(int(PLENGTH), "-") + "\n")

    # END of setup
    # -> ensemble.conformers contains all conformers with their info from input (sorted by CREST energy if possible)
    # -> settings are updated with cml args
    # -> output data is reloaded if wanted
    # -> working directory is set

    return ensemble, parts_config


def cleanup_run(cwd: str | Path, complete: bool = False):
    """
    Delete all unneeded files.

    :param cwd: Current working directory
    :param complete: If true, delete all files including ensembles and outputs
    """
    import shutil

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
    for subdir, _, files in os.walk(cwd):
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
    """
    Print a comparison table of the final rankings.

    :param comparison: Dictionary of comparison data
    """
    from tabulate import tabulate
    from ..utilities import h1

    if len(comparison) > 1:
        printf(h1("FINAL RANKING COMPARISON"))

        headers = ["CONF#"]

        headers.extend([f"ΔGtot {part.capitalize()}" for part in comparison])

        units = [
            "",
        ]

        units.extend(["[kcal/mol]" for _ in headers[1:]])

        confs = [conf for conf in list(comparison.values())[-1]]

        printmap = {"CONF#": lambda confname: confname}

        def make_printmap(h: str):
            return lambda confname: f"{comparison[h.split()[1].lower()][confname]:.2f}"

        for header in headers[1:]:
            printmap[header] = make_printmap(header)

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
