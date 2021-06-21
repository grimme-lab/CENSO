"""
CENSO run code:
- reading commandline input --> cml()
- parsing remote configuration file, reading conformer ensemble,
  checking parameters and creating or reading enso.json conformer information
  --> enso_startup()
- run cheap_screeing --> part0()
- run prescreening --> part1()
- run optimization --> part2()
- run refinement   --> part3()
- run nmrproperties --> part4() or
- run optical_rotation --part5()
"""
from os import getcwd
from time import perf_counter
import sys
from traceback import print_exc
from .cfg import PLENGTH, DESCR, __version__
from .inputhandling import cml, internal_settings
from .setupcenso import enso_startup
from .cheapscreening import part0
from .prescreening import part1
from .optimization import part2
from .refinement import part3
from .nmrproperties import part4
from .opticalrotation import part5
from .utilities import print
from .tutorial import interactiv_doc


def main(argv=None):
    """
    Execute the CENSO code.
    """
    # get commandline input:
    args = cml(DESCR, internal_settings(), argv)
    if args.version:
        print(__version__)
        sys.exit(0)

    if args.tutorial:
        interactiv_doc()
        sys.exit(0)

    # setup conformers and process input: cml >> configfile > internal defaults
    args, config, conformers, ensembledata = enso_startup(getcwd(), args)

    # RUNNING PART0
    if config.part0:
        tic = perf_counter()
        try:
            config, conformers, store_confs, ensembledata = part0(
                config, conformers, ensembledata
            )
        except Exception as error:
            print("ERROR in part0!")
            print("\nThe error-message is {}\n".format(error))
            print("Traceback for debugging:".center(PLENGTH, "*"))
            print_exc()
            print("".center(PLENGTH, "*"))
            print("Going to exit!")
            sys.exit(1)
        toc = perf_counter()
        ensembledata.part_info["part0"] = toc - tic
        ensembledata.previous_part_info["part0"] += ensembledata.part_info["part0"]
        print(f"Ran part0 in {ensembledata.part_info['part0']:0.4f} seconds")

    # RUNNING PART1
    if config.part1:
        tic = perf_counter()
        try:
            store_confs
        except NameError:
            store_confs = []
        try:
            config, conformers, store_confs, ensembledata = part1(
                config, conformers, store_confs, ensembledata
            )
        except Exception as error:
            print("ERROR in part1!")
            print("\nThe error-message is {}\n".format(error))
            print("Traceback for debugging:".center(PLENGTH, "*"))
            print_exc()
            print("".center(PLENGTH, "*"))
            print("Going to exit!")
            sys.exit(1)
        toc = perf_counter()
        ensembledata.part_info["part1"] = toc - tic
        ensembledata.previous_part_info["part1"] += ensembledata.part_info["part1"]
        ensembledata.previous_part_info["part1_firstsort"] += ensembledata.part_info[
            "part1_firstsort"
        ]
        print(f"Ran part1 in {ensembledata.part_info['part1']:0.4f} seconds")

    # RUNNING PART2
    if config.part2:
        tic = perf_counter()
        try:
            store_confs
        except NameError:
            store_confs = []
        try:
            config, conformers, store_confs, ensembledata = part2(
                config, conformers, store_confs, ensembledata
            )
        except Exception as error:
            print("ERROR in part2!")
            print("\nThe error-message is {}\n".format(error))
            print("Traceback for debugging:".center(PLENGTH, "*"))
            print_exc()
            print("".center(PLENGTH, "*"))
            print("Going to exit!")
            sys.exit(1)
        toc = perf_counter()
        ensembledata.part_info["part2"] = toc - tic
        ensembledata.previous_part_info["part2"] += ensembledata.part_info["part2"]
        ensembledata.previous_part_info["part2_opt"] += ensembledata.part_info[
            "part2_opt"
        ]
        print(f"Ran part2 in {ensembledata.part_info['part2']:0.4f} seconds")

    # RUNNING PART3
    if config.part3:
        tic = perf_counter()
        try:
            store_confs
        except NameError:
            store_confs = []
        try:
            config, conformers, store_confs, ensembledata = part3(
                config, conformers, store_confs, ensembledata
            )
        except Exception as error:
            print("ERROR in part3!")
            print("\nThe error-message is {}\n".format(error))
            print("Traceback for debugging:".center(PLENGTH, "*"))
            print_exc()
            print("".center(PLENGTH, "*"))
            print("Going to exit!")
            sys.exit(1)
        toc = perf_counter()
        ensembledata.part_info["part3"] = toc - tic
        ensembledata.previous_part_info["part3"] += ensembledata.part_info["part3"]
        print(f"Ran part3 in {ensembledata.part_info['part3']:0.4f} seconds")

    # RUNNING PART4
    if config.part4:
        tic = perf_counter()
        try:
            store_confs
        except NameError:
            store_confs = []
        try:
            config, conformers, store_confs, ensembledata = part4(
                config, conformers, store_confs, ensembledata
            )
        except Exception as error:
            print("ERROR in part4!")
            print("\nThe error-message is {}\n".format(error))
            print("Traceback for debugging:".center(PLENGTH, "*"))
            print_exc()
            print("".center(PLENGTH, "*"))
            print("Going to exit!")
            sys.exit(1)
        toc = perf_counter()
        ensembledata.part_info["part4"] = toc - tic
        ensembledata.previous_part_info["part4"] += ensembledata.part_info["part4"]
        print(f"Ran part4 in {ensembledata.part_info['part4']:0.4f} seconds")

    # RUNNING PART5
    if config.optical_rotation:
        tic = perf_counter()
        try:
            store_confs
        except NameError:
            store_confs = []
        try:
            config, conformers, store_confs, ensembledata = part5(
                config, conformers, store_confs, ensembledata
            )
        except Exception as error:
            print("ERROR in part5!")
            print("\nThe error-message is {}\n".format(error))
            print("Traceback for debugging:".center(PLENGTH, "*"))
            print_exc()
            print("".center(PLENGTH, "*"))
            print("Going to exit!")
            sys.exit(1)
        toc = perf_counter()
        ensembledata.part_info["part5"] = toc - tic
        ensembledata.previous_part_info["part5"] += ensembledata.part_info["part5"]
        print(f"Ran part5 in {ensembledata.part_info['part5']:0.4f} seconds")

    # save current data to jsonfile
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in conformers]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    if args.create_si:
        config.create_SI(ensembledata)

    # END of CENSO
    timings = 0.0
    prev_timings = 0.0
    if len(str(config.nconf)) > 5:
        conflength = len(str(config.nconf))
    else:
        conflength = 5

    try:
        tmp = []
        tmp_prev = []
        if config.part0:
            tmp.append(ensembledata.part_info["part0"])
            tmp_prev.append(ensembledata.previous_part_info["part0"])
        if config.part1:
            tmp.append(ensembledata.part_info["part1"])
            tmp_prev.append(ensembledata.previous_part_info["part1"])
            tmp_prev.append(ensembledata.previous_part_info["part1_firstsort"])
        if config.part2:
            tmp.append(ensembledata.part_info["part2"])
            tmp_prev.append(ensembledata.previous_part_info["part2"])
            tmp.append(ensembledata.part_info["part2_opt"])
            tmp_prev.append(ensembledata.previous_part_info["part2_opt"])
        if config.part3:
            tmp.append(ensembledata.part_info["part3"])
            tmp_prev.append(ensembledata.previous_part_info["part3"])
        if config.part4:
            tmp.append(ensembledata.part_info["part4"])
            tmp_prev.append(ensembledata.previous_part_info["part4"])
        if config.optical_rotation:
            tmp.append(ensembledata.part_info["part5"])
            tmp_prev.append(ensembledata.previous_part_info["part5"])
        timelen = max([len(f"{float(value):.2f}") for value in tmp]) + 2
        prev_timelen = max([len(f"{float(value):.2f}") for value in tmp_prev]) + 2
        if timelen < 7:
            timelen = 7
    except Exception:
        timelen = 20
        prev_timelen = 20

    print(
        f"\n\n{'Part':20}: {'#conf':>{conflength}}    {'time': >{timelen}}      time (including restarts)"
    )
    print("".ljust(int(PLENGTH / 1.4), "-"))
    print(
        f"{'Input':20}: {ensembledata.nconfs_per_part['starting']:{conflength}}    {'-':^{timelen+2}}    {'-':^{timelen+2}}"
    )
    if config.part0:
        print(
            f"{'Part0_all':20}: {ensembledata.nconfs_per_part['part0']:{conflength}}    {ensembledata.part_info['part0']:{timelen}.2f} s  {ensembledata.previous_part_info['part0']:>{prev_timelen}.2f} s"
        )
        timings += ensembledata.part_info["part0"]
        prev_timings += ensembledata.previous_part_info["part0"]
    if config.part1:
        print(
            f"{'Part1_initial_sort':20}: {ensembledata.nconfs_per_part['part1_firstsort']:{conflength}}    {ensembledata.part_info['part1_firstsort']:{timelen}.2f} s  {ensembledata.previous_part_info['part1_firstsort']:>{prev_timelen}.2f} s"
        )
        print(
            f"{'Part1_all':20}: {ensembledata.nconfs_per_part['part1_firstsort']:{conflength}}    {ensembledata.part_info['part1']:{timelen}.2f} s  {ensembledata.previous_part_info['part1']:>{prev_timelen}.2f} s"
        )
        timings += ensembledata.part_info["part1"]
        prev_timings += ensembledata.previous_part_info["part1"]
    if config.part2:
        print(
            f"{'Part2_opt':20}: {ensembledata.nconfs_per_part['part2_opt']:{conflength}}    {ensembledata.part_info['part2_opt']:{timelen}.2f} s  {ensembledata.previous_part_info['part2_opt']:>{prev_timelen}.2f} s"
        )
        print(
            f"{'Part2_all':20}: {ensembledata.nconfs_per_part['part2']:{conflength}}    {ensembledata.part_info['part2']:{timelen}.2f} s  {ensembledata.previous_part_info['part2']:>{prev_timelen}.2f} s"
        )
        timings += ensembledata.part_info["part2"]
        prev_timings += ensembledata.previous_part_info["part2"]
    if config.part3:
        print(
            f"{'Part3_all':20}: {ensembledata.nconfs_per_part['part3']:{conflength}}    {ensembledata.part_info['part3']:{timelen}.2f} s  {ensembledata.previous_part_info['part3']:>{prev_timelen}.2f} s"
        )
        timings += ensembledata.part_info["part3"]
        prev_timings += ensembledata.previous_part_info["part3"]
    if config.part4:
        print(
            f"{'Part4':20}: {ensembledata.nconfs_per_part['part4']:{conflength}}    {ensembledata.part_info['part4']:{timelen}.2f} s  {ensembledata.previous_part_info['part4']:>{prev_timelen}.2f} s"
        )
        timings += ensembledata.part_info["part4"]
        prev_timings += ensembledata.previous_part_info["part4"]
    if config.optical_rotation:
        print(
            f"{'Part5':20}: {ensembledata.nconfs_per_part['part5']:{conflength}}    {ensembledata.part_info['part5']:{timelen}.2f} s  {ensembledata.previous_part_info['part5']:>{prev_timelen}.2f} s"
        )
        timings += ensembledata.part_info["part5"]
        prev_timings += ensembledata.previous_part_info["part5"]
    print("".ljust(int(PLENGTH / 1.4), "-"))
    print(
        f"{'All parts':20}: {'-':>{conflength}}    {timings:{timelen}.2f} s  {prev_timings:{prev_timelen}.2f} s"
    )
    print("\nCENSO all done!")
    return 0
