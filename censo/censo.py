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
- run ... --> part6()
"""
from os import getcwd
from time import perf_counter
import sys
from traceback import print_exc
from censo.cfg import PLENGTH, DESCR, __version__
from censo.inputhandling import cml
from censo.screening import part1
from censo.optimization import part2
from censo.refinement import part3
from censo.nmrproperties import part4
from censo.opticalrotation import part5
from censo.utilities import print
from censo.tutorial import interactiv_doc
from censo.core import CensoCore
from censo.prescreening import Prescreening

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
def main(argv=None):
    """
    Execute the CENSO code.
    """
    # parse command line into args
    args = cml(DESCR, argv)
    if args.version:
        print(__version__)
        sys.exit(0)

    if args.tutorial:
        interactiv_doc()
        sys.exit(0)

    # initialize blank core
    core = CensoCore.factory(getcwd(), args)

    # read input and setup conformers
    core.read_input()

    # read paths for external programs (definition in rcfile?)
    core.read_program_paths()

    ### END of setup
    # -> core.conformers contains all conformers with their info from input (sorted)
    # -> core.ensembledata is blank

    ### default: part1 - 3
    # TODO - reduce copy/paste code with list of functions which is iterated over
    run = [Prescreening, ]

    for part in run:
        print(f"Ran {part.name} in {part(core).run()} seconds!")
        
    # RUNNING PART0
    # cheap prescreening
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
