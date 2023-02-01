"""
CENSO run code:

"""
from os import getcwd
import os
import shutil
from time import perf_counter
import sys
from traceback import print_exc
from censo_test.cfg import PLENGTH, DESCR, ASSETS_PATH, __version__
from censo_test.inputhandling import cml
from censo_test.screening import part1
from censo_test.optimization import part2
from censo_test.refinement import part3
from censo_test.nmrproperties import part4
from censo_test.opticalrotation import part5
from censo_test.utilities import print
from censo_test.tutorial import interactiv_doc
from censo_test.core import CensoCore
from censo_test.prescreening import Prescreening
from censo_test.storage import CensoStorage
from censo_test.settings import InternalSettings

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
    
    # first, check program integrity
    # TODO - proper implementation?
    # FIXME - InternalSettings module is loaded before integrity is verified!!
    if not os.path.isdir(ASSETS_PATH):
        raise FileNotFoundError(ASSETS_PATH)
        
    # get setup info
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
        
    # setup storage with args (basically only rcpath and ensemblepath)
    storage = CensoStorage(args, cwd)
    
    # setup internal settings with args
    settings = InternalSettings(storage)
    
    # initialize core linked to storage
    core = CensoCore.factory(storage)

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