"""
Contains enso_startup for the initialization of all parameters set for the
suseqent calculation.
"""
import os
import sys
import json
from collections import OrderedDict
from .cfg import CODING, PLENGTH, DESCR, WARNLEN, censo_solvent_db, __version__, NmrRef
from .inputhandling import config_setup, internal_settings
from .datastructure import MoleculeData
from .qm_job import QmJob
from .utilities import (
    mkdir_p,
    do_md5,
    t2x,
    move_recursively,
    get_energy_from_ensemble,
    frange,
    print,
)
from .ensembledata import EnsembleData


def enso_startup(cwd, args):
    """
    1) read cml input,
    2) print header
    3) read or create enso control file '.censorc'
    4) read or write flags.dat control file
    5) check for crest_conformers.xyz
    6) check program settings
    7) read or write enso.json
    """

    print(DESCR)
    config = config_setup(path=os.path.abspath(cwd))

    if args.cleanup:
        print("Cleaning up the directory from unneeded files!")
        config.cleanup_run()
        print("Removed files and going to exit!")
        sys.exit(0)
    elif args.cleanup_all:
        print("Cleaning up the directory from ALL unneeded files!")
        config.cleanup_run(True)
        print("Removed files and going to exit!")
        sys.exit(0)

    configfname = ".censorc"
    if args.writeconfig:
        # check if .censorc in local or home dir and ask user it the program
        # paths should be copied
        tmp = None
        newconfigfname = "censorc_new"
        if os.path.isfile(os.path.join(config.cwd, configfname)):
            tmp = os.path.join(config.cwd, configfname)
        elif os.path.isfile(os.path.join(os.path.expanduser("~"), configfname)):
            tmp = os.path.join(os.path.expanduser("~"), configfname)
        if tmp is not None:
            print(f"An existing .censorc has been found in {tmp}")
            print(
                f"Do you want to copy existing program path information to the "
                f"new remote configuration file?"
            )
            print("Please type 'yes' or 'no':")
            user_input = input()
            if user_input.strip() in ("y", "yes"):
                config.read_program_paths(tmp)
                config.write_rcfile(
                    os.path.join(config.cwd, newconfigfname), usepaths=True
                )
                print("")
            elif user_input.strip() in ("n", "no"):
                config.write_rcfile(os.path.join(config.cwd, newconfigfname))
        else:
            config.write_rcfile(os.path.join(config.cwd, newconfigfname))
        print(
            "A new ensorc was written into the current directory file: "
            f"{newconfigfname}!\nYou have to adjust the settings to your needs"
            " and it is mandatory to correctly set the program paths!\n"
            "Additionally move the file to the correct filename: '.censorc'\n"
            "and place it either in your /home/$USER/ or current directory.\n"
            "\nAll done!"
        )
        sys.exit(0)

    if args.inprcpath:
        try:
            tmp_path = os.path.abspath(os.path.expanduser(args.inprcpath))
            if os.path.isfile(tmp_path):
                config.configpath = tmp_path
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            print(
                f"{'ERROR:':{WARNLEN}}Could not find the configuration file: {configfname}."
            )
            print("\nGoing to exit!")
            sys.exit(1)
    elif (
        args.restart
        and os.path.isfile(os.path.join(config.cwd, "enso.json"))
        and os.path.isfile(
            config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
            .get("settings", {})
            .get("configpath", "")
        )
    ):
        # read configpath from previous enso.json if restart is requested
        config.configpath = (
            config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
            .get("settings", {})
            .get("configpath", "")
        )
    elif os.path.isfile(os.path.join(config.cwd, configfname)):
        # local configuration file before remote configuration file
        config.configpath = os.path.join(config.cwd, configfname)
    elif os.path.isfile(os.path.join(os.path.expanduser("~"), configfname)):
        # remote configuration file
        config.configpath = os.path.join(os.path.expanduser("~"), configfname)
    else:
        print(
            f"{'ERROR:':{WARNLEN}}Could not find the configuration file: {configfname}.\n"
            f"{'':{WARNLEN}}The file has to be either in /home/$USER/ or the current "
            "working directory!\n"
            "The configurationfile can be otherwise directly referenced using: "
            "'censo -inprc /path/to/.censorc'"
        )
        print("\nGoing to exit!")
        sys.exit(1)

    ### solvent database adjustable by user
    censo_assets_path = os.path.expanduser("~/.censo_assets")
    if not os.path.isdir(censo_assets_path):
        mkdir_p(censo_assets_path)
    solvent_user_path = os.path.expanduser(
        os.path.join("~/.censo_assets/", "censo_solvents.json")
    )
    if os.path.isfile(solvent_user_path):
        config.save_infos.append(
            "Reading file: {}\n".format(os.path.basename(solvent_user_path))
        )
        try:
            with open(solvent_user_path, "r", encoding=CODING, newline=None) as inp:
                censo_solvent_db.update(json.load(inp, object_pairs_hook=OrderedDict))
        except (ValueError, TypeError, FileNotFoundError):
            print(
                f"{'ERROR:':{WARNLEN}}Your censo_solvents.json file in {solvent_user_path} is corrupted!\n"
                f"{'':{WARNLEN}}You can delete your corrupted file and a new censo_solvents.json will be created on the next start of CENSO."
            )
            raise

    else:
        with open(solvent_user_path, "w") as out:
            json.dump(censo_solvent_db, out, indent=4, sort_keys=True)
        config.save_infos.append(
            "Creating file: {}".format(os.path.basename(solvent_user_path))
        )
    ### END solvent database adjustable by user
    ### NMR reference shielding constant database adjustable by user
    if not os.path.isdir(censo_assets_path):
        mkdir_p(censo_assets_path)
    nmr_ref_user_path = os.path.expanduser(
        os.path.join("~/.censo_assets/", "censo_nmr_ref.json")
    )
    if not os.path.isfile(nmr_ref_user_path):
        with open(nmr_ref_user_path, "w") as out:
            tmp = NmrRef()
            json.dump(tmp, out, default=NmrRef.NMRRef_to_dict, indent=4, sort_keys=True)
        config.save_infos.append(
            "Creating file: {}".format(os.path.basename(nmr_ref_user_path))
        )
    ### END NMR reference shielding constant database adjustable by user

    ### if restart read all settings from previous run (enso.json)
    if args.restart and os.path.isfile(os.path.join(config.cwd, "enso.json")):
        tmp = config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
        previous_settings = tmp.get("settings")
        for key, value in previous_settings.items():
            if vars(args).get(key, "unKn_own") == "unKn_own":
                # print(key, 'not_known')
                continue
            if getattr(args, key, "unKn_own") is None:
                setattr(args, key, value)
    ### END if restart

    if config.configpath:
        # combine args und comandline
        # check if startread in file:
        startread = "$CRE SORTING SETTINGS:"
        with open(config.configpath, "r") as myfile:
            try:
                data = myfile.readlines()
                censorc_version = "0.0.0"
                for line in data:
                    if "$VERSION" in line:
                        censorc_version = line.split(":")[1]
                if int(censorc_version.split(".")[1]) < int(
                    __version__.split(".")[1]
                ) or int(censorc_version.split(".")[0]) < int(
                    __version__.split(".")[0]
                ):
                    print(
                        f"{'ERROR:':{WARNLEN}}There has been an API break and you have to "
                        f"create a new .censorc.\n{'':{WARNLEN}}E.g. 'censo -newconfig'"
                    )
                    print(
                        f"{'INFORMATION:':{WARNLEN}}Due to the API break data in ~/.censo_assets/ might need updating.\n"
                        f"{'':{WARNLEN}}Therefore delete the files {'censo_nmr_ref.json'} and {'censo_solvents.json'} so "
                        "that they can be re-created upon the next censo call."
                    )
                    sys.exit(1)
                myfile.seek(0)  # reset reader
            except (ValueError, KeyError, AttributeError) as e:
                print(e)
                print(
                    f"{'ERROR:':{WARNLEN}}Please create a new .censorc --> 'censo -newconfig'"
                )
                sys.exit(1)
            if not startread in myfile.read():
                print(
                    f"{'ERROR:':{WARNLEN}}You are using a corrupted .censorc. Create a new one!"
                )
                sys.exit(1)
        # combine commandline and .censorc
        config.read_config(config.configpath, startread, args)

    if args.copyinput:
        config.read_program_paths(config.configpath)
        config.write_censo_inp(config.cwd)
        print(
            "The file censo.inp with the current settings has been written to "
            "the current working directory."
        )
        print("\nGoing to exit!")
        sys.exit(1)

    # read ensemble input file (e.g. crest_conformers.xyz)
    if os.path.isfile(os.path.join(config.cwd, "enso.json")):
        tmp = config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
        if "ensemble_info" in tmp and args.inp is None:
            inpfile = os.path.basename(tmp["ensemble_info"].get("filename"))
            if os.path.isfile(inpfile):
                args.inp = inpfile
            if args.debug:
                print(f"Using Input file from: {inpfile}")
    if args.inp is None:
        args.inp = "crest_conformers.xyz"
    if os.path.isfile(args.inp):
        config.ensemblepath = args.inp
        # identify coord or xyz trajectory
        config.md5 = do_md5(config.ensemblepath)
        with open(config.ensemblepath, "r", encoding=CODING, newline=None) as inp:
            foundcoord = False
            for line in inp:
                if "$coord" in line:
                    foundcoord = True
                    break
        if foundcoord:
            _, config.nat = t2x(
                config.ensemblepath, writexyz=True, outfile="converted.xyz"
            )
            config.ensemblepath = os.path.join(config.cwd, "converted.xyz")
            config.maxconf = 1
            config.nconf = 1
        else:
            with open(
                config.ensemblepath, "r", encoding=CODING, newline=None
            ) as infile:
                try:
                    config.nat = int(infile.readline().strip().split()[0])
                    filelen = 1
                except (ValueError, TypeError, IndexError) as e:
                    print(
                        f"{'ERROR:':{WARNLEN}}Could not get the number of atoms or the "
                        "number of conformers from the inputfile "
                        f"{os.path.basename(args.inp)}"
                    )
                    sys.exit(1)
                for line in infile:
                    filelen += 1
                try:
                    config.maxconf = int(filelen / (config.nat + 2))
                    if filelen % (config.nat + 2) != 0:
                        raise ValueError
                except ValueError:
                    print(
                        f"{'ERROR:':{WARNLEN}}Could not get the number of atoms or the "
                        "number of conformers from the inputfile "
                        f"{os.path.basename(args.inp)}"
                    )
                    sys.exit(1)
    else:
        print(f"{'ERROR:':{WARNLEN}}The input file can not be found!")
        sys.exit(1)

    # determine number of conformers:
    if args.nconf:
        if args.nconf > config.maxconf:
            config.nconf = config.maxconf
        else:
            config.nconf = args.nconf
    elif isinstance(config.nconf, int):
        pass
        # keep from .censorc
    else:
        config.nconf = config.maxconf

    # check settings-combination and show error:
    error_logical = config.check_logic()

    if getattr(args, "onlyread", None) is not None:
        if getattr(args, "onlyread") == "on":
            setattr(config, "onlyread", True)
        elif getattr(args, "onlyread") == "off":
            setattr(config, "onlyread", False)

    # printing parameters
    config.print_parameters(cmlcall=sys.argv)
    config.read_program_paths(config.configpath)
    requirements = config.needed_external_programs()
    error_logical = config.processQMpaths(requirements, error_logical)
    # print errors
    if config.save_errors:
        print("")
        print("".ljust(PLENGTH, "*"))
        for _ in list(config.save_errors):
            print(config.save_errors.pop(0))
        print("".ljust(PLENGTH, "*"))
    if config.onlyread:
        config.save_errors.append(
            f"{'WARNING:':{WARNLEN}}Using the option ``readonly`` to re-read"
            + f" data from old outputs. This option is experimental therefore "
            + f"check results carefully.\n"
            + f"{'':{WARNLEN}}It can only succeed if exactly the same "
            + f"input commands (.censorc + command line) are supplied!"
        )
    # END print error

    if error_logical and not args.debug and args.checkinput:
        print(
            f"\n{'ERROR:':{WARNLEN}}CENSO can not continue due to input errors!\n"
            f"{'':{WARNLEN}}Fix errors and run censo -checkinput again!"
        )
        print("\nGoing to exit!")
        sys.exit(1)
    elif error_logical and not args.debug:
        print(f"\n{'ERROR:':{WARNLEN}}CENSO can not continue due to input errors!")
        print("\nGoing to exit!")
        sys.exit(1)

    if not error_logical or args.debug:
        print("\n" + "".ljust(PLENGTH, "-"))
        print(" Processing data from previous run (enso.json)".center(PLENGTH, " "))
        print("".ljust(PLENGTH, "-") + "\n")
        # read enso.json
        if os.path.isfile(os.path.join(config.cwd, "enso.json")):
            config.jsonpath = os.path.join(config.cwd, "enso.json")
            save_data = config.read_json(config.jsonpath)
            # Check if settings and "ensemble_info" are present else end!
            if "settings" not in save_data or "ensemble_info" not in save_data:
                print(
                    f"{'ERROR:':{WARNLEN}}important information for restarting missing from "
                    f"{config.jsonpath}!"
                )
                print("\nGoing to exit!")
                sys.exit(1)
            # create previous run data
            previousrun = config_setup(internal_settings)
            for item in save_data["settings"].keys():
                setattr(previousrun, item, save_data["settings"].get(item))
            # end previous run data
            if getattr(previousrun, "fixed_temperature", "unKn_own") == "unKn_own":
                print(
                    f"{'ERROR:':{WARNLEN}}important information ('fixed_temperature') for restarting missing from "
                    f"{config.jsonpath}!"
                )
                print("\nGoing to exit!")
                sys.exit(1)
            else:
                if isinstance(
                    getattr(previousrun, "fixed_temperature", "unKn_own"), float
                ):
                    setattr(
                        config,
                        "fixed_temperature",
                        (getattr(previousrun, "fixed_temperature")),
                    )
                else:
                    print(
                        f"{'ERROR:':{WARNLEN}}Can't read fixed_temperature information necessary for restart!"
                    )
                    print("\nGoing to exit!")
                    sys.exit(1)
            ### test concerning temperature changes for providing Boltzmann weights
            ### at different temperatures
            prev_temp = getattr(previousrun, "temperature", None)
            if prev_temp != getattr(config, "temperature"):
                print(
                    f"{'INFORMATION:':{WARNLEN}}The temperature has been changed from (previous) {prev_temp} K to (current) {getattr(config, 'temperature')} K."
                )
                # check multitemp on in previous
                if getattr(previousrun, "multitemp", None):
                    # has trange been changed ?
                    prev_trange = getattr(previousrun, "trange", [])
                    if not isinstance(prev_trange, list) or prev_trange != getattr(
                        config, "trange"
                    ):
                        print(
                            f"{'ERROR:':{WARNLEN}}The temperature range has been changed together with the temperature which is not allowed."
                        )
                        print("\nGoing to exit!")
                        sys.exit(1)
                    # is the temperature in trange
                    if float(getattr(config, "temperature")) in [
                        i
                        for i in frange(prev_trange[0], prev_trange[1], prev_trange[2])
                    ]:
                        # have to set everything to correct new temperature!
                        print(
                            f"{'INFORMATION:':{WARNLEN}}Temperature change is possible!"
                        )
                    else:
                        print(
                            f"{'ERROR:':{WARNLEN}}Temperature has not been previously calculated, e.g. not in temperature range (trange). Change of temperature not possible!"
                        )
                        print(
                            f"{'INFORMATION:':{WARNLEN}}Temperatures in trange are: {[i for i in frange(prev_trange[0], prev_trange[1], prev_trange[2])]}"
                        )
                        print("\nGoing to exit!")
                        sys.exit(1)
                else:
                    print(
                        f"{'ERROR:':{WARNLEN}}Temperature change is requested, but no previous other temperatures have been calculated!"
                    )
                    print("\nGoing to exit!")
                    sys.exit(1)
            ### END test temperatures

            if config.md5 != previousrun.md5:
                print(
                    f"{'WARNING:':{WARNLEN}}The inputfile containing all conformers was "
                    "changed, compared to the previous run!"
                )
            for flag in config.restart_unchangeable:
                if getattr(config, flag, "None") != getattr(previousrun, flag, "None2"):
                    print(
                        f"{'ERROR:':{WARNLEN}}setting {flag} was changed from "
                        f"{getattr(config, flag, 'None')} to {getattr(previousrun, flag, 'None')}!"
                    )
                    error_logical = True
            if (
                getattr(config, "evaluate_rrho", "None")
                != getattr(previousrun, "evaluate_rrho", "None2")
            ) and getattr(config, "part2", "None"):
                print(
                    f"{'ERROR:':{WARNLEN}}setting {'evaluate_rrho'} can not be changed "
                    f"in geometry optimization!\n"
                )
                error_logical = True
            if error_logical and not args.debug:
                print(
                    f"{'ERROR:':{WARNLEN}}All flags which are concerned with geometry "
                    f"optimization \n{'':{WARNLEN}}(func, prog, ancopt, sm, solv, chrg, "
                    "unpaired) are not allowed to be changed!\n"
                    f"{'':{WARNLEN}}If you want to change these settings, "
                    "start from scratch in a new folder!"
                )
                print("\nGoing to exit!")
                sys.exit(1)
            if not args.checkinput:
                move_recursively(config.cwd, os.path.basename(config.jsonpath))

            # Check if flags have been changed between two runs and adjust data
            # e.g. reset or load previously calculated
            for flag in config.restart_changeable.keys():
                if getattr(config, flag, "None") != getattr(previousrun, flag, "None2"):
                    print(
                        f"{'WARNING:':{WARNLEN}}setting {flag} was changed from "
                        f"{getattr(previousrun, flag, 'None')} to "
                        f"{getattr(config, flag, 'None')}!"
                    )
                    config.restart_changeable[flag] = True
                    if (
                        flag == "multitemp"
                        and getattr(config, flag, "None")
                        and not getattr(previousrun, flag, "None")
                    ):
                        # multitemp only reset if not calculated before!
                        #  --> off --> on
                        print(
                            f"{'WARNING:':{WARNLEN}}{flag} is requested and the different "
                            "temperatures have not been evaluated in the\n"
                            f"{'':WARNLEN}previous run! Resetting calculations concerning trange!"
                        )
                    elif (
                        flag == "multitemp"
                        and not getattr(config, flag, "None")
                        and getattr(previousrun, flag, "None")
                    ):
                        # multitemp only reset if not calculated before!
                        #  --> off --> on
                        config.restart_changeable[flag] = False
                    if flag == "trange":
                        # if temp in trange has not been calculated reset!
                        prev_t = getattr(previousrun, flag)
                        prev_trange = [
                            i for i in frange(prev_t[0], prev_t[1], prev_t[2])
                        ]
                        cur_t = getattr(config, flag)
                        current_trange = [
                            i for i in frange(cur_t[0], cur_t[1], cur_t[2])
                        ]
                        for temp in current_trange:
                            if temp not in prev_trange:
                                config.restart_changeable[flag] = True
                                break
                            else:
                                config.restart_changeable[flag] = False

            conformers = []
            for conf in save_data.keys():
                if conf == "ensemble_info":
                    ensembledata = EnsembleData(
                        id=save_data[conf].get("id"),
                        filename=save_data[conf].get("filename"),
                        part_info=save_data[conf].get("part_info"),
                        previous_part_info=save_data[conf].get(
                            "previous_part_info", EnsembleData().previous_part_info
                        ),
                        avGcorrection=save_data[conf].get("avGcorrection"),
                        comment=save_data[conf].get("comment"),
                        bestconf=save_data[conf].get("bestconf"),
                        nconfs_per_part=save_data[conf].get("nconfs_per_part"),
                    )
                    ensembledata.nconfs_per_part["starting"] = config.nconf
                elif conf not in ("settings", "ensemble_info"):
                    for info in vars(MoleculeData(0)).keys():
                        if save_data[conf].get(info, "xXx") == "xXx":
                            print(
                                f"{'WARNING:':{WARNLEN}}Missing data {info} from enso.json! "
                                "Default is added."
                            )
                    molecule = QmJob(
                        save_data[conf].get("id"),
                        chrg=save_data[conf].get("chrg"),
                        uhf=save_data[conf].get("uhf"),
                        xtb_energy=save_data[conf].get("xtb_energy"),
                        xtb_energy_unbiased=save_data[conf].get("xtb_energy_unbiased"),
                        xtb_free_energy=save_data[conf].get("xtb_free_energy"),
                        rel_xtb_energy=save_data[conf].get("rel_xtb_energy"),
                        rel_xtb_free_energy=save_data[conf].get("rel_xtb_free_energy"),
                        sym=save_data[conf].get("sym", getattr(MoleculeData(0), "sym")),
                        linear=save_data[conf].get(
                            "linear", getattr(MoleculeData(0), "linear")
                        ),
                        symnum=save_data[conf].get(
                            "symnum",
                            MoleculeData(0)._get_sym_num(
                                sym=save_data[conf].get(
                                    "sym", getattr(MoleculeData(0), "sym")
                                ),
                                linear=save_data[conf].get(
                                    "linear", getattr(MoleculeData(0), "linear")
                                ),
                            ),
                        ),
                        gi=save_data[conf].get("gi", getattr(MoleculeData(0), "sym")),
                        removed=save_data[conf].get(
                            "removed", getattr(MoleculeData(0), "removed")
                        ),
                        temperature_info=save_data[conf].get(
                            "temperature_info",
                            getattr(MoleculeData(0), "temperature_info"),
                        ),
                        cheap_prescreening_sp_info=save_data[conf].get(
                            "cheap_prescreening_sp_info",
                            getattr(MoleculeData(0), "cheap_prescreening_sp_info"),
                        ),
                        cheap_prescreening_gsolv_info=save_data[conf].get(
                            "cheap_prescreening_gsolv_info",
                            getattr(MoleculeData(0), "cheap_prescreening_gsolv_info"),
                        ),
                        prescreening_sp_info=save_data[conf].get(
                            "prescreening_sp_info",
                            getattr(MoleculeData(0), "prescreening_sp_info"),
                        ),
                        lowlevel_sp_info=save_data[conf].get(
                            "lowlevel_sp_info",
                            getattr(MoleculeData(0), "lowlevel_sp_info"),
                        ),
                        highlevel_sp_info=save_data[conf].get(
                            "highlevel_sp_info",
                            getattr(MoleculeData(0), "highlevel_sp_info"),
                        ),
                        prescreening_grrho_info=save_data[conf].get(
                            "prescreening_grrho_info",
                            getattr(MoleculeData(0), "prescreening_grrho_info"),
                        ),
                        lowlevel_grrho_info=save_data[conf].get(
                            "lowlevel_grrho_info",
                            getattr(MoleculeData(0), "lowlevel_grrho_info"),
                        ),
                        lowlevel_hrrho_info=save_data[conf].get(
                            "lowlevel_hrrho_info",
                            getattr(MoleculeData(0), "lowlevel_hrrho_info"),
                        ),
                        highlevel_grrho_info=save_data[conf].get(
                            "highlevel_grrho_info",
                            getattr(MoleculeData(0), "highlevel_grrho_info"),
                        ),
                        highlevel_hrrho_info=save_data[conf].get(
                            "highlevel_hrrho_info",
                            getattr(MoleculeData(0), "highlevel_hrrho_info"),
                        ),
                        prescreening_gsolv_info=save_data[conf].get(
                            "prescreening_gsolv_info",
                            getattr(MoleculeData(0), "prescreening_gsolv_info"),
                        ),
                        lowlevel_gsolv_info=save_data[conf].get(
                            "lowlevel_gsolv_info",
                            getattr(MoleculeData(0), "lowlevel_gsolv_info"),
                        ),
                        lowlevel_gsolv_compare_info=save_data[conf].get(
                            "lowlevel_gsolv_compare_info",
                            getattr(MoleculeData(0), "lowlevel_gsolv_compare_info"),
                        ),
                        highlevel_gsolv_info=save_data[conf].get(
                            "highlevel_gsolv_info",
                            getattr(MoleculeData(0), "highlevel_gsolv_info"),
                        ),
                        optimization_info=save_data[conf].get(
                            "optimization_info",
                            getattr(MoleculeData(0), "optimization_info"),
                        ),
                        nmr_coupling_info=save_data[conf].get(
                            "nmr_coupling_info",
                            getattr(MoleculeData(0), "nmr_coupling_info"),
                        ),
                        nmr_shielding_info=save_data[conf].get(
                            "nmr_shielding_info",
                            getattr(MoleculeData(0), "nmr_shielding_info"),
                        ),
                        part_info=save_data[conf].get(
                            "part_info", getattr(MoleculeData(0), "part_info")
                        ),
                        comment=save_data[conf].get(
                            "comment", getattr(MoleculeData(0), "comment")
                        ),
                        optical_rotation_info=save_data[conf].get(
                            "optical_rotation_info",
                            getattr(MoleculeData(0), "optical_rotation_info"),
                        ),
                    )

                    # adjust to restart changeable data:
                    for key, value in config.restart_changeable.items():
                        if value and key == "multitemp":
                            molecule.reset_range_info(
                                trange=[
                                    i
                                    for i in frange(
                                        config.trange[0],
                                        config.trange[1],
                                        config.trange[2],
                                    )
                                ]
                            )
                        elif value and key == "trange":
                            molecule.reset_range_info(
                                trange=[
                                    i
                                    for i in frange(
                                        config.trange[0],
                                        config.trange[1],
                                        config.trange[2],
                                    )
                                ]
                            )
                        elif value and key in (
                            "part1_gfnv",
                            "part2_gfnv",
                            "part3_gfnv",
                        ):
                            exc = {
                                "part1_gfnv": "prescreening_grrho_info",
                                "part2_gfnv": "lowlevel_grrho_info",
                                "part3_gfnv": "highlevel_grrho_info",
                            }
                            if getattr(config, key) != getattr(previousrun, key):
                                # save calculated to
                                molecule.save_prev(
                                    exc[key], getattr(molecule, exc[key]).get("method")
                                )
                                # load new if available
                                method, _ = config.get_method_name(
                                    "rrhoxtb",
                                    gfn_version=getattr(config, key),
                                    bhess=config.bhess,
                                )
                                molecule.load_prev(exc[key], method)
                        elif value and key in ("smgsolv1", "smgsolv2", "smgsolv3"):
                            exc = {
                                "smgsolv1": "prescreening_gsolv_info",
                                "smgsolv2": "lowlevel_gsolv_info",
                                "smgsolv3": "highlevel_gsolv_info",
                            }
                            exc_implicit = {
                                "smgsolv1": "prescreening_sp_info",
                                "smgsolv2": "lowlevel_sp_info",
                                "smgsolv3": "highlevel_sp_info",
                            }
                            exc2 = {
                                "smgsolv1": "part1_gfnv",
                                "smgsolv2": "part2_gfnv",
                                "smgsolv3": "part3_gfnv",
                            }
                            if getattr(config, key) != getattr(previousrun, key):
                                # save additive gsolv calculated to
                                molecule.save_prev(
                                    exc[key], getattr(molecule, exc[key]).get("method")
                                )
                                # Gsolv for implicit solvation included in E
                                # save energy calculated to
                                molecule.save_prev(
                                    exc_implicit[key],
                                    getattr(molecule, exc_implicit[key]).get("method"),
                                )
                                # load new if available
                                # method naming -->
                                if key in ("smgsolv1", "smgsolv2"):
                                    func = config.func
                                    basis = config.basis
                                elif key in ("smgsolv3",):
                                    func = config.func3
                                    basis = config.basis3
                                if getattr(config, key) == "smd_gsolv":
                                    e_method, gsolv_method = config.get_method_name(
                                        "smd_gsolv",
                                        func=func,
                                        basis=basis,
                                        sm=getattr(config, key),
                                    )
                                elif getattr(config, key) in (
                                    "cosmors",
                                    "cosmors-fine",
                                ):
                                    e_method, gsolv_method = config.get_method_name(
                                        "cosmors",
                                        sm=getattr(config, key),
                                        func=func,
                                        basis=basis,
                                    )
                                elif getattr(config, key) in (
                                    "alpb_gsolv",
                                    "gbsa_gsolv",
                                ):
                                    e_method, gsolv_method = config.get_method_name(
                                        getattr(config, key),
                                        sm=getattr(config, key),
                                        func=func,
                                        basis=basis,
                                        gfn_version=getattr(config, exc2[key]),
                                    )
                                elif getattr(config, key) in (
                                    "cpcm",
                                    "cosmo",
                                    "smd",
                                    "dcosmors",
                                ):
                                    # Gsolv for implicit solvation included in E
                                    # need to reset gsolv (--> gsolv_method2 has to be "incl. in E")
                                    e_method, gsolv_method = config.get_method_name(
                                        "sp_implicit",
                                        func=func,
                                        basis=basis,
                                        sm=getattr(config, key),
                                    )
                                else:
                                    print("UNEXPECTED")
                                    e_method = ""
                                    gsolv_method = ""
                                molecule.load_prev(exc_implicit[key], e_method)
                                molecule.load_prev(exc[key], gsolv_method)
                        elif value and key in (
                            "func_or",
                            "basis_or",
                            "freq_or",
                            "func_or_scf",
                        ):
                            # save calculated to
                            molecule.save_prev(
                                "optical_rotation_info",
                                getattr(molecule, "optical_rotation_info").get(
                                    "method"
                                ),
                            )
                            # load new if available
                            method, _ = config.get_method_name(
                                "opt-rot",
                                prog=config.prog,
                                basis=config.basis_or,
                                func=config.func_or,
                                solvent=config.solvent,
                                sm="cosmo",
                                func2=config.func_or_scf,
                            )
                            molecule.load_prev("optical_rotation_info", method)
                        #####
                        elif value and key in ("func_j", "basis_j", "sm4_j"):
                            # reset to default, load is not available
                            molecule.nmr_coupling_info = getattr(
                                MoleculeData(0), "nmr_coupling_info"
                            )
                        elif value and key in ("func_s", "basis_s", "sm4_s"):
                            # reset to default, load is not available
                            molecule.nmr_shielding_info = getattr(
                                MoleculeData(0), "nmr_shielding_info"
                            )
                        #####
                        elif value and key in ("func3", "basis3"):
                            # save calculated to
                            molecule.save_prev(
                                "highlevel_sp_info",
                                getattr(molecule, "highlevel_sp_info").get("method"),
                            )
                            # load new if available
                            if config.smgsolv3 in ("cpcm", "smd", "cosmo", "dcosmors"):
                                expected = "sp_implicit"
                            else:
                                expected = "sp"
                            method, _ = config.get_method_name(
                                expected,
                                prog=config.prog3,
                                basis=config.basis3,
                                func=config.func3,
                                sm=config.smgsolv3,
                            )
                            molecule.load_prev("highlevel_sp_info", method)
                    # finally add molecule to list
                    conformers.append(molecule)
            # if nconf is increased add new conformers!
            newconfs = []
            for i in range(1, config.nconf + 1):
                considered = False
                for conf in list(conformers):
                    if conf.id == i:
                        considered = True
                        break
                if not considered:
                    print(f"Adding CONF{i} as new conformer!")
                    newconfs.append(QmJob(i))
            if newconfs:
                conformers.extend(newconfs)
                get_energy_from_ensemble(config.ensemblepath, config, conformers)
        elif not args.checkinput:
            # enso.json does not exist, create new conformers
            # don't create enso.json on checkinput
            # this is essentially the first start of CENSO
            print(
                f"{'INFORMATION:':{WARNLEN}}No restart information exists and is created during this run!\n"
            )
            conformers = []
            for i in range(1, config.nconf + 1):
                conformers.append(QmJob(i))
            # read energy from input_file and calculate rel_energy
            get_energy_from_ensemble(config.ensemblepath, config, conformers)
            # set fixed temperature for this and all following runs
            config._set_fixed_temperature()
            ensembledata = EnsembleData()
            ensembledata.filename = args.inp
            ensembledata.nconfs_per_part["starting"] = config.nconf
            config.write_json(
                config.cwd,
                [i.provide_runinfo() for i in conformers] + [ensembledata],
                config.provide_runinfo(),
            )

    if (args.checkinput and not error_logical) or (args.debug and args.checkinput):
        print("\nInput check is finished. The ENSO program can be executed.\n")
        sys.exit(0)
    if not conformers:
        print(f"{'ERROR:':{WARNLEN}}No conformers are considered!")
        print("\nGoing to exit!")
        sys.exit(1)

    # formatting information:
    config.lenconfx = max([len(str(i.id)) for i in conformers])
    conformers.sort(key=lambda x: int(x.id))
    return args, config, conformers, ensembledata
