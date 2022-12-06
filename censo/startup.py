"""
Contains enso_startup for the initialization of all parameters set for the
subsequent calculation.
"""
from argparse import Namespace
import os
import sys

from censo.core import CensoCore
from censo.cfg import (
    PLENGTH,
    WARNLEN,
    __version__,
)
from censo.datastructure import MoleculeData
from censo.qm_job import QmJob
from censo.utilities import (
    do_md5,
    t2x,
    move_recursively,
    get_energy_from_ensemble,
    frange,
    print,
    timeit
)
from censo.ensembledata import EnsembleData

core: CensoCore
cwd: str
args: Namespace

@timeit
@CensoCore.check_instance
def enso_startup():
    """
    4) read or write flags.dat control file
    5) check for crest_conformers.xyz
    6) check program settings
    7) read or write enso.json
    """

    global core, cwd, args, part_settings

    core = CensoCore.core()
    cwd = core.cwd
    args = core.args

    # check settings-combination and show error:
    # TODO - show all errors collectively
    try:
        core.internal_settings.check_logic()
    except KeyboardInterrupt:
        sys.exit(1)
    except Exception as error:
        if config.save_errors:
            print("")
            print("".ljust(PLENGTH, "*"))
            for _ in list(config.save_errors):
                print(config.save_errors.pop(0))
            print("".ljust(PLENGTH, "*"))
            print("")
            print("")
            print("Resulting python error:", error)
        else:
            print(error)
        sys.exit(1)


    """ if getattr(args, "onlyread", None) is not None:
        if getattr(args, "onlyread") == "on":
            setattr(config, "onlyread", True)
        elif getattr(args, "onlyread") == "off":
            setattr(config, "onlyread", False) """

    # printing parameters
    config.print_parameters(cmlcall=sys.argv)

    # TODO - check for the paths of needed programs
    
    """ # print errors
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
    # END print error """

    # TODO - debug stuff
    """ if error_logical and not args.debug and args.checkinput:
        print(
            f"\n{'ERROR:':{WARNLEN}}CENSO can not continue due to input errors!\n"
            f"{'':{WARNLEN}}Fix errors and run censo -checkinput again!"
        )
        print("\nGoing to exit!")
        sys.exit(1)
    elif error_logical and not args.debug:
        print(f"\n{'ERROR:':{WARNLEN}}CENSO can not continue due to input errors!")
        print("\nGoing to exit!")
        sys.exit(1) """

    # BEGIN - normal run without errors
    #if not error_logical or args.debug:
    if True:
        print("\n" + "".ljust(PLENGTH, "-"))
        print(" Processing data from previous run (enso.json)".center(PLENGTH, " "))
        print("".ljust(PLENGTH, "-") + "\n")
        # BEGIN - restart stuff
        # read enso.json
        """ if os.path.isfile(os.path.join(config.cwd, "enso.json")):
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
                    f"{'INFORMATION:':{WARNLEN}}The temperature has been changed "
                    f"from (previous) {prev_temp} K to (current) {getattr(config, 'temperature')} K."
                )
                # check multitemp on in previous
                if getattr(previousrun, "multitemp", None):
                    # has trange been changed ?
                    prev_trange = getattr(previousrun, "trange", [])
                    if not isinstance(prev_trange, list) or prev_trange != getattr(
                        config, "trange"
                    ):
                        print(
                            f"{'ERROR:':{WARNLEN}}The temperature range has been "
                            f"changed together with the temperature which is not allowed."
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
                            f"{'ERROR:':{WARNLEN}}Temperature has not been previously "
                            f"calculated, e.g. not in temperature range (trange). "
                            f"Change of temperature not possible!"
                        )
                        print(
                            f"{'INFORMATION:':{WARNLEN}}Temperatures in trange "
                            f"are: {[i for i in frange(prev_trange[0], prev_trange[1], prev_trange[2])]}"
                        )
                        print("\nGoing to exit!")
                        sys.exit(1)
                else:
                    print(
                        f"{'ERROR:':{WARNLEN}}Temperature change is requested, "
                        f"but no previous other temperatures have been calculated!"
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
                    if flag == 'basis' and 'automatic' in (getattr(config, flag, 'None'), getattr(previousrun, flag, 'None')):
                        print(
                            f"{'INFORMATION:':{WARNLEN}}setting {flag} was changed from "
                            f"{getattr(config, flag, 'None')} to {getattr(previousrun, flag, 'None')}!"
                        )                 
                    else:
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
                        # TODO - add UV/Vis
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
                get_energy_from_ensemble(config.ensemblepath, config, conformers) """
        # END - restart stuff
        
        # was elif:
        if not args.checkinput:
            # if input not to be checked dump conformer info etc. into new json
            # enso.json does not exist, create new conformers
            # don't create enso.json on checkinput
            # this is essentially the first start of CENSO
            print(
                f"{'INFORMATION:':{WARNLEN}}No restart information exists and is created during this run!\n"
            )
            conformers = []
            # create QmJobs of rank i
            for i in range(1, core.internal_settings.runinfo["nconf"] + 1):
                conformers.append(QmJob(i))
            # read energy from input_file and calculate rel_energy
            get_energy_from_ensemble(
                core.ensemble_path, 
                core.internal_settings.runinfo["nat"], 
                core.internal_settings.runinfo["nconf"], 
                core.internal_settings.runinfo["maxconf"],
                conformers
            )
            
            ensembledata = EnsembleData()
            ensembledata.filename = args.inp
            ensembledata.nconfs_per_part["starting"] = core.internal_settings.runinfo["nconf"]
            
            # TODO
            try:
                core.write_json(
                    [i.provide_runinfo() for i in conformers] + [ensembledata]
                )
            except Exception:
                """something TODO"""

    if args.checkinput:
        print("\nInput check is finished. The CENSO program can be executed.\n")
        sys.exit(0)
        
    if not conformers:
        print(f"{'ERROR:':{WARNLEN}}No conformers are considered!")
        print("\nGoing to exit!")
        sys.exit(1)

    # formatting information:
    config.lenconfx = max([len(str(i.id)) for i in conformers])
    conformers.sort(key=lambda x: int(x.id))
    return conformers, ensembledata