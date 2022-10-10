"""
REFINEMENT == Part3
designed to yield high level free energies on any conformers (DFT or non-DFT 
optimized).
"""
from multiprocessing import JoinableQueue as Queue
import shutil
import os
import sys
from .cfg import PLENGTH, CODING, AU2KCAL, DIGILEN, WARNLEN, qm_prepinfo
from .utilities import (
    check_for_folder,
    print_block,
    new_folders,
    last_folders,
    frange,
    calc_boltzmannweights,
    printout,
    move_recursively,
    write_trj,
    check_tasks,
    print,
    print_errors,
    ensemble2coord,
    conf_in_interval,
)
from .parallel import run_in_parallel
from .orca_job import OrcaJob
from .tm_job import TmJob
from .datastructure import MoleculeData


def part3(config, conformers, store_confs, ensembledata):
    """
    Refinement of the ensemble, at high level DFT (possibly with implicit solvation)
    Calculate low level free energies with COSMO-RS single-point and gsolv
    Input:
    - config [config_setup object] contains all settings
    - conformers [list of molecule_data objects] each conformer is represented
    Return:
    -> config
    -> conformers
    """
    save_errors = []
    if config.progress:
        print("#>>># CENSO: Starting part3", file=sys.stderr)
    print("\n" + "".ljust(PLENGTH, "-"))
    print("CRE REFINEMENT - PART3".center(PLENGTH, " "))
    print("".ljust(PLENGTH, "-") + "\n")
    # print flags for part3
    info = []
    info.append(["prog3", "program for part3"])
    info.append(["part3_threshold", "Boltzmann sum threshold G_thr(3)"])
    info.append(["func3", "functional for part3"])
    info.append(["basis3", "basis set for part3"])
    if config.solvent != "gas":
        info.append(["solvent", "solvent"])
        info.append(["smgsolv3", "solvent model"])
    info.append(["temperature", "temperature"])
    if config.multitemp:
        info.append(["multitemp", "evalulate at different temperatures"])
        info.append(
            [
                "printoption",
                "temperature range",
                [
                    i
                    for i in frange(
                        config.trange[0], config.trange[1], config.trange[2]
                    )
                ],
            ]
        )
    info.append(["prog_rrho", "program for mRRHO contribution"])
    if config.prog_rrho == "xtb":
        info.append(["part3_gfnv", "GFN version for mRRHO and/or GBSA_Gsolv"])
        if config.bhess:
            info.append(
                ["bhess", "Apply constraint to input geometry during mRRHO calculation"]
            )

    max_len_digilen = 0
    for item in info:
        if item[0] == "justprint":
            if "short-notation" in item[1]:
                tmp = len(item[1]) - len("short-notation:")
            else:
                tmp = len(item[1])
        else:
            tmp = len(item[1])
        if tmp > max_len_digilen:
            max_len_digilen = tmp
    max_len_digilen += 1
    if max_len_digilen < DIGILEN:
        max_len_digilen = DIGILEN

    optionsexchange = {True: "on", False: "off"}
    for item in info:
        if item[0] == "justprint":
            print(item[1:][0])
        else:
            if item[0] == "printoption":
                option = item[2]
            else:
                option = getattr(config, item[0])
            if option is True or option is False:
                option = optionsexchange[option]
            elif isinstance(option, list):
                option = [str(i) for i in option]
                if len(str(option)) > 40:
                    length = 0
                    reduced = []
                    for i in option:
                        length += len(i) + 2
                        if length < 40:
                            reduced.append(i)
                    reduced.append("...")
                    option = reduced
                    length = 0
                option = ", ".join(option)
            print(
                "{}: {:{digits}} {}".format(
                    item[1], "", option, digits=max_len_digilen - len(item[1])
                )
            )
    print("")
    # end print

    # TODO - fix the mess
    calculate = []  # has to be calculated in this run, list of QmJobs
    prev_calculated = []  # was already calculated in a previous run
    try:
        store_confs
    except NameError:
        store_confs = []  # stores all confs which are sorted out!

    # setup queues
    q = Queue()
    resultq = Queue()

    if config.prog3 == "tm":
        job = TmJob
    elif config.prog3 == "orca":
        job = OrcaJob

    for conf in list(conformers):
        if conf.removed:
            store_confs.append(conformers.pop(conformers.index(conf)))
            print(f"CONF{conf.id} is removed as requested by the user!")
            continue
        if conf.id > config.nconf:
            store_confs.append(conformers.pop(conformers.index(conf)))
            continue
        if conf.highlevel_sp_info["info"] == "failed":
            conf = conformers.pop(conformers.index(conf))
            store_confs.append(conf)
            print(f"Calculation of CONF{conf.id} failed in the previous run!")
        elif conf.highlevel_sp_info["info"] == "not_calculated":
            # has to be calculated now
            conf = conformers.pop(conformers.index(conf))
            calculate.append(conf)
        elif conf.highlevel_sp_info["info"] == "prep-failed":
            # has to be retried now
            conf = conformers.pop(conformers.index(conf))
            calculate.append(conf)
        elif conf.highlevel_sp_info["info"] == "calculated":
            conf = conformers.pop(conformers.index(conf))
            if config.solvent != "gas":
                # check if solvation calculation is calculated as well
                if conf.highlevel_gsolv_info["info"] == "failed":
                    store_confs.append(conf)
                    print(
                        f"Calculation of the solvation contribution for CONF"
                        f"{conf.id} failed in the previous run!"
                    )
                elif conf.highlevel_gsolv_info["info"] == "not_calculated":
                    calculate.append(conf)
                elif conf.highlevel_gsolv_info["info"] == "prep-failed":
                    # retry
                    calculate.append(conf)
                elif conf.highlevel_gsolv_info["info"] == "calculated":
                    conf.job["success"] = True
                    prev_calculated.append(conf)
                else:
                    print("UNEXPECTED BEHAVIOUR")
            elif config.solvent == "gas":
                conf.job["success"] = True
                prev_calculated.append(conf)
    if not calculate and not prev_calculated:
        print(f"{'ERROR:':{WARNLEN}}No conformers left!")
        print("Going to exit!")
        sys.exit(1)

    geometries_from_input = False
    # check if calculated on unoptimized geometries:
    if any(
        [
            conf.optimization_info["info"] == "not_calculated"
            for conf in calculate + prev_calculated
        ]
    ):
        if config.part2:
            print_errors(
                f"{'ERROR:':{WARNLEN}}Calculating (free) energies on "
                f"DFT unoptimized geometries!\n"
                f"{'':{WARNLEN}}Even though part2 is calculated!\n"
                f"{'':{WARNLEN}}Calculation on mixture of optimized "
                f"and unoptimized geometries is not advised!",
                save_errors,
            )
            print("Going to exit!")
            sys.exit(1)
        else:
            print_errors(
                f"{'WARNING:':{WARNLEN}}Calculating (free) energies "
                f"on DFT unoptimized geometries!",
                save_errors,
            )
            geometries_from_input = True
    # SI geometry:
    if geometries_from_input:
        ensembledata.si["part3"]["Geometry"] = "GFNn-xTB (input geometry)"
    else:
        ensembledata.si["part3"]["Geometry"], _ = config.get_method_name(
            "xtbopt",
            func=config.func,
            basis=config.basis,
            solvent=config.solvent,
            sm=config.sm2,
        )
        ensembledata.si["part3"]["Geometry"] += f" @optlevel: {config.optlevel2}"
    #

    if config.solvent == "gas":
        print("\nCalculating single-point energies!")
    else:
        print(
            "\nCalculating single-point energies and solvation contribution (G_solv)!"
        )

    instruction = {
        "func": config.func3,
        "basis": getattr(config, "basis3", "def2-TZVPD"),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": config.smgsolv3,
        "temperature": config.temperature,
        "gfn_version": config.part3_gfnv,
        "copymos": "",
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
        "onlyread": config.onlyread,
    }
    tmp_SI = None # ???
    if config.multitemp:
        instruction["trange"] = [
            i for i in frange(config.trange[0], config.trange[1], config.trange[2])
        ]
    else:
        instruction["trange"] = []
    if config.solvent == "gas":
        instruction["jobtype"] = "sp"
        instruction["prepinfo"] = ["high"]
        instruction["method"], _ = config.get_method_name(
            "sp", func=instruction["func"], basis=instruction["basis"]
        )
        folder = "part3"
        name = "highlevel single-point"
    else:
        if config.smgsolv3 in config.smgsolv_3:
            # additive GSolv
            # COSMORS
            if "cosmors" in config.smgsolv3 and config.smgsolv3 != "dcosmors":
                job = TmJob
                exc_fine = {"cosmors": "normal", "cosmors-fine": "fine"}
                tmp = {
                    "jobtype": "cosmors",
                    "prepinfo": ["high"],
                    "cosmorssetup": config.external_paths["cosmorssetup"],
                    "cosmorsparam": exc_fine.get(config.smgsolv3, "normal"),
                    "ctd-param": config.cosmorsparam,
                }
                if config.vapor_pressure:
                    tmp["vapor_pressure"] = True
                instruction.update(tmp)
                instruction["method"], instruction["method2"] = config.get_method_name(
                    "cosmors",
                    func=instruction["func"],
                    basis=instruction["basis"],
                    sm=instruction["sm"],
                )
                folder = "part3/COSMO"
                name = "highlevel COSMO-RS"
            # GBSA-Gsolv / ALPB-Gsolv
            elif instruction["sm"] in ("gbsa_gsolv", "alpb_gsolv"):
                # do DFT gas phase sp and additive Gsolv
                instruction["jobtype"] = instruction["sm"]
                instruction["method"], instruction["method2"] = config.get_method_name(
                    instruction["jobtype"],
                    func=instruction["func"],
                    basis=instruction["basis"],
                    sm=instruction["sm"],
                    gfn_version=instruction["gfn_version"],
                )
                if (
                    conf.highlevel_sp_info["info"] == "calculated"
                    and conf.highlevel_sp_info["method"] == instruction["method"]
                ):
                    # do not calculate gas phase sp again!
                    tmp_SI = ["high"]
                    instruction["prepinfo"] = []
                else:
                    instruction["prepinfo"] = ["high"]
                name = "highlevel " + str(instruction["sm"]).upper()
                folder = "part3"
            # SMD-Gsolv
            elif instruction["sm"] == "smd_gsolv":
                job = OrcaJob
                instruction["prepinfo"] = ["high"]
                instruction["jobtype"] = instruction["sm"]
                instruction["method"], instruction["method2"] = config.get_method_name(
                    "smd_gsolv",
                    func=instruction["func"],
                    basis=instruction["basis"],
                    sm=instruction["sm"],
                )
                name = "highlevel" + str(instruction["sm"]).upper()
                folder = "part3"
        else:
            # with implicit solvation
            instruction["jobtype"] = "sp_implicit"
            instruction["prepinfo"] = ["high"]
            instruction["method"], instruction["method2"] = config.get_method_name(
                "sp_implicit",
                func=instruction["func"],
                basis=instruction["basis"],
                sm=instruction["sm"],
            )
            name = "high level single-point"
            folder = "part3"

    if prev_calculated:
        check_for_folder(config.cwd, [i.id for i in prev_calculated], folder)
        print("The calculation was performed before for:")
        print_block(["CONF" + str(i.id) for i in prev_calculated])
    pl = config.lenconfx + 4 + len(str("/" + folder))

    check = {True: "was successful", False: "FAILED"}
    if geometries_from_input:
        # geoms folder needed for coord best and trj
        geoms_folder = folder
    if calculate:
        # make new folder:
        save_errors, store_confs, calculate = new_folders(
            config.cwd, calculate, folder, save_errors, store_confs
        )
        # write coord to folder
        if config.smgsolv3 in ("cosmors", "cosmors-fine"):
            cp_to = ["part3", folder]
        else:
            cp_to = [folder]
        if geometries_from_input:
            print_errors(
                f"{'INFORMATION:':{WARNLEN}}Taking geometries from ensemble input (unoptimized)!",
                save_errors,
            )
        for folder in cp_to:
            if geometries_from_input:
                # working on DFT unoptimized geometries
                calculate, store_confs, save_errors = ensemble2coord(
                    config, folder, calculate, store_confs, save_errors
                )
            else:
                # need to copy optimized coord to folder
                for conf in list(calculate):
                    tmp1 = os.path.join(
                        config.cwd, "CONF" + str(conf.id), config.func, "coord"
                    )
                    tmp2 = os.path.join("CONF" + str(conf.id), folder, "coord")
                    try:
                        shutil.copy(tmp1, tmp2)
                    except FileNotFoundError:
                        print_errors(
                            f"{'ERROR:':{WARNLEN}}can't copy optimized "
                            f"geometry of CONF{conf.id}!",
                            save_errors,
                        )
                        store_confs.append(calculate.pop(calculate.index(conf)))

        if config.solvent == "gas":
            print("The high level single-point is now calculated for:")
        else:
            print("The high level gsolv calculation is now calculated for:")
        print_block(["CONF" + str(i.id) for i in calculate])
        # parallel calculation:
        calculate = run_in_parallel(
            config,
            q,
            resultq,
            job,
            config.maxthreads,
            config.omp,
            calculate,
            instruction,
            config.balance,
            folder,
        )
        # check if too many calculations failed

        for conf in list(calculate):
            if instruction["jobtype"] in ("sp", "sp_implicit"):
                line = (
                    f"{name} calculation {check[conf.job['success']]}"
                    f" for {last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.highlevel_sp_info["info"] = "failed"
                    conf.highlevel_sp_info["method"] = instruction["method"]
                    conf.part_info["part3"] = "refused"
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.highlevel_sp_info["energy"] = conf.job["energy"]
                    conf.highlevel_sp_info["info"] = "calculated"
                    conf.highlevel_sp_info["method"] = instruction["method"]
                    if instruction["jobtype"] == "sp_implicit":
                        conf.highlevel_gsolv_info["range"] = {
                            conf.job["temperature"]: 0.0
                        }
                        conf.highlevel_gsolv_info["info"] = "calculated"
                        conf.highlevel_gsolv_info["method"] = instruction["method2"]
            elif instruction["jobtype"] in (
                "cosmors",
                "smd_gsolv",
                "gbsa_gsolv",
                "alpb_gsolv",
            ):
                line = (
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 3):>{pl}}: "
                    f"{conf.job['energy2']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.part_info["part3"] = "refused"
                    conf.highlevel_sp_info["info"] = "failed"
                    conf.highlevel_sp_info["method"] = instruction["method"]
                    conf.highlevel_gsolv_info["info"] = "failed"
                    conf.highlevel_gsolv_info["method"] = instruction["method2"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    if (
                        instruction["jobtype"] in ("gbsa_gsolv", "alpb_gsolv")
                        and conf.highlevel_sp_info["info"] == "calculated"
                        and conf.highlevel_sp_info["method"] == instruction["method"]
                    ):
                        conf.job["energy"] = conf.highlevel_sp_info["energy"]
                    conf.highlevel_sp_info["energy"] = conf.job["energy"]
                    conf.highlevel_sp_info["info"] = "calculated"
                    conf.highlevel_sp_info["method"] = instruction["method"]
                    conf.highlevel_gsolv_info["gas-energy"] = conf.job["energy"]
                    conf.highlevel_gsolv_info["info"] = "calculated"
                    conf.highlevel_gsolv_info["method"] = instruction["method2"]
                    conf.highlevel_gsolv_info["range"] = conf.job["erange1"]
            else:
                print(
                    f'UNEXPECTED BEHAVIOUR: {conf.job["success"]} {conf.job["jobtype"]}'
                )
        # save current data to jsonfile
        config.write_json(
            config.cwd,
            [i.provide_runinfo() for i in calculate]
            + [i.provide_runinfo() for i in prev_calculated]
            + [i.provide_runinfo() for i in store_confs]
            + [ensembledata],
            config.provide_runinfo(),
        )
        check_tasks(calculate, config.check)
    else:
        print("No conformers are considered additionally.")
    # adding conformers calculated before:
    if prev_calculated:
        # adding conformers calculated before:
        for conf in list(prev_calculated):
            conf.job["workdir"] = os.path.normpath(
                os.path.join(config.cwd, "CONF" + str(conf.id), folder)
            )
            if instruction["jobtype"] in ("sp", "sp_implicit"):
                print(
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.highlevel_sp_info['energy']:>.8f}"
                )
            elif instruction["jobtype"] in (
                "cosmors",
                "smd_gsolv",
                "gbsa_gsolv",
                "alpb_gsolv",
            ):
                print(
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 3):>{pl}}: "
                    f"{getattr(conf, 'highlevel_gsolv_info').get('range',{}).get(config.temperature, 0.0):>.8f}"
                )
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))

    # create data for possible SI generation:-----------------------------------
    # Energy:
    ensembledata.si["part3"]["Energy"] = instruction["method"]
    if "DOGCP" in instruction["prepinfo"]:
        ensembledata.si["part3"]["Energy"] += " + GCP"
    # Energy_settings:
    try:
        if job == TmJob:
            if tmp_SI is not None:
                tmp = " ".join(qm_prepinfo["tm"][tmp_SI[0]]).replace("-", "")
            else:
                tmp = " ".join(qm_prepinfo["tm"][instruction["prepinfo"][0]]).replace(
                    "-", ""
                )
        elif job == OrcaJob:
            if tmp_SI is not None:
                tmp = " ".join(qm_prepinfo["orca"][tmp_SI[0]])
            else:
                tmp = " ".join(qm_prepinfo["orca"][instruction["prepinfo"][0]])
        ensembledata.si["part3"]["Energy_settings"] = tmp
    except (TypeError, IndexError):
        pass
    # GmRRHO is done below!
    # Solvation:
    if config.solvent == "gas":
        ensembledata.si["part3"]["G_solv"] = "gas-phase"
    else:
        if instruction["jobtype"] in ("cosmors", "cosmors-fine", "cosmors-normal"):
            ensembledata.si["part3"]["G_solv"] = (
                instruction["method2"] + " param= " + instruction["ctd-param"]
            )
        else:
            ensembledata.si["part3"]["G_solv"] = instruction["method2"]
    # Geometry is done above:
    # QM-CODE:
    ensembledata.si["part3"]["main QM code"] = str(config.prog3).upper()
    # Threshold:
    ensembledata.si["part3"][
        "Threshold"
    ] = f"Boltzmann sum threshold: {config.part3_threshold} %"
    # END SI generation --------------------------------------------------------

    for conf in calculate:
        conf.reset_job_info()

    if not calculate:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)
    # ***************************************************************************
    if config.evaluate_rrho:
        instruction_rrho = {
            "jobtype": "rrhoxtb",
            "func": getattr(config, "part3_gfnv"),
            "gfn_version": getattr(config, "part3_gfnv"),
            "temperature": config.temperature,
            "charge": config.charge,
            "unpaired": config.unpaired,
            "solvent": config.solvent,
            "bhess": config.bhess,
            "sm_rrho": config.sm_rrho,
            "rmsdbias": config.rmsdbias,
            "cwd": config.cwd,
            "consider_sym": config.consider_sym,
            "energy": 0.0,
            "energy2": 0.0,
            "success": False,
            "progpath": config.external_paths["xtbpath"],
            "imagthr": config.imagthr,
            "sthr": config.sthr,
            "scale": config.scale,
            "onlyread": config.onlyread,
        }
        folder_rrho = "rrho_part3"
        instruction_rrho["method"], _ = config.get_method_name(
            "rrhoxtb",
            bhess=config.bhess,
            gfn_version=instruction_rrho["gfn_version"],
            sm=instruction_rrho["sm_rrho"],
            solvent=instruction_rrho["solvent"],
        )

        # GmRRHO for SI information:
        if config.evaluate_rrho:
            ensembledata.si["part3"]["G_mRRHO"] = instruction_rrho["method"]
            if "bhess" in ensembledata.si["part3"]["G_mRRHO"]:
                ensembledata.si["part3"]["G_mRRHO"] += " SPH"
        else:
            ensembledata.si["part3"]["G_mRRHO"] = "not included"
        # END SI

        if config.multitemp:
            instruction_rrho["trange"] = [
                i for i in frange(config.trange[0], config.trange[1], config.trange[2])
            ]
        else:
            instruction_rrho["trange"] = []

        # check if calculated
        for conf in list(calculate):
            if conf.removed:
                store_confs.append(calculate.pop(calculate.index(conf)))
                print(f"CONF{conf.id} is removed as requested by the user!")
                continue
            if conf.highlevel_grrho_info["info"] == "calculated":
                conf = calculate.pop(calculate.index(conf))
                conf.job["success"] = True
                conf.sym = conf.highlevel_grrho_info.get("sym", conf.sym)
                conf.linear = conf.highlevel_grrho_info.get("linear", conf.linear)
                conf.symnum = conf.highlevel_grrho_info.get("symnum", conf.symnum)
                prev_calculated.append(conf)
            elif conf.highlevel_grrho_info["info"] == "failed":
                conf = calculate.pop(calculate.index(conf))
                conf.part_info["part3"] = "refused"
                store_confs.append(conf)
                print(f"Calculation of CONF{conf.id} failed in the previous run!")
            elif conf.highlevel_grrho_info["info"] in ("not_calculated", "prep-failed"):
                # stay in calculate (e.g not_calculated or prep-failed)
                # check if method has been calculated in part2 and part2 is active in this run
                # (otherwise possible mixup of optimized and unoptimized geometries)
                if (
                    instruction_rrho["method"] == conf.lowlevel_grrho_info["method"]
                    and config.part2
                ):
                    # has been calculated before, just copy
                    conf.job["success"] = True
                    attributes = vars(MoleculeData(0)).get("highlevel_grrho_info")
                    tmp = {}
                    for key in attributes.keys():
                        if key != "prev_methods":
                            tmp[key] = getattr(conf, "lowlevel_grrho_info").get(key)
                    getattr(conf, "highlevel_grrho_info").update(tmp)
                    # if rrho is taken from part2 still make folder for rrho_part3
                    save_errors, store_confs, calculate = new_folders(
                        config.cwd,
                        [conf],
                        folder_rrho,
                        save_errors,
                        store_confs,
                        silent=True,
                    )
                    prev_calculated.append(calculate.pop(calculate.index(conf)))
                elif (
                    instruction_rrho["method"]
                    in conf.lowlevel_grrho_info["prev_methods"].keys()
                    and config.part2
                ):
                    # has been calculated before, just copy
                    conf.job["success"] = True
                    conf.load_prev(
                        "lowlevel_grrho_info",
                        instruction_rrho["method"],
                        saveto="highlevel_grrho_info",
                    )
                    # if rrho is taken from part2 still make folder for rrho_part3
                    save_errors, store_confs, calculate = new_folders(
                        config.cwd,
                        [conf],
                        folder_rrho,
                        save_errors,
                        store_confs,
                        silent=True,
                    )
                    prev_calculated.append(calculate.pop(calculate.index(conf)))
                else:
                    # stays in calculate and has to be calculated now
                    pass
            else:
                print("UNEXPECTED BEHAVIOUR")
        if not calculate and not prev_calculated:
            print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
            print("Going to exit!")
            sys.exit(1)
        # do the rrho stuff:
        if config.solvent == "gas":
            if geometries_from_input:
                print("\nCalculating high level G_mRRHO on input geometry!")
            else:
                print("\nCalculating high level G_mRRHO on DFT geometry!")
        else:
            if geometries_from_input:
                print(
                    "\nCalculating high level G_mRRHO with implicit solvation "
                    "on input geometry!"
                )
            else:
                print(
                    "\nCalculating high level G_mRRHO with implicit solvation "
                    "on DFT geometry!"
                )
        if prev_calculated:
            check_for_folder(config.cwd, [i.id for i in prev_calculated], folder_rrho)
            print("The G_mRRHO calculation was performed before for:")
            print_block(["CONF" + str(i.id) for i in prev_calculated])
        pl = config.lenconfx + 4 + len(str("/" + folder_rrho))

        if calculate:
            print("The G_mRRHO calculation is now performed for:")
            print_block(["CONF" + str(i.id) for i in calculate])
            # create folders:
            save_errors, store_confs, calculate = new_folders(
                config.cwd, calculate, folder_rrho, save_errors, store_confs
            )
            # write coord to folder
            if geometries_from_input:
                # working on DFT unoptimized geometries
                calculate, store_confs, save_errors = ensemble2coord(
                    config, folder_rrho, calculate, store_confs, save_errors
                )
            else:
                # copy optimized geoms to folder
                for conf in list(calculate):
                    try:
                        tmp_from = os.path.join(
                            config.cwd, "CONF" + str(conf.id), config.func
                        )
                        tmp_to = os.path.join(
                            config.cwd, "CONF" + str(conf.id), folder_rrho
                        )
                        shutil.copy(
                            os.path.join(tmp_from, "coord"),
                            os.path.join(tmp_to, "coord"),
                        )
                    except shutil.SameFileError:
                        pass
                    except FileNotFoundError:
                        if not os.path.isfile(os.path.join(tmp_from, "coord")):
                            print(
                                f"{'ERROR:':{WARNLEN}}while copying the coord file from {tmp_from}! "
                                "The corresponding file does not exist."
                            )
                        elif not os.path.isdir(tmp_to):
                            print(
                                f"{'ERROR:':{WARNLEN}}Could not create folder {tmp_to}!"
                            )
                        print(f"{'ERROR:':{WARNLEN}}Removing conformer {conf.name}!")
                        conf.highlevel_grrho_info["info"] = "prep-failed"
                        store_confs.append(calculate.pop(calculate.index(conf)))
                        save_errors.append(
                            f"CONF{conf.id} was removed, because IO failed!"
                        )
            # parallel execution (calculate: list of qm_jobs)
            calculate = run_in_parallel(
                config,
                q,
                resultq,
                job,
                config.maxthreads,
                config.omp,
                calculate,
                instruction_rrho,
                config.balance,
                folder_rrho,
            )
            check = {True: "was successful", False: "FAILED"}
            # check if too many calculations failed

            # TODO - grab uvvis results if calculated
            ###
            for conf in list(calculate):
                print(
                    f"The G_mRRHO calculation @ {conf.job['symmetry']} "
                    f"{check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f} "
                    f"""S_rot(sym)= {conf.calc_entropy_sym(
                        config.temperature,
                        symnum=conf.job['symnum']):>.7f} """
                    f"""using= {conf.get_mrrho(
                        config.temperature,
                        rrho='direct_input',
                        consider_sym=config.consider_sym,
                        symnum=conf.job['symnum'],
                        direct_input=conf.job["energy"]):>.7f}"""
                )
                if not conf.job["success"]:
                    conf.part_info["part3"] = "refused"
                    conf.highlevel_grrho_info["info"] = "failed"
                    conf.highlevel_grrho_info["method"] = instruction_rrho["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.sym = conf.job["symmetry"]
                    conf.linear = conf.job["linear"]
                    conf.highlevel_grrho_info["sym"] = conf.job["symmetry"]
                    conf.highlevel_grrho_info["linear"] = conf.job["linear"]
                    conf.highlevel_grrho_info["symnum"] = conf.job["symnum"]
                    conf.highlevel_grrho_info["rmsd"] = conf.job["rmsd"]
                    conf.highlevel_grrho_info["info"] = "calculated"
                    conf.highlevel_grrho_info["method"] = instruction_rrho["method"]
                    conf.highlevel_grrho_info["range"] = conf.job["erange1"]
                    conf.highlevel_hrrho_info["range"] = conf.job["erange2"]
                    conf.highlevel_hrrho_info["info"] = "calculated"
                    conf.highlevel_hrrho_info["method"] = instruction_rrho["method"]
                    conf.uvvis_info["nroots"] = conf.job["nroots"]
                    conf.uvvis_info["excitations"] = conf.job["excitations"]
            # save current data to jsonfile (get settings via config_setup.provide_runinfo())
            config.write_json(
                config.cwd,
                [i.provide_runinfo() for i in calculate]
                + [i.provide_runinfo() for i in prev_calculated]
                + [i.provide_runinfo() for i in store_confs]
                + [ensembledata],
                config.provide_runinfo(),
            )
            check_tasks(calculate, config.check)
        else:
            print("No conformers are considered additionally.")
        # adding conformers calculated before:
        if prev_calculated:
            for conf in list(prev_calculated):
                conf.job["workdir"] = os.path.normpath(
                    os.path.join(config.cwd, "CONF" + str(conf.id), folder_rrho)
                )
                print(
                    f"The G_mRRHO calculation @ {conf.sym} "
                    f"{check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"""{conf.get_mrrho(
                        config.temperature,
                        rrho='highlevel_grrho_info',
                        consider_sym=True):>.8f} """
                    f"S_rot(sym)= {conf.calc_entropy_sym(config.temperature):>.7f}"
                    f""" using= {conf.get_mrrho(
                        config.temperature,
                        rrho='highlevel_grrho_info',
                        consider_sym=config.consider_sym):>.7f}"""
                )
                calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
        if not calculate:
            print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
            print("Going to exit!")
            sys.exit(1)
        # save current data to jsonfile
        config.write_json(
            config.cwd,
            [i.provide_runinfo() for i in calculate]
            + [i.provide_runinfo() for i in prev_calculated]
            + [i.provide_runinfo() for i in store_confs]
            + [ensembledata],
            config.provide_runinfo(),
        )
    # printout for part3 -------------------------------------------------------
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("* Gibbs free energies of part3 *".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "xtb_energy"),
        lambda conf: getattr(conf, "rel_xtb_energy"),
        lambda conf: getattr(conf, "highlevel_sp_info")["energy"],
        lambda conf: getattr(conf, "highlevel_gsolv_info")
        .get("range", {})
        .get(config.temperature, 0.0),
        lambda conf: conf.get_mrrho(
            config.temperature, "highlevel_grrho_info", config.consider_sym
        ),
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "rel_free_energy"),
        lambda conf: getattr(conf, "bm_weight") * 100,
    ]
    columnheader = [
        "CONF#",
        "E(GFNn-xTB)",
        "ΔE(GFNn-xTB)",
        "E [Eh]",
        "Gsolv [Eh]",
        "GmRRHO [Eh]",
        "Gtot",
        "ΔGtot",
        "Boltzmannweight",
    ]
    columndescription = [
        "",
        "[a.u.]",
        "[kcal/mol]",
        "",
        "",
        "",
        "[Eh]",
        "[kcal/mol]",
        f"  % at {config.temperature:.2f} K",
    ]
    columndescription2 = ["", "", "", "", "", "", "", "", ""]
    columnformat = [
        "",
        (12, 7),
        (5, 2),
        (12, 7),
        (12, 7),
        (12, 7),
        (12, 7),
        (5, 3),
        (5, 2),
    ]
    if config.solvent == "gas":
        columndescription[3] = instruction["method"]
    elif config.solvent != "gas":
        columndescription[3] = instruction["method"]
        columndescription[4] = instruction["method2"]
    if config.evaluate_rrho:
        columndescription[5] = instruction_rrho["method"]
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho:
            # ignore rrho in printout
            columncall.pop(5)
            columnheader.pop(5)
            columndescription.pop(5)
            columndescription2.pop(5)
            columnformat.pop(5)
        if config.solvent == "gas":
            columncall.pop(4)
            columnheader.pop(4)
            columndescription.pop(4)
            columndescription2.pop(4)
            columnformat.pop(4)

    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "highlevel_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "highlevel_gsolv_info"
        e = "highlevel_sp_info"
        conf.calc_free_energy(
            e=e,
            solv=solv,
            rrho=rrho,
            t=config.temperature,
            consider_sym=config.consider_sym,
        )

    calculate = calc_boltzmannweights(calculate, "free_energy", config.temperature)
    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
    except ValueError:
        raise
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part3.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
        columndescription2=columndescription2,
    )
    # end printout for part3
    conf_in_interval(calculate)
    # --------------------------------------------------------------------------
    for conf in calculate:
        if conf.free_energy == minfree:
            ensembledata.bestconf["part3"] = conf.id

    ################################################################################
    # TODO - plot uvvis spectra
    # calculate average G correction
    print("\nCalculating Boltzmann averaged free energy of ensemble!\n")
    avGcorrection = {"avG": {}, "avE": {}, "avGsolv": {}, "avGRRHO": {}}
    if config.multitemp:
        trange = [
            i for i in frange(config.trange[0], config.trange[1], config.trange[2])
        ]
    else:
        trange = [config.temperature]
    # calculate Boltzmannweights
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho and config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avG(T) /a.u.':>14} "
                f"{'avGcorrection(T)* /a.u.':>22}"
            )
        elif not config.evaluate_rrho:
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGsolv(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
                f"{'avGcorrection(T)* /a.u.':>22}"
            )
        elif config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGmRRHO(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
                f"{'avGcorrection(T)* /a.u.':>22}"
            )
    else:
        line = (
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
            f"{'avGmRRHO(T) /a.u.':>16} {'avGsolv(T) /a.u.':>16} "
            f"{'avG(T) /a.u.':>14}"
            f" {'avGcorrection(T)* /a.u.':>22}"
        )
    print(line)
    print("".ljust(int(PLENGTH), "-"))
    for temperature in trange:
        # get free energy at (T)
        for conf in calculate:
            if not config.evaluate_rrho:
                rrho = None
            else:
                rrho = "highlevel_grrho_info"
            if config.solvent == "gas":
                solv = None
            else:
                solv = "highlevel_gsolv_info"
            e = "highlevel_sp_info"
            conf.calc_free_energy(
                e=e,
                solv=solv,
                rrho=rrho,
                t=temperature,
                consider_sym=config.consider_sym,
            )

        calculate = calc_boltzmannweights(calculate, "free_energy", temperature)
        avG = 0.0
        avE = 0.0
        avGRRHO = 0.0

        avGsolv = 0.0
        for conf in calculate:
            avG += conf.bm_weight * conf.free_energy
            avE += conf.bm_weight * conf.highlevel_sp_info["energy"]
            avGRRHO += conf.bm_weight * conf.get_mrrho(
                temperature, "highlevel_grrho_info", config.consider_sym
            )
            avGsolv += conf.bm_weight * conf.highlevel_gsolv_info["range"].get(
                temperature, 0.0
            )

        avGcorrection["avG"][temperature] = avG
        avGcorrection["avE"][temperature] = avE
        avGcorrection["avGRRHO"][temperature] = avGRRHO
        avGcorrection["avGsolv"][temperature] = avGsolv
        # printout:
        if not config.evaluate_rrho or config.solvent == "gas":
            if not config.evaluate_rrho and config.solvent == "gas":
                line = (
                    f"{temperature:^15} {avE:>14.7f}  {avG:>14.7f} "
                    f"{ensembledata.avGcorrection.get('avGcorrection',{}).get(temperature, 0.0):>22.7f}"
                )
            elif not config.evaluate_rrho:
                line = (
                    f"{temperature:^15} {avE:>14.7f} {avGsolv:>16.7f} "
                    f"{avG:>14.7f} "
                    f"{ensembledata.avGcorrection.get('avGcorrection',{}).get(temperature, 0.0):>22.7f}"
                )
            elif config.solvent == "gas":
                line = (
                    f"{temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                    f"{avG:>14.7f} "
                    f"{ensembledata.avGcorrection.get('avGcorrection',{}).get(temperature, 0.0):>22.7f}"
                )
        else:
            line = (
                f"{temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                f"{avGsolv:>16.7f} {avG:>14.7f} "
                f"{ensembledata.avGcorrection.get('avGcorrection',{}).get(temperature, 0.0):>22.7f}"
            )
        if temperature == config.temperature:
            print(line, "    <<==part3==")
        else:
            print(line)
    print("".ljust(int(PLENGTH), "-"))
    print(
        "avGcorrection(T)* = Correction to free energy, which can be added (manually) to avG(T)."
    )
    print("                    If only a small ensemble is evaluated in part3 ")
    print("                    this can incorporate information of higher lying ")
    print("                    conformers of a more complete ensemble from part2.")
    print("")

    ################################################################################
    # -----------------------------Trange Ouput END------------------------------
    # reset boltzmannweights to correct temperature
    # get free energy at (T)
    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "highlevel_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "highlevel_gsolv_info"
        e = "highlevel_sp_info"
        conf.calc_free_energy(
            e=e,
            solv=solv,
            rrho=rrho,
            t=config.temperature,
            consider_sym=config.consider_sym,
        )
    calculate = calc_boltzmannweights(calculate, "free_energy", config.temperature)

    # SORTING for the next part:
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("Conformers considered further".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    # evaluate conformer consideration based on Boltzmann-population
    calculate.sort(reverse=True, key=lambda x: float(x.bm_weight))
    sumup = 0.0
    for conf in list(calculate):
        sumup += conf.bm_weight
        if sumup >= (config.part3_threshold / 100):
            if conf.bm_weight < (1 - (config.part3_threshold / 100)):
                mol = calculate.pop(calculate.index(conf))
                mol.part_info["part3"] = "refused"
                store_confs.append(mol)
            else:
                conf.part_info["part3"] = "passed"
        else:
            conf.part_info["part3"] = "passed"

    ensembledata.nconfs_per_part["part3"] = len(calculate)

    if calculate:
        print(
            f"\nConformers that are below the Boltzmann threshold G_thr(3) of {config.part3_threshold}%:"
        )
        print_block(["CONF" + str(i.id) for i in calculate])

    # save current data to jsonfile
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    if geometries_from_input:
        folder = geoms_folder
    else:
        folder = config.func

    # write ensemble
    move_recursively(config.cwd, "enso_ensemble_part3.xyz")
    kwargs = {"energy": "xtb_energy", "rrho": "highlevel_grrho_info"}
    write_trj(
        sorted(calculate, key=lambda x: float(x.free_energy)),
        config.cwd,
        "enso_ensemble_part3.xyz",
        folder,
        config.nat,
        "free_energy",
        config.temperature,
        config.consider_sym,
        **kwargs,
    )

    # write coord.enso_best
    for conf in calculate:
        if conf.id == ensembledata.bestconf["part3"]:
            # copy the lowest conformer to file coord.enso_best
            with open(
                os.path.join("CONF" + str(conf.id), folder, "coord"),
                "r",
                encoding=CODING,
                newline=None,
            ) as f:
                coord = f.readlines()
            with open(
                os.path.join(config.cwd, "coord.enso_best"), "w", newline=None
            ) as best:
                best.write(
                    "$coord  # {}   {}   !CONF{} \n".format(
                        conf.free_energy,
                        conf.get_mrrho(
                            config.temperature,
                            rrho="highlevel_grrho_info",
                            consider_sym=config.consider_sym,
                        ),
                        conf.id,
                    )
                )
                for line in coord[1:]:
                    if "$" in line:  # stop at $end ...
                        break
                    best.write(line)
                best.write("$end \n")

    # reset
    for conf in calculate:
        conf.free_energy = 0.0
        conf.rel_free_energy = 0.0
        conf.bm_weight = 0.0
        conf.reset_job_info()

    if save_errors:
        print("\n***---------------------------------------------------------***")
        print("Printing most relevant errors again, just for user convenience:")
        for _ in list(save_errors):
            print(save_errors.pop())
        print("***---------------------------------------------------------***")
    tmp = int((PLENGTH - len("END of Part3")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part3" + "".rjust(tmp, "<"))
    if config.progress:
        print("#>>># CENSO: Finished part3", file=sys.stderr)
    return config, calculate, store_confs, ensembledata
