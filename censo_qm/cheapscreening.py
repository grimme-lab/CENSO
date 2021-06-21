"""
prescreening == part0, calculate cheap free energy on GFNn-xTB input geometry
The idea is to improve on the description of E with a very fast DFT method.
"""
import os
import sys
from multiprocessing import JoinableQueue as Queue
from .cfg import PLENGTH, DIGILEN, AU2KCAL, WARNLEN, qm_prepinfo, dfa_settings
from .parallel import run_in_parallel
from .orca_job import OrcaJob
from .tm_job import TmJob
from .utilities import (
    check_for_folder,
    print_block,
    new_folders,
    last_folders,
    ensemble2coord,
    printout,
    check_tasks,
    print,
    print_errors,
    calc_boltzmannweights,
    conf_in_interval,
)


def part0(config, conformers, ensembledata):
    """
    Cheap prescreening of the ensemble, with single-points on combined ensemble
    geometries.
    Input:
    - config [config_setup object] contains all settings
    - conformers [list of molecule_data objects] each conformer is represented
    - ensembledata -> instance for saving ensemble (not conf) related data
    Return:
    -> config
    -> conformers
    -> store_confs
    """
    save_errors = []
    if config.progress:
        print("#>>># CENSO: Starting part0", file=sys.stderr)
    print("\n" + "".ljust(PLENGTH, "-"))
    print("CRE CHEAP-PRESCREENING - PART0".center(PLENGTH, " "))
    print("".ljust(PLENGTH, "-") + "\n")
    # print flags for part1
    info = []
    info.append(["prog", "program"])
    info.append(["func0", "functional for part0"])
    info.append(["basis0", "basis set for part0"])
    info.append(["part0_threshold", "threshold g_thr(0)"])
    info.append(["nconf", "starting number of considered conformers"])
    info.append(["printoption", "temperature", config.fixed_temperature])

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
                option = ", ".join(option)
            print(
                "{}: {:{digits}} {}".format(
                    item[1], "", option, digits=max_len_digilen - len(item[1])
                )
            )
    print("")
    # end print

    calculate = []  # has to be calculated in this run
    prev_calculated = []  # was already calculated in a previous run
    store_confs = []  # stores all confs which are sorted out!

    if config.temperature != config.fixed_temperature:
        print(
            f"{'INFORMATION:':{WARNLEN}}The temperature in part0 is kept fixed at {config.fixed_temperature} to keep consistent"
        )
        print(f"{'':{WARNLEN}}data in case of restarts.\n")

    print("Calculating efficient gas-phase single-point energies:")

    # setup queues
    q = Queue()
    resultq = Queue()

    folder = "part0_sp"

    if config.prog == "tm":
        job = TmJob
    elif config.prog == "orca":
        job = OrcaJob

    for conf in list(conformers):
        conf = conformers.pop(conformers.index(conf))
        if conf.removed:
            store_confs.append(conf)
            print(f"CONF{conf.id} is removed as requested by the user.")
            continue
        if conf.id > config.nconf:
            store_confs.append(conf)
            continue
        if conf.cheap_prescreening_sp_info["info"] == "not_calculated":
            calculate.append(conf)
        elif conf.cheap_prescreening_sp_info["info"] == "failed":
            store_confs.append(conf)
            print(f"Calculation of CONF{conf.id} failed in the previous run.")
        elif conf.cheap_prescreening_sp_info["info"] == "calculated":
            conf.job["success"] = True
            prev_calculated.append(conf)
        else:
            print(f"{'ERROR:':{WARNLEN}}UNEXPECTED BEHAVIOUR")
    if not calculate and not prev_calculated:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
    if prev_calculated:
        check_for_folder(config.cwd, [i.id for i in prev_calculated], folder)
        print("The efficient gas-phase single-point was calculated before for:")
        print_block(["CONF" + str(i.id) for i in prev_calculated])
    pl = config.lenconfx + 4 + len(str("/" + folder))

    if config.solvent != "gas":
        instruction = {
            "jobtype": "alpb_gsolv",
            "func": config.func0,
            "basis": getattr(
                config,
                "basis0",
                dfa_settings.composite_method_basis.get(config.func0, "def2-SV(P)"),
            ),
            "charge": config.charge,
            "unpaired": config.unpaired,
            "solvent": config.solvent,
            "sm": config.sm_rrho,
            "gfn_version": config.part0_gfnv,
            "energy": 0.0,
            "energy2": 0.0,
            "success": False,
            "temperature": config.fixed_temperature,
            "onlyread": config.onlyread,
        }

    elif config.solvent == "gas":
        instruction = {
            "jobtype": "sp",
            "func": config.func0,
            "basis": getattr(
                config,
                "basis0",
                dfa_settings.composite_method_basis.get(config.func0, "def2-SV(P)"),
            ),
            "charge": config.charge,
            "unpaired": config.unpaired,
            "solvent": "gas",
            "sm": "gas-phase",
            "energy": 0.0,
            "energy2": 0.0,
            "success": False,
            "onlyread": config.onlyread,
        }

    if config.prog == "tm":
        instruction["prepinfo"] = ["clear", "-grid", "1", "-scfconv", "5", "DOGCP"]
        ensembledata.si["part0"]["Energy_settings"] = "scfconv 5, grid 1"

    elif config.prog == "orca":
        instruction["prepinfo"] = ["low", "DOGCP"]
        ensembledata.si["part0"]["Energy_settings"] = " ".join(
            qm_prepinfo["orca"][instruction["prepinfo"][0]]
        )

    instruction["method"], instruction["method2"], = config.get_method_name(
        instruction["jobtype"],
        func=instruction["func"],
        basis=instruction["basis"],
        sm=instruction["sm"],
        solvent=instruction["solvent"],
        prog=config.prog,
        gfn_version=instruction.get("gfn_version", ""),
    )

    name = "efficient gas-phase single-point"
    # folder = "part0_sp"
    check = {True: "was successful", False: "FAILED"}
    if calculate:
        print(f"The {name} is calculated for:")
        print_block(["CONF" + str(i.id) for i in calculate])
        # create folders:
        save_errors, store_confs, calculate = new_folders(
            config.cwd, calculate, folder, save_errors, store_confs
        )
        # write coord to folder
        calculate, store_confs, save_errors = ensemble2coord(
            config, folder, calculate, store_confs, save_errors
        )
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

        for conf in list(calculate):
            if instruction["jobtype"] == "alpb_gsolv":
                line = (
                    f"The {name} {check[conf.job['success']]}"
                    f" for {last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"E(DFT) = {conf.job['energy']:>.8f}"
                    f" Gsolv = {conf.job['energy2']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.cheap_prescreening_sp_info["info"] = "failed"
                    conf.cheap_prescreening_sp_info["method"] = conf.job["method"]
                    conf.cheap_prescreening_gsolv_info["info"] = "failed"
                    conf.cheap_prescreening_gsolv_info["method"] = conf.job["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.cheap_prescreening_sp_info["energy"] = conf.job["energy"]
                    conf.cheap_prescreening_sp_info["info"] = "calculated"
                    conf.cheap_prescreening_sp_info["method"] = conf.job["method"]
                    # conf.cheap_prescreening_gsolv_info["energy"] = conf.job["energy2"]
                    conf.cheap_prescreening_gsolv_info["range"] = {
                        config.fixed_temperature: conf.job["energy2"]
                    }
                    conf.cheap_prescreening_gsolv_info["info"] = "calculated"
                    conf.cheap_prescreening_gsolv_info["method"] = conf.job["method2"]
                    conf.cheap_prescreening_gsolv_info["gas-energy"] = conf.job[
                        "energy_xtb_gas"
                    ]
                    conf.cheap_prescreening_gsolv_info["solv-energy"] = conf.job[
                        "energy_xtb_solv"
                    ]
            elif instruction["jobtype"] == "sp":
                line = (
                    f"The {name} calculation {check[conf.job['success']]}"
                    f" for {last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"E(DFT) = {conf.job['energy']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.cheap_prescreening_sp_info["info"] = "failed"
                    conf.cheap_prescreening_sp_info["method"] = conf.job["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.cheap_prescreening_sp_info["energy"] = conf.job["energy"]
                    conf.cheap_prescreening_sp_info["info"] = "calculated"
                    conf.cheap_prescreening_sp_info["method"] = conf.job["method"]
            else:
                print_errors(
                    f"{'ERROR:':{WARNLEN}}UNEXPECTED BEHAVIOUR: {conf.job['success']} {conf.job['jobtype']}",
                    save_errors,
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

    if prev_calculated:
        # adding conformers calculated before:
        for conf in list(prev_calculated):
            conf.job["workdir"] = os.path.normpath(
                os.path.join(config.cwd, "CONF" + str(conf.id), config.func)
            )
            if instruction["jobtype"] == "alpb_gsolv":
                print(
                    f"The {name} {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"E(DFT) = {conf.cheap_prescreening_sp_info['energy']:>.8f}"
                    f" Gsolv = {conf.cheap_prescreening_gsolv_info['range'].get(config.fixed_temperature, 0.0):>.8f}"
                )
            elif instruction["jobtype"] == "sp":
                print(
                    f"The {name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"E(DFT) = {conf.cheap_prescreening_sp_info['energy']:>.8f}"
                )
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
    for conf in calculate:
        conf.reset_job_info()
    if not calculate:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)

    # create data for possible SI generation:-----------------------------------
    # Energy:
    ensembledata.si["part0"]["Energy"] = instruction["method"]
    if "DOGCP" in instruction["prepinfo"]:
        ensembledata.si["part0"]["Energy"] += " + GCP"
    # Energy_settings are set above!
    # GmRRHO:
    ensembledata.si["part0"]["G_mRRHO"] = "not included"
    # Solvation:
    if config.solvent == "gas":
        ensembledata.si["part0"]["G_solv"] = "gas-phase"
    else:
        ensembledata.si["part0"]["G_solv"] = instruction["method2"]
    # Geometry:
    ensembledata.si["part0"]["Geometry"] = "GFNn-xTB (input geometry)"
    # QM-CODE:
    ensembledata.si["part0"]["main QM code"] = str(config.prog).upper()
    # Threshold:
    ensembledata.si["part0"]["Threshold"] = f"{config.part0_threshold} kcal/mol"
    # END SI generation --------------------------------------------------------

    # ***************************************************************************
    # sorting by E
    # (remove high lying conformers above part0_threshold)
    print("\n" + "".ljust(int(PLENGTH), "-"))
    print(
        "Removing high lying conformers by improved energy description".center(
            int(PLENGTH), " "
        )
    )
    print("".ljust(int(PLENGTH), "-") + "\n")

    if config.solvent != "gas":
        solvation = "cheap_prescreening_gsolv_info"
    else:
        solvation = None
    rrho = None
    energy = "cheap_prescreening_sp_info"
    for conf in calculate:
        conf.calc_free_energy(
            e=energy,
            solv=solvation,
            rrho=rrho,
            t=config.fixed_temperature,
            consider_sym=config.consider_sym,
        )
    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
    except ValueError:
        raise
    for conf in calculate:
        if conf.free_energy == minfree:
            ensembledata.bestconf["part0"] = conf.id
            lowest_e = conf.cheap_prescreening_sp_info["energy"]
            # lowest_gsolv = conf.cheap_prescreening_gsolv_info["energy"]
            lowest_gsolv = (
                getattr(conf, "cheap_prescreening_gsolv_info")
                .get("range", {})
                .get(config.fixed_temperature, 0.0)
            )
    if config.solvent != "gas":
        try:
            minfree_gfnx = min(
                [
                    i.cheap_prescreening_gsolv_info["solv-energy"]
                    for i in calculate
                    if i is not None
                ]
            )
        except ValueError:
            raise
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
        if config.solvent != "gas":
            conf.tmp_rel_xtb = (
                conf.cheap_prescreening_gsolv_info["solv-energy"] - minfree_gfnx
            ) * AU2KCAL
        conf.tmp_rel_e = (
            conf.cheap_prescreening_sp_info["energy"] - lowest_e
        ) * AU2KCAL
        conf.tmp_rel_gsolv = (
            getattr(conf, "cheap_prescreening_gsolv_info")
            .get("range", {})
            .get(config.fixed_temperature, 0.0)
            - lowest_gsolv
        ) * AU2KCAL

    try:
        maxreldft = max([i.rel_free_energy for i in calculate if i is not None])
    except ValueError:
        print_errors(
            f"{'ERROR:':{WARNLEN}}No conformer left or error in maxreldft!", save_errors
        )
    # print sorting
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "cheap_prescreening_gsolv_info")["solv-energy"],
        lambda conf: getattr(conf, "tmp_rel_xtb"),
        lambda conf: getattr(conf, "cheap_prescreening_sp_info")["energy"],
        lambda conf: getattr(conf, "cheap_prescreening_gsolv_info")
        .get("range", {})
        .get(config.fixed_temperature, 0.0),
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "tmp_rel_e"),
        lambda conf: getattr(conf, "tmp_rel_gsolv"),
        lambda conf: getattr(conf, "rel_free_energy"),
    ]
    columnheader = [
        "CONF#",
        f"E [Eh]",
        f"ΔE [kcal/mol]",
        "E [Eh]",
        "Gsolv [Eh]",
        "gtot",
        "ΔE(DFT)",
        "ΔGsolv",
        "Δgtot",
    ]
    columndescription = [
        "",
        "",
        "",
        "",
        "[Eh]",
        "[Eh]",
        "[kcal/mol]",
        "[kcal/mol]",
        "[kcal/mol]",
    ]
    columnformat = [
        "",
        (12, 7),
        (5, 2),
        (12, 7),
        (12, 7),
        (12, 7),
        (5, 2),
        (5, 2),
        (5, 2),
    ]
    columndescription[1] = f"{config.part0_gfnv.upper()}-xTB[{config.sm_rrho}]"
    columndescription[2] = f"{config.part0_gfnv.upper()}-xTB[{config.sm_rrho}]"
    columndescription[3] = instruction["method"]
    columndescription[4] = instruction["method2"]
    if config.solvent == "gas":
        columncall[1] = lambda conf: getattr(conf, "xtb_energy")
        columncall[2] = lambda conf: getattr(conf, "rel_xtb_energy")
        columnheader[1] = "G(GFNn-xTB)"
        columnheader[2] = "ΔG(GFNn-xTB)"
        columndescription[1] = "[Eh]"
        columndescription[2] = "[kcal/mol]"
        columndescription[4] = "gas-phase"

    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part0.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
    )
    print("".ljust(int(PLENGTH), "-"))
    # --------------------------------------------------------------------------
    calculate = calc_boltzmannweights(
        calculate, "free_energy", config.fixed_temperature
    )
    conf_in_interval(calculate, full_free_energy=False)
    # --------------------------------------------------------------------------

    # write to enso.json
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )
    # ---------------------------------------------------------------------------
    # sorting
    if maxreldft > config.part0_threshold:
        print("\n" + "".ljust(int(PLENGTH / 2), "-"))
        print("Conformers considered further".center(int(PLENGTH / 2), " "))
        print("".ljust(int(PLENGTH / 2), "-") + "\n")
        for conf in list(calculate):
            if conf.rel_free_energy <= config.part0_threshold:
                conf.part_info["part0"] = "passed"
            else:
                conf.part_info["part0"] = "refused"
                store_confs.append(calculate.pop(calculate.index(conf)))
        if calculate:
            print(
                f"These conformers are below the {config.part0_threshold:.3f} "
                f"kcal/mol g_thr(0) threshold.\n"
            )
            print_block(["CONF" + str(i.id) for i in calculate])
        else:
            print_errors(
                f"{'ERROR:':{WARNLEN}}There are no more conformers left!", save_errors
            )
    else:
        for conf in list(calculate):
            conf.part_info["part0"] = "passed"
        print(
            "\nAll relative (free) energies are below the initial g_thr(0) threshold "
            f"of {config.part0_threshold} kcal/mol.\nAll conformers are "
            "considered further."
        )
    ensembledata.nconfs_per_part["part0"] = len(calculate)

    ################################################################################
    # calculate average G correction
    print(
        "\nCalculating Boltzmann averaged (free) energy of ensemble on input "
        "geometries (not DFT optimized)!\n"
    )
    # calculate Boltzmannweights
    print(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} " f"{'avG(T) /a.u.':>14} ")
    print("".ljust(int(PLENGTH), "-"))

    calculate = calc_boltzmannweights(
        calculate, "free_energy", config.fixed_temperature
    )
    avG = 0.0
    avE = 0.0
    for conf in calculate:
        avG += conf.bm_weight * conf.free_energy
        avE += conf.bm_weight * conf.cheap_prescreening_sp_info["energy"]
    # printout:
    print(f"{config.fixed_temperature:^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0==")
    print("".ljust(int(PLENGTH), "-"))
    print("")
    ################################################################################

    # reset
    for conf in calculate:
        conf.free_energy = 0.0
        conf.rel_free_energy = None
        conf.bm_weight = 0.0
        conf.tmp_rel_xtb = 0.0
        conf.tmp_rel_e = 0.0
        conf.tmp_rel_gsolv = 0.0
        conf.reset_job_info()

    # write to enso.json
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    if save_errors:
        print("\n***---------------------------------------------------------***")
        print("Printing most relevant errors again, just for user convenience:")
        for _ in list(save_errors):
            print(save_errors.pop())
        print("***---------------------------------------------------------***")

    tmp = int((PLENGTH - len("END of Part0")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part0" + "".rjust(tmp, "<"))
    if config.progress:
        print("#>>># CENSO: Finished part0", file=sys.stderr)
    return config, calculate, store_confs, ensembledata
