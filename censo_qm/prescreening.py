"""
prescreening == part1, calculate free energy on GFNn-xTB input geometry
idea is to improve on E and (Gsolv)
"""
import os
import sys
import math
from multiprocessing import JoinableQueue as Queue
from .cfg import PLENGTH, DIGILEN, AU2KCAL, CODING, censo_solvent_db
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
    move_recursively,
    write_trj,
    check_tasks,
    calc_std_dev,
    spearman,
    print,
    print_errors,
    calc_boltzmannweights,
)


def part1(config, conformers, store_confs, ensembledata):
    """
    Prescreening of the ensemble, with single-points on combined ensemble
    geometries.
    Calculate low level free energies with COSMO-RS single-point and gsolv
    contribution and GFNFF-bhess thermostatistical contribution.
    Input:
    - config [conifg_setup object] contains all settings
    - conformers [list of molecule_data objects] each conformer is represented
    - ensembledata -> instance for saving ensemble (not conf) related data
    Return:
    -> config
    -> conformers
    -> store_confs
    """
    save_errors = []
    print("\n" + "".ljust(PLENGTH, "-"))
    print("CRE PRESCREENING - PART1".center(PLENGTH, " "))
    print("".ljust(PLENGTH, "-") + "\n")
    # print flags for part1
    info = []
    info.append(["prog", "program"])
    info.append(["func", "functional for part1 and 2"])
    info.append(["basis", "basis set for part1 and 2"])
    if config.solvent != "gas":
        info.append(["solvent", "Solvent"])
        info.append(["smgsolv1", "solvent model for Gsolv contribution"])
    info.append(["part1_threshold", "threshold"])
    info.append(
        ["printoption", "starting number of considered conformers", len(conformers)]
    )
    info.append(["evaluate_rrho", "calculate mRRHO contribution"])
    if config.evaluate_rrho:
        info.append(["prog_rrho", "program for mRRHO contribution"])
        if config.prog_rrho == "xtb":
            info.append(["part1_gfnv", "GFN version for mRRHO and/or GBSA_Gsolv"])
            if config.bhess:
                info.append(
                    [
                        "bhess",
                        "Apply constraint to input geometry during mRRHO calculation",
                    ]
                )
    info.append(["temperature", "temperature"])
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
                    item[1], "", option, digits=DIGILEN - len(item[1])
                )
            )
    print("")
    # end print

    calculate = []  # has to be calculated in this run
    prev_calculated = []  # was already calculated in a previous run
    try:
        store_confs
    except NameError:
        store_confs = []  # stores all confs which are sorted out!

    if config.solvent == "gas":
        print("Calculating single-point energies:")
    else:
        print("Calculating single-point energies and solvation contribution (G_solv):")

    # setup queues
    q = Queue()
    resultq = Queue()

    if config.prog == "tm":
        job = TmJob
    elif config.prog == "orca":
        job = OrcaJob

    for conf in list(conformers):
        if conf.removed:
            store_confs.append(conformers.pop(conformers.index(conf)))
            print(f"CONF{conf.id} is removed as requested by the user.")
            continue
        if conf.id > config.nconf:
            store_confs.append(conformers.pop(conformers.index(conf)))
            continue
        if conf.prescreening_sp_info["info"] == "not_calculated":
            conf = conformers.pop(conformers.index(conf))
            calculate.append(conf)
        elif conf.prescreening_sp_info["info"] == "failed":
            conf = conformers.pop(conformers.index(conf))
            store_confs.append(conf)
            print(f"Calculation of CONF{conf.id} failed in the previous run.")
        elif conf.prescreening_sp_info["info"] == "calculated":
            conf = conformers.pop(conformers.index(conf))
            if config.solvent != "gas":
                # check if solvation calculation calculated as well!
                if conf.prescreening_gsolv_info["info"] == "failed":
                    store_confs.append(conf)
                    print(
                        f"Calculation of the solvation contribution for CONF"
                        f"{conf.id} failed in the previous run."
                    )
                elif (
                    conf.prescreening_gsolv_info["info"] == "not_calculated"
                    and config.smgsolv1 in config.smgsolv_1
                ):
                    # additive solvation
                    calculate.append(conf)
                else:
                    # implicit solvation
                    conf.job["success"] = True
                    prev_calculated.append(conf)
            elif config.solvent == "gas":
                conf.job["success"] = True
                prev_calculated.append(conf)
    if not calculate and not prev_calculated:
        print("ERROR: No conformers left!")
    if prev_calculated:
        check_for_folder(config.cwd, [i.id for i in prev_calculated], config.func)
        print("The prescreening_single-point was calculated before for:")
        print_block(["CONF" + str(i.id) for i in prev_calculated])
    pl = config.lenconfx + 4 + len(str("/" + config.func))

    instruction = {
        "prepinfo": ["low+"],  # TM: m4 scfconv 6
        "func": config.func,
        "basis": getattr(
            config, "basis", config.func_basis_default.get(config.func, "def2-mTZVPP")
        ),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": config.smgsolv1,
        "omp": config.omp,
        "temperature": config.temperature,
        "gfn_version": config.part1_gfnv,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
    }

    if config.solvent == "gas":
        instruction["jobtype"] = "sp"
        instruction["method"], _ = config.get_method_name(
            "sp", func=instruction["func"], basis=instruction["basis"]
        )
        name = "prescreening_single-point"
        folder = instruction["func"]
        if config.prog == "orca":
            instruction["progpath"] = config.external_paths["orcapath"]
    else:
        if config.smgsolv1 in config.smgsolv_1:
            # additive Gsolv
            # COSMORS
            if config.smgsolv1 != "dcosmors" and "cosmors" in config.smgsolv1:
                job = TmJob
                exc_fine = {"cosmors": "normal", "cosmors-fine": "fine"}
                tmp = {
                    "jobtype": "cosmors",
                    "cosmorssetup": config.external_paths["cosmorssetup"],
                    "cosmorsparam": exc_fine.get(config.smgsolv1, "normal"),
                    "cosmothermversion": config.external_paths["cosmothermversion"],
                    "ctd-param": config.cosmorsparam,
                }
                instruction.update(tmp)
                instruction["method"], instruction["method2"] = config.get_method_name(
                    "cosmors",
                    func=instruction["func"],
                    basis=instruction["basis"],
                    sm=instruction["sm"],
                )
                name = "prescreening COSMO-RS"
                folder = str(instruction["func"]) + "/COSMO"
            # GBSA-Gsolv / ALPB-Gsolv
            elif instruction["sm"] in ("gbsa_gsolv", "alpb_gsolv"):
                # do DFT gas phase sp and additive Gsolv
                instruction["jobtype"] = instruction["sm"]
                if config.prog == "orca":
                    instruction["progpath"] = config.external_paths["orcapath"]
                instruction["xtb_driver_path"] = config.external_paths["xtbpath"]
                instruction["method"], instruction["method2"] = config.get_method_name(
                    instruction["jobtype"],
                    func=instruction["func"],
                    basis=instruction["basis"],
                    sm=instruction["sm"],
                    gfn_version=instruction["gfn_version"],
                )
                if (
                    conf.prescreening_sp_info["info"] == "calculated"
                    and conf.prescreening_sp_info["method"] == instruction["method"]
                ):
                    # do not calculate gas phase sp again!
                    instruction["energy"] = conf.prescreening_sp_info["energy"]
                    instruction["prepinfo"] = []
                name = "prescreening " + str(instruction["sm"]).upper()
                folder = str(instruction["func"]) + "/Gsolv"
            # SMD-Gsolv
            elif instruction["sm"] == "smd_gsolv":
                job = OrcaJob
                instruction["jobtype"] = instruction["sm"]
                instruction["progpath"] = config.external_paths["orcapath"]
                instruction["method"], instruction["method2"] = config.get_method_name(
                    "smd_gsolv",
                    func=instruction["func"],
                    basis=instruction["basis"],
                    sm=instruction["sm"],
                )
                name = "prescreening" + str(instruction["sm"]).upper()
                folder = str(instruction["func"]) + "/Gsolv"
        else:
            # with implicit solvation
            instruction["jobtype"] = "sp_implicit"
            # instruction["prepinfo"] = ["low+"]
            if config.prog == "orca":
                instruction["progpath"] = config.external_paths["orcapath"]
            instruction["method"], instruction["method2"] = config.get_method_name(
                "sp_implicit",
                func=instruction["func"],
                basis=instruction["basis"],
                sm=instruction["sm"],
            )
            name = "prescreening_single-point"
            folder = instruction["func"]

    check = {True: "was successful", False: "FAILED"}
    if calculate:
        print(f"The {name} is calculated for:")
        print_block(["CONF" + str(i.id) for i in calculate])
        # create folders:
        save_errors, store_confs, calculate = new_folders(
            config.cwd, calculate, config.func, save_errors, store_confs
        )
        # write coord to folder
        calculate, store_confs, save_errors = ensemble2coord(
            config, config.func, calculate, store_confs, save_errors
        )
        if config.solvent != "gas":
            if folder != str(config.func):
                # create the COSMO folder
                save_errors, store_confs, calculate = new_folders(
                    config.cwd, calculate, folder, save_errors, store_confs
                )
                # write the coord file to the COSMO folder
                calculate, store_confs, save_errors = ensemble2coord(
                    config, folder, calculate, store_confs, save_errors
                )

        # parallel calculation:
        calculate = run_in_parallel(
            config, q, resultq, job, config.maxthreads, calculate, instruction, folder
        )

        for conf in list(calculate):
            if instruction["jobtype"] == "sp":
                line = (
                    f"{name} calculation {check[conf.job['success']]}"
                    f" for {last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.prescreening_sp_info["info"] = "failed"
                    conf.prescreening_sp_info["method"] = conf.job["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.prescreening_sp_info["energy"] = conf.job["energy"]
                    conf.prescreening_sp_info["info"] = "calculated"
                    conf.prescreening_sp_info["method"] = conf.job["method"]
            elif instruction["jobtype"] == "sp_implicit":
                line = (
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.prescreening_sp_info["info"] = "failed"
                    conf.prescreening_sp_info["method"] = conf.job["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.prescreening_sp_info["energy"] = conf.job["energy"]
                    conf.prescreening_sp_info["info"] = "calculated"
                    conf.prescreening_sp_info["method"] = conf.job["method"]
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
                    conf.prescreening_sp_info["info"] = "failed"
                    conf.prescreening_sp_info["method"] = conf.job["method"]
                    conf.prescreening_gsolv_info["info"] = "failed"
                    conf.prescreening_gsolv_info["method"] = conf.job["method2"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.prescreening_sp_info["energy"] = conf.job["energy"]
                    conf.prescreening_sp_info["info"] = "calculated"
                    conf.prescreening_sp_info["method"] = conf.job["method"]
                    conf.prescreening_gsolv_info["gas-energy"] = conf.job["energy"]
                    conf.prescreening_gsolv_info["energy"] = conf.job["energy2"]
                    conf.prescreening_gsolv_info["info"] = "calculated"
                    conf.prescreening_gsolv_info["method"] = conf.job["method2"]
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
    if prev_calculated:
        # adding conformers calculated before:
        for conf in list(prev_calculated):
            conf.job["workdir"] = os.path.normpath(
                os.path.join(config.cwd, "CONF" + str(conf.id), config.func)
            )
            if instruction["jobtype"] in ("sp", "sp_implicit"):
                print(
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.prescreening_sp_info['energy']:>.8f}"
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
                    f"{conf.prescreening_gsolv_info['energy']:>.8f}"
                )
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
    for conf in calculate:
        conf.reset_job_info()
    if not calculate:
        print_errors("ERROR: No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)

    # ***************************************************************************
    # first sorting by E or Gsolv
    # (remove high lying conformers above part1_threshold + 1.5 kcal/mol)
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("Removing high lying conformers".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")

    for conf in calculate:
        rrho = None
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        e = "prescreening_sp_info"
        conf.calc_free_energy(e=e, solv=solv, rrho=rrho)
    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
    except ValueError:
        raise
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
    try:
        maxreldft = max([i.rel_free_energy for i in calculate if i is not None])
    except ValueError:
        print_errors("ERROR: No conformer left or Error in maxreldft!", save_errors)
    # print sorting
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "xtb_energy"),
        lambda conf: getattr(conf, "rel_xtb_energy"),
        lambda conf: getattr(conf, "prescreening_sp_info")["energy"],
        lambda conf: getattr(conf, "prescreening_gsolv_info")["energy"],
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "rel_free_energy"),
    ]
    columnheader = [
        "CONF#",
        "E(GFNn-xTB)",
        "ΔE(GFNn-xTB)",
        "E [Eh]",
        "Gsolv [Eh]",
        "Gtot",
        "ΔGtot",
    ]
    columndescription = ["", "[a.u.]", "[kcal/mol]", "", "", "[Eh]", "[kcal/mol]"]
    columndescription2 = ["", "", "", "", "", "", "", ""]
    columnformat = ["", (12, 7), (5, 2), (12, 7), (12, 7), (12, 7), (5, 2)]

    if config.solvent == "gas":
        columnheader[5] = "Etot"
        columnheader[6] = "ΔEtot"
        columndescription[3] = instruction["method"]
        # ignore gsolv in printout
        columncall.pop(4)
        columnheader.pop(4)
        columndescription.pop(4)
        columnformat.pop(4)
    elif config.solvent != "gas":
        # energy
        columndescription[3] = instruction["method"]
        # gsolv
        columndescription[4] = instruction["method2"]

    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part1preG.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
        columndescription2=columndescription2,
    )

    if maxreldft > (config.part1_threshold):
        print("\n" + "".ljust(int(PLENGTH / 2), "-"))
        print("Conformers considered further".center(int(PLENGTH / 2), " "))
        print("".ljust(int(PLENGTH / 2), "-") + "\n")
        for conf in list(calculate):
            if conf.rel_free_energy > (config.part1_threshold):
                store_confs.append(calculate.pop(calculate.index(conf)))
        if calculate:
            print(f"Below the threshold of {config.part1_threshold} kcal/mol.\n")
            print_block(["CONF" + str(i.id) for i in calculate])
        else:
            print("Error: There are no more conformers left!")
    else:
        print(
            "\nAll relative (free) energies are below the threshold "
            f"of ({config.part1_threshold} kcal/mol.\nAll conformers are "
            "considered further."
        )
    ensembledata.nconfs_per_part["part1_firstsort"] = len(calculate)
    # reset
    for conf in calculate:
        conf.free_energy = 0.0
        conf.rel_free_energy = None
    print("".ljust(int(PLENGTH / 2), "-"))
    # ***************************************************************************
    if config.evaluate_rrho:
        # check if prescreening rrho has been calculated
        if config.solvent == "gas":
            print("\nCalculating prescreening G_mRRHO!")
        else:
            print("\nCalculating prescreening G_mRRHO with implicit solvation!")

        for conf in list(calculate):
            if conf.prescreening_grrho_info["info"] == "not_calculated":
                pass
            elif conf.prescreening_grrho_info["info"] == "failed":
                conf = calculate.pop(calculate.index(conf))
                conf.__class__ = job
                store_confs.append(conf)
                print(f"Calculation of CONF{conf.id} failed in the previous run!")
            elif conf.prescreening_grrho_info["info"] == "calculated":
                conf = calculate.pop(calculate.index(conf))
                conf.__class__ = job
                conf.job["success"] = True
                prev_calculated.append(conf)

        if not calculate and not prev_calculated:
            print_errors("ERROR: No conformers left!", save_errors)
            print("Going to exit!")
            sys.exit(1)
        folderrho = "rrho_part1"
        if prev_calculated:
            check_for_folder(config.cwd, [i.id for i in prev_calculated], folderrho)
            print("The prescreening G_mRRHO calculation was performed before for:")
            print_block(["CONF" + str(i.id) for i in prev_calculated])
        pl = config.lenconfx + 4 + len(str("/" + folderrho))
        instruction_prerrho = {
            "jobtype": "rrhoxtb",
            "func": getattr(config, "part1_gfnv"),
            "gfn_version": getattr(config, "part1_gfnv"),
            "temperature": config.temperature,
            "charge": config.charge,
            "unpaired": config.unpaired,
            "solvent": config.solvent,
            "omp": config.omp,
            "progpath": config.external_paths["xtbpath"],
            "bhess": config.bhess,
            "consider_sym": config.consider_sym,
            "sm_rrho": config.sm_rrho,
            "rmsdbias": config.rmsdbias,
            "cwd": config.cwd,
            "copymos": "",
            "sym": "c1",
            "multiTemp": False,
            "energy": 0.0,
            "energy2": 0.0,
            "success": False,
            "imagthr": config.imagthr,
            "sthr": config.sthr,
            "scale":config.scale,
        }

        instruction_prerrho["method"], _ = config.get_method_name(
            "rrhoxtb",
            bhess=config.bhess,
            gfn_version=instruction_prerrho["gfn_version"],
            sm=instruction_prerrho["sm_rrho"],
            solvent=instruction_prerrho["solvent"],
        )
        if calculate:
            print("The prescreening G_mRRHO calculation is now performed for:")
            print_block(["CONF" + str(i.id) for i in calculate])
            # create folders:
            save_errors, store_confs, calculate = new_folders(
                config.cwd, calculate, folderrho, save_errors, store_confs
            )
            # write coord to folder
            calculate, store_confs, save_errors = ensemble2coord(
                config, folderrho, calculate, store_confs, save_errors
            )
            calculate = run_in_parallel(
                config,
                q,
                resultq,
                job,
                config.maxthreads,
                calculate,
                instruction_prerrho,
                folderrho,
            )
            check = {True: "was successful", False: "FAILED"}
            # check if too many calculations failed

            ###
            for conf in list(calculate):
                print(
                    f"The prescreening G_mRRHO calculation @ {conf.job['symmetry']} "
                    f"{check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f}"
                )
                if not conf.job["success"]:
                    conf.prescreening_grrho_info["info"] = "failed"
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.sym = conf.job["symmetry"]
                    conf.prescreening_grrho_info["rmsd"] = conf.job["rmsd"]
                    conf.prescreening_grrho_info["energy"] = conf.job["energy"]
                    conf.prescreening_grrho_info["info"] = "calculated"
                    conf.prescreening_grrho_info["method"] = instruction_prerrho[
                        "method"
                    ]

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
            for conf in list(prev_calculated):
                conf.job["workdir"] = os.path.normpath(
                    os.path.join(config.cwd, "CONF" + str(conf.id), folderrho)
                )
                print(
                    f"The prescreening G_mRRHO calculation @ {conf.sym} "
                    f"{check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.prescreening_grrho_info['energy']:>.8f}"
                )
                calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
        if not calculate:
            print_errors("ERROR: No conformers left!", save_errors)
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

    # # printout for part1 -------------------------------------------------------
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("* Gibbs free energies of part1 *".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "xtb_free_energy"),
        lambda conf: getattr(conf, "rel_xtb_free_energy"),
        lambda conf: getattr(conf, "prescreening_sp_info")["energy"],
        lambda conf: getattr(conf, "prescreening_gsolv_info")["energy"],
        lambda conf: getattr(conf, "prescreening_grrho_info")["energy"],
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "rel_free_energy"),
    ]
    columnheader = [
        "CONF#",
        "G(GFNn-xTB)",
        "ΔG(GFNn-xTB)",
        "E [Eh]",
        "Gsolv [Eh]",
        "GmRRHO [Eh]",
        "Gtot",
        "ΔGtot",
    ]
    columndescription = [
        "",  # CONFX
        "[a.u.]",  # xtb energy
        "[kcal/mol]",  # rel xtb_energy
        str(config.func),  # E
        "",  # GSolv
        "",
        "[Eh]",  # Gtot
        "[kcal/mol]",  # rel Gtot
    ]
    columndescription2 = ["", "", "", "", "", "", "", ""]
    columnformat = ["", (12, 7), (5, 2), (12, 7), (12, 7), (12, 7), (12, 7), (5, 2)]
    if config.solvent == "gas":
        # Energy
        columndescription[3] = instruction["method"]
    elif config.solvent != "gas":
        # Energy
        columndescription[3] = instruction["method"]
        # Gsolv
        columndescription[4] = instruction["method2"]
    if config.evaluate_rrho:
        columndescription[5] = instruction_prerrho["method"]  # Grrho

    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho:
            # ignore rrho in printout
            columncall.pop(5)
            columnheader.pop(5)
            columndescription.pop(5)
            columnformat.pop(5)
        if config.solvent == "gas":
            columncall.pop(4)
            columnheader.pop(4)
            columndescription.pop(4)
            columnformat.pop(4)

    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "prescreening_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        conf.calc_free_energy(e="prescreening_sp_info", solv=solv, rrho=rrho)
        conf.xtb_free_energy = conf.calc_free_energy(
            e="xtb_energy", solv=None, rrho=rrho, out=True
        )

    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
        minfree_xtb = min([i.xtb_free_energy for i in calculate if i is not None])
    except ValueError:
        raise ValueError
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
        conf.rel_xtb_free_energy = (conf.xtb_free_energy - minfree_xtb) * AU2KCAL
    try:
        maxreldft = max([i.rel_free_energy for i in calculate if i is not None])
    except ValueError:
        print_errors("ERROR: No conformer left or Error in maxreldft!", save_errors)
    # print sorting
    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part1.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
        columndescription2=columndescription2,
    )
    # --------------------------------------------------------------------------
    for conf in calculate:
        if conf.free_energy == minfree:
            ensembledata.bestconf["part1"] = conf.id

    # write to enso.json
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    # ***************************************************************************
    # Fuzzy or smart sorting
    # increase the individual threshold for conformers with GRRHO differing from
    # the mean GmRRHO
    if len(calculate) == 1:
        std_dev = 0.0
    else:
        std_dev = calc_std_dev(
            [
                conf.prescreening_grrho_info["energy"] * AU2KCAL
                for conf in calculate
                if conf.prescreening_grrho_info["energy"] is not None
            ]
        )
    max_fuzzy = 1
    fuzzythr = max_fuzzy * (1 - math.exp(-1 * 5 * (std_dev ** 2)))
    print(
        "\nAdditional global 'fuzzy-threshold' based on the standard deviation of (G_mRRHO):"
    )
    print(f"Std_dev(G_mRRHO) = {std_dev:.3f} kcal/mol")
    print(f"Fuzzythreshold   = {fuzzythr:.3f} kcal/mol")
    print(
        f"Final sorting threshold = {config.part1_threshold:.3f} + "
        f"{fuzzythr:.3f} = {config.part1_threshold + fuzzythr:.3f} kcal/mol"
    )
    for conf in calculate:
        conf.prescreening_grrho_info["fuzzythr"] = fuzzythr

    # spearman between DFT and DFT + RRHO
    if config.evaluate_rrho and len(calculate) > 1:
        for conf in calculate:
            rrho = None
            if config.solvent == "gas":
                solv = None
            else:
                solv = "prescreening_gsolv_info"
            e = "prescreening_sp_info"
            conf.calc_free_energy(e=e, solv=solv, rrho=rrho)
        try:
            minfree = min([i.free_energy for i in calculate if i is not None])
        except ValueError:
            raise ValueError
        without_RRHO = []
        calculate.sort(key=lambda x: int(x.id))
        for conf in calculate:
            without_RRHO.append((conf.free_energy - minfree) * AU2KCAL)
        for conf in calculate:
            conf.free_energy = 0.0
        for conf in calculate:
            rrho = "prescreening_grrho_info"
            if config.solvent == "gas":
                solv = None
            else:
                solv = "prescreening_gsolv_info"
            e = "prescreening_sp_info"
            conf.calc_free_energy(e=e, solv=solv, rrho=rrho)
        try:
            minfree = min([i.free_energy for i in calculate if i is not None])
        except ValueError:
            raise ValueError
        with_RRHO = []
        calculate.sort(key=lambda x: int(x.id))
        for conf in calculate:
            with_RRHO.append((conf.free_energy - minfree) * AU2KCAL)
        for conf in calculate:
            conf.free_energy = 0.0
        if config.solvent != "gas":
            print(
                f"Spearman correlation coefficient between (E + Solv) "
                f"and (E + Solv + mRRHO) = {spearman(without_RRHO, with_RRHO):.3f}"
            )
        else:
            print(
                f"Spearman correlation coefficient between (E) "
                f"and (E + mRRHO) = {spearman(without_RRHO, with_RRHO):.3f}"
            )

    # sorting
    if maxreldft > config.part1_threshold:
        print("\n" + "".ljust(int(PLENGTH / 2), "-"))
        print("Conformers considered further".center(int(PLENGTH / 2), " "))
        print("".ljust(int(PLENGTH / 2), "-") + "\n")
        for conf in list(calculate):
            if conf.rel_free_energy <= config.part1_threshold:
                conf.part_info["part1"] = "passed"
            elif conf.rel_free_energy <= (
                config.part1_threshold + conf.prescreening_grrho_info["fuzzythr"]
            ):
                print(f"Considered CONF{conf.id} because of increased fuzzythr.")
                conf.part_info["part1"] = "passed"
                continue
            else:
                conf.part_info["part1"] = "refused"
                store_confs.append(calculate.pop(calculate.index(conf)))
        if calculate:
            print(
                f"These conformers are below the {config.part1_threshold+fuzzythr:.3f} "
                f"kcal/mol threshold.\n"
            )
            print_block(["CONF" + str(i.id) for i in calculate])
        else:
            print_errors("Error: There are no more conformers left!", save_errors)
    else:
        for conf in list(calculate):
            conf.part_info["part1"] = "passed"
        print(
            "\nAll relative (free) energies are below the initial threshold "
            f"of {config.part1_threshold} kcal/mol.\nAll conformers are "
            "considered further."
        )
    ensembledata.nconfs_per_part["part1"] = len(calculate)
    # save current data to jsonfile
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    # free energy:
    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "prescreening_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        conf.calc_free_energy(e="prescreening_sp_info", solv=solv, rrho=rrho)

    # write coord.enso_best
    for conf in calculate:
        if conf.id == ensembledata.bestconf["part1"]:
            # copy the lowest optimized conformer to file coord.enso_best
            with open(
                os.path.join("CONF" + str(conf.id), config.func, "coord"),
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
                        conf.prescreening_grrho_info["energy"],
                        conf.id,
                    )
                )
                for line in coord[1:]:
                    if "$" in line:  # stop at $end ...
                        break
                    best.write(line)
                best.write("$end \n")

    ################################################################################
    # calculate average G correction
    print(
        "\nCalculating Boltzmann averaged free energy of ensemble on "
        f"input geometries (not DFT optimized)!\n"
    )
    # calculate Boltzmannweights
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho and config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avG(T) /a.u.':>14} "
            )
        elif not config.evaluate_rrho:
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGsolv(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
            )
        elif config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGmRRHO(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
            )
    else:
        line = (
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
            f"{'avGmRRHO(T) /a.u.':>16} {'avGsolv(T) /a.u.':>16} "
            f"{'avG(T) /a.u.':>14}"
        )
    print(line)
    print("".ljust(int(PLENGTH), "-"))
    # get free energy at (T)
    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "prescreening_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        e = "prescreening_sp_info"
        conf.calc_free_energy(e=e, solv=solv, rrho=rrho)

    calculate = calc_boltzmannweights(calculate, "free_energy", config.temperature)
    avG = 0.0
    avE = 0.0
    avGRRHO = 0.0
    avGsolv = 0.0
    for conf in calculate:
        avG += conf.bm_weight * conf.free_energy
        avE += conf.bm_weight * conf.prescreening_sp_info["energy"]
        avGRRHO += conf.bm_weight * conf.prescreening_grrho_info["energy"]
        avGsolv += conf.bm_weight * conf.prescreening_gsolv_info["energy"]

    # printout:
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho and config.solvent == "gas":
            line = f"{config.temperature:^15} {avE:>14.7f}  {avG:>14.7f} "
        elif not config.evaluate_rrho:
            line = (
                f"{config.temperature:^15} {avE:>14.7f} {avGsolv:>16.7f} "
                f"{avG:>14.7f} "
            )
        elif config.solvent == "gas":
            line = (
                f"{config.temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                f"{avG:>14.7f} "
            )
    else:
        line = (
            f"{config.temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
            f"{avGsolv:>16.7f} {avG:>14.7f} "
        )
    print(line, "    <<==part1==")
    print("".ljust(int(PLENGTH), "-"))
    print("")

    ################################################################################

    print("\nCalculating unbiased GFNn-xTB energy")
    instruction_gfn = {
        "jobtype": "xtb_sp",
        "func": getattr(config, "part1_gfnv"),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "progpath": config.external_paths["xtbpath"],
        "sm": config.sm_rrho,
        "rmsdbias": config.rmsdbias,
        "omp": config.omp,
        "temperature": config.temperature,
        "gfn_version": config.part1_gfnv,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
    }
    if calculate:
        folder_gfn = "GFN_unbiased"
        save_errors, store_confs, calculate = new_folders(
            config.cwd, calculate, folder_gfn, save_errors, store_confs
        )
        # write coord to folder
        calculate, store_confs, save_errors = ensemble2coord(
            config, folder_gfn, calculate, store_confs, save_errors
        )
        calculate = run_in_parallel(
            config,
            q,
            resultq,
            job,
            config.maxthreads,
            calculate,
            instruction_gfn,
            folder_gfn,
        )
        for conf in list(calculate):
            if not conf.job["success"]:
                conf.xtb_energy_unbiased = conf.xtb_energy
            else:
                conf.xtb_energy_unbiased = conf.job["energy"]
    # write ensemble
    move_recursively(config.cwd, "enso_ensemble_part1.xyz")
    if config.evaluate_rrho:
        kwargs = {"energy": "xtb_energy_unbiased", "rrho": "prescreening_grrho_info"}
    else:
        kwargs = {"energy": "xtb_energy_unbiased"}
    write_trj(
        sorted(calculate, key=lambda x: float(x.free_energy)),
        config.cwd,
        "enso_ensemble_part1.xyz",
        config.func,
        config.nat,
        "free_energy",
        **kwargs,
    )

    # reset
    for conf in calculate:
        conf.free_energy = 0.0
        conf.rel_free_energy = None
        conf.bm_weight = 0.0
        conf.reset_job_info()
    if save_errors:
        print(
            "\n***---------------------------------------------------------***",
            file=sys.stderr,
        )
        print(
            "Printing most relevant errors again, just for user convenience:",
            file=sys.stderr,
        )
        for _ in list(save_errors):
            print(save_errors.pop(), file=sys.stderr)
        print(
            "***---------------------------------------------------------***",
            file=sys.stderr,
        )

    tmp = int((PLENGTH - len("END of Part1")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part1" + "".rjust(tmp, "<"))
    return config, calculate, store_confs, ensembledata
