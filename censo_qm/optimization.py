"""
Optimization == part2
performing optimization of the CRE and provide low level free energies.
"""
from multiprocessing import JoinableQueue as Queue
import shutil
import time
import os
import sys
from copy import deepcopy
from .cfg import PLENGTH, CODING, AU2KCAL, DIGILEN, WARNLEN
from .utilities import (
    check_for_folder,
    print_block,
    new_folders,
    last_folders,
    ensemble2coord,
    frange,
    calc_boltzmannweights,
    spearman,
    printout,
    move_recursively,
    write_trj,
    crest_routine,
    check_tasks,
    print,
    print_errors,
    calc_weighted_std_dev,
)
from .orca_job import OrcaJob
from .tm_job import TmJob
from .parallel import run_in_parallel


def part2(config, conformers, store_confs, ensembledata):
    """
    Optimization of the ensemble, at DFT level (possibly with implicit solvation)
    Calculate low level free energies with COSMO-RS single-point and gsolv 
    contribution and GFNFF-bhess thermostatistical contribution on DFT optimized
    geometries
    Input:
    - config [conifg_setup object] contains all settings 
    - conformers [list of molecule_data objects] each conformer is represented
    Return:
    -> config
    -> conformers
    """
    save_errors = []
    if config.progress:
        print("#>>># CENSO: Starting part2", file=sys.stderr)
    print("\n" + "".ljust(PLENGTH, "-"))
    print("CRE OPTIMIZATION - PART2".center(PLENGTH, " "))
    print("".ljust(PLENGTH, "-") + "\n")
    # print flags for part2
    info = []
    info.append(["prog", "program"])
    info.append(["func", "functional for part2"])
    info.append(["basis", "basis set for part2"])
    info.append(["ancopt", "using the xTB-optimizer for optimization"])
    if config.opt_spearman:
        info.append(["opt_spearman", "using the new ensemble optimizer"])
        info.append(
            ["opt_limit", "optimize all conformers below this G_thr(opt,2) threshold"]
        )
        info.append(["printoption", "Spearman threshold", f"{config.spearmanthr:.3f}"])
        info.append(["optcycles", "number of optimization iterations"])
        if config.func == "r2scan-3c":
            info.append(["radsize", "radsize"])
    if config.ancopt and config.optlevel2 is not None:
        info.append(["optlevel2", "optimization level in part2"])
    if config.solvent != "gas":
        info.append(["solvent", "solvent"])
        info.append(["sm2", "solvent model applied in the optimization"])
        if config.smgsolv2 not in (None, "sm"):
            info.append(["smgsolv2", "solvent model for Gsolv contribution"])
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
    info.append(["part2_threshold", "Boltzmann sum threshold G_thr(2) for sorting in part2"])
    info.append(["evaluate_rrho", "calculate mRRHO contribution"])
    if config.evaluate_rrho:
        info.append(["prog_rrho", "program for mRRHO contribution"])
        if config.prog_rrho == "xtb":
            info.append(["part2_gfnv", "GFN version for mRRHO and/or GBSA_Gsolv"])
            if config.bhess:
                info.append(
                    [
                        "bhess",
                        "Apply constraint to input geometry during mRRHO calculation",
                    ]
                )
    max_len_digilen = 0
    for item in info:
        if item[0] == 'justprint':
            if "short-notation" in item[1]:
                tmp = len(item[1]) -len('short-notation:')
            else:
                tmp = len(item[1])
        else:
            tmp = len(item[1])
        if tmp > max_len_digilen:
            max_len_digilen = tmp
    max_len_digilen +=1
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

    calculate = []  # has to be calculated in this run
    prev_calculated = []  # was already calculated in a previous run
    try:
        store_confs
    except NameError:
        store_confs = []  # stores all confs which are sorted out!

    if config.solvent == "gas":
        print("Optimizing geometries at DFT level!")
    else:
        print("Optimizing geometries at DFT level with implicit solvation!")

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
            print(f"CONF{conf.id} is removed as requested by the user!")
            continue
        if conf.id > config.nconf:
            store_confs.append(conformers.pop(conformers.index(conf)))
            continue
        if conf.optimization_info["info"] == "not_calculated":
            conf = conformers.pop(conformers.index(conf))
            calculate.append(conf)
        elif conf.optimization_info["info"] == "failed":
            conf = conformers.pop(conformers.index(conf))
            store_confs.append(conf)
            print(f"Optimization of CONF{conf.id} failed in the previous run!")
        elif conf.optimization_info["info"] == "prep-failed":
            print(
                f"Preparation for the optimization of CONF{conf.id} failed in the "
                "previous run and is tried again!"
            )
            conf = conformers.pop(conformers.index(conf))
        elif conf.optimization_info["info"] == "calculated":
            conf = conformers.pop(conformers.index(conf))
            if conf.optimization_info["convergence"] == "converged":
                conf.job["success"] = True
                conf.job["ecyc"] = conf.optimization_info["ecyc"]
                if conf.optimization_info.get("cregen_sort", "pass") == "pass":
                    prev_calculated.append(conf)
                elif conf.optimization_info.get("cregen_sort", "pass") == "removed":
                    print(
                        f"CONF{conf.id} has been sorted out by CREGEN in a previous run."
                    )
                    store_confs.append(conf)
            else:
                # "not_converged" or "stopped_before_converged"
                store_confs.append(conf)
    if not calculate and not prev_calculated:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)
    if prev_calculated:
        check_for_folder(config.cwd, [i.id for i in prev_calculated], config.func)
        print("The optimization was performed before for:")
        print_block(["CONF" + str(i.id) for i in prev_calculated])
    pl = config.lenconfx + 4 + len(str("/" + config.func))

    instruction_prep = {
        "jobtype": "prep",
        "func": config.func,
        "basis": getattr(config, "basis", config.func_basis_default[config.func]),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": config.sm2,
        "optlevel": config.optlevel2,
        "copymos": "",
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
    }

    # INSTRUCTION OPT !!!!
    instruction_opt = {
        "func": config.func,
        "basis": getattr(config, "basis", config.func_basis_default[config.func]),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "fullopt": True,  # output to opt-part2.out
        "converged": False,
        "hlow": config.hlow,
        "sm": config.sm2,
        "optcycles": config.optcycles,
        "optlevel": config.optlevel2,
        "multiTemp": False,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
    }

    instruction_rrho_crude = {
        "jobtype": "rrhoxtb",
        "func": getattr(config, "part2_gfnv"),
        "gfn_version": getattr(config, "part2_gfnv"),
        "temperature": config.temperature,
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "progpath": config.external_paths["xtbpath"],
        "bhess": config.bhess,
        "sm_rrho": config.sm_rrho,
        "rmsdbias": config.rmsdbias,
        "cwd": config.cwd,
        "consider_sym": config.consider_sym,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
        "imagthr": config.imagthr,
        "sthr": config.sthr,
        "scale":config.scale,
    }
    instruction_rrho_crude["method"], _ = config.get_method_name(
        "rrhoxtb",
        bhess=config.bhess,
        gfn_version=instruction_rrho_crude["gfn_version"],
        sm=instruction_rrho_crude["sm_rrho"],
        solvent=instruction_rrho_crude["solvent"],
    )

    # Set optlevel and scfconv stuff ---------------------------------------
    # r2scan-3c has additional settings in tm_job._prep_cefine!
    if config.optlevel2 in ("crude", "sloppy", "loose"):
        instruction_prep["prepinfo"] = ["low"]
        instruction_opt["prepinfo"] = ["low"]
    elif config.optlevel2 == "lax":
        instruction_prep["prepinfo"] = ["low"]
        instruction_opt["prepinfo"] = ["low"]
    elif config.optlevel2 == "normal":
        instruction_prep["prepinfo"] = ["low+"]
        instruction_opt["prepinfo"] = ["low+"]
    elif config.optlevel2 in ("tight", "vtight", "extreme"):
        instruction_prep["prepinfo"] = ["high"]
        instruction_opt["prepinfo"] = ["high"]
    else:
        instruction_prep["prepinfo"] = ["low+"]
        instruction_opt["prepinfo"] = ["low+"]
    # -----------------------------------------------------------------------
    if config.ancopt:
        instruction_opt["jobtype"] = "xtbopt"
        instruction_opt["xtb_driver_path"] = config.external_paths["xtbpath"]
    else:
        instruction_opt["jobtype"] = "opt"

    if config.func == "r2scan-3c":
        instruction_prep["prepinfo"].extend(["-radsize", str(config.radsize)])

    instruction_opt["method"], _ = config.get_method_name(
        instruction_opt["jobtype"],
        func=instruction_opt["func"],
        basis=instruction_opt["basis"],
        solvent=instruction_opt["solvent"],
        sm=instruction_opt["sm"],
    )

    check = {True: "was successful", False: "FAILED"}
    if calculate:
        print("The optimization is calculated for:")
        print_block(["CONF" + str(i.id) for i in calculate])
        # create folders:
        save_errors, store_confs, calculate = new_folders(
            config.cwd, calculate, config.func, save_errors, store_confs
        )
        # write coord to folder
        calculate, store_confs, save_errors = ensemble2coord(
            config, config.func, calculate, store_confs, save_errors
        )

        # parallel prep execution
        calculate = run_in_parallel(
            config,
            q,
            resultq,
            job,
            config.maxthreads,
            config.omp,
            calculate,
            instruction_prep,
            config.balance,
            config.func,
        )
        # check if too many calculations failed

        for conf in list(calculate):
            if instruction_prep["jobtype"] == "prep":
                line = (
                    f"Preparation in {last_folders(conf.job['workdir'], 2):>{pl}} "
                    f"{check[conf.job['success']]}."
                )
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.optimization_info["info"] = "prep-failed"
                    store_confs.append(calculate.pop(calculate.index(conf)))

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
    # reset
    for conf in calculate:
        conf.reset_job_info()
    # ***************************************************************************
    # NEW ENSEMBLE OPTIMIZER:
    if calculate:
        print("Starting optimizations".center(70, "*"))
        ### settings in instruction_opt are overwriting conf.job everytime,(while loop)
        ### therefore dont write information which has to be reaccessed to it!

        run = 1
        timings = []  # time per cycle
        cycle_spearman = []  # spearmanthr used in evaluation per cycle
        nconf_cycle = []  # number of conformers at end of each cycle

        do_increase = 0.6
        if config.opt_limit * do_increase >= 1.5:
            ewin_increase = config.opt_limit * do_increase
            print(
                f"\nStarting threshold is set to {config.opt_limit} + "
                f"{do_increase*100} % = {config.opt_limit + ewin_increase} kcal/mol\n"
            )
        else:
            ewin_increase = 1.5
            print(
                f"\nStarting threshold is set to {config.opt_limit} + "
                f"{ewin_increase} kcal/mol = {config.opt_limit + ewin_increase} kcal/mol\n"
            )
        ewin_initial = config.opt_limit + ewin_increase
        ewin = config.opt_limit + ewin_increase

        print(f"Lower limit is set to G_thr(opt,2) = {config.opt_limit} kcal/mol\n")
        lower_limit = config.opt_limit
        maxecyc_prev = 1
        maxecyc = 0
        converged_run1 = []
        if config.nat > 200:
            # stopcycle = don't optimize more than stopcycle cycles
            stopcycle = config.nat * 2
        else:
            stopcycle = 200
        if config.opt_spearman:
            while calculate:
                tic = time.perf_counter()
                print(f"CYCLE {str(run)}".center(70, "*"))
                # if len(calculate) == 1:

                # if stopcycle - maxecyc <= 0:
                #     limit = config.optcycles
                # else:
                #     limit = stopcycle - maxecyc
                # instruction_opt["optcycles"] = limit

                # calculate batch of optimizations
                calculate = run_in_parallel(
                    config,
                    q,
                    resultq,
                    job,
                    config.maxthreads,
                    config.omp,
                    calculate,
                    instruction_opt,
                    config.balance,
                    config.func,
                )

                # check if optimization crashed
                for conf in list(calculate):
                    if not conf.job["success"]:
                        print(f"removing CONF{conf.id} because optimization crashed.")
                        conf.optimization_info["info"] = "failed"
                        conf.optimization_info["convergence"] = "not_converged"
                        conf.optimization_info["method"] = instruction_opt["method"]
                        store_confs.append(calculate.pop(calculate.index(conf)))
                    elif conf.job["success"]:
                        if conf.job["converged"]:
                            # don't optimize further:
                            print(
                                f"Geometry optimization converged for: "
                                f"CONF{conf.id} within {conf.job['cycles']:>3} cycles"
                            )
                            conf.optimization_info["info"] = "calculated"
                            conf.optimization_info["energy"] = conf.job["energy"]
                            conf.optimization_info["cycles"] = conf.job["cycles"]
                            conf.optimization_info["ecyc"] = conf.job["ecyc"]
                            conf.optimization_info["decyc"] = conf.job["decyc"]
                            conf.optimization_info["convergence"] = "converged"
                            conf.optimization_info["method"] = instruction_opt["method"]
                            if run == 1:
                                converged_run1.append(conf.id)
                            else:
                                prev_calculated.append(
                                    calculate.pop(calculate.index(conf))
                                )
                        else:
                            conf.optimization_info["energy"] = conf.job["energy"]
                            # optimization cycles didn't result in convergence
                            # optimize further
                if not calculate:
                    toc = time.perf_counter()
                    timings.append(toc - tic)
                    cycle_spearman.append("")
                    nconf_cycle.append(len(calculate) + len(prev_calculated))
                    break
                if run == 1 and calculate and config.evaluate_rrho:
                    # run GmRRHO on crudely optimized geometry
                    folder_rrho_crude = os.path.join(config.func, "rrho_crude")
                    # create folders:
                    save_errors, store_confs, calculate = new_folders(
                        config.cwd,
                        calculate,
                        folder_rrho_crude,
                        save_errors,
                        store_confs,
                    )
                    # copy optimized geoms to folder
                    for conf in list(calculate):
                        try:
                            tmp_from = os.path.join(
                                config.cwd, "CONF" + str(conf.id), config.func
                            )
                            tmp_to = os.path.join(
                                config.cwd, "CONF" + str(conf.id), folder_rrho_crude
                            )
                            shutil.copy(
                                os.path.join(tmp_from, "coord"),
                                os.path.join(tmp_to, "coord"),
                            )
                        except shutil.SameFileError:
                            pass
                        except FileNotFoundError:
                            if not os.path.isfile(os.path.join(tmp_from, "coord")):
                                print_errors(
                                    f"{'ERROR:':{WARNLEN}}while copying the coord file from {tmp_from}! "
                                    "The corresponding file does not exist.", save_errors
                                )
                            elif not os.path.isdir(tmp_to):
                                print_errors(
                                    f"{'ERROR:':{WARNLEN}}Could not create folder {tmp_to}!", save_errors
                                )
                            print_errors(f"{'ERROR:':{WARNLEN}}Removing conformer {conf.name}!", save_errors)
                            conf.lowlevel_grrho_info["info"] = "prep-failed"
                            store_confs.append(calculate.pop(calculate.index(conf)))
                    # parallel execution:
                    calculate = run_in_parallel(
                        config,
                        q,
                        resultq,
                        job,
                        config.maxthreads,
                        config.omp,
                        calculate,
                        instruction_rrho_crude,
                        config.balance,
                        folder_rrho_crude,
                    )
                    check = {True: "was successful", False: "FAILED"}
                    # check if too many calculations failed
                    ###
                    for conf in list(calculate):
                        print(
                            f"The G_mRRHO calculation on crudely optimized DFT "
                            f"geometry @ {conf.job['symmetry']} "
                            f"{check[conf.job['success']]} for "
                            f"{last_folders(conf.job['workdir'], 3):>{pl}}: "
                            f"{conf.job['energy']:>.8f}"
                        )
                        if not conf.job["success"]:
                            store_confs.append(calculate.pop(calculate.index(conf)))
                        else:
                            conf.optimization_info["energy_rrho"] = conf.job["energy"]
                            conf.optimization_info[
                                "method_rrho"
                            ] = instruction_rrho_crude["method"]
                            conf.optimization_info["info_rrho"] = "calculated"

                    for conf in list(calculate):
                        if conf.id in converged_run1:
                            prev_calculated.append(calculate.pop(calculate.index(conf)))
                    if not calculate:
                        toc = time.perf_counter()
                        timings.append(toc - tic)
                        cycle_spearman.append("")
                        nconf_cycle.append(len(calculate) + len(prev_calculated))
                        break
                if run >= 2 and config.crestcheck:
                    # do sorting with cregen!
                    calculate, prev_calculated, store_confs = crest_routine(
                        config,
                        calculate,
                        config.func,
                        store_confs,
                        prev_calculated=prev_calculated,
                    )
                    if not calculate:
                        toc = time.perf_counter()
                        timings.append(toc - tic)
                        cycle_spearman.append("")
                        nconf_cycle.append(len(calculate) + len(prev_calculated))
                        break
                maxecyc = max([len(conf.job["ecyc"]) for conf in calculate])
                print(f"Max number of performed iterations: {maxecyc}")
                if len(calculate + prev_calculated) == 1:
                    # can't do spearman with only one conf
                    run_spearman = False
                elif len(calculate) > 1 and run > 1:
                    run_spearman = True
                else:
                    run_spearman = True
                gesc = True  # gESC with already good sorting
                if run == 1 and gesc:
                    # only evaluate spearman starting from second cycle
                    print("Spearman rank evaluation is performed in the next cycle.")
                    cycle_spearman.append("")
                    run_spearman = False
                # if run == 1 and not gesc:
                #     # only evaluate spearman starting from second cycle
                #     print("Spearman rank evaluation is performed in the next cycle.")
                #     cycle_spearman.append("")
                #     run_spearman = False
                #     run += 1
                #     toc = time.perf_counter()
                #     timings.append(toc - tic)
                #     nconf_cycle.append(len(calculate) + len(prev_calculated))
                #     print(f"CYCLE {run} performed in {toc -tic:0.4f} seconds")
                #     continue

                # lists of equal lenght:
                for conf in sorted(
                    calculate + prev_calculated, key=lambda x: int(x.id)
                ):
                    if len(conf.job["ecyc"]) < maxecyc:
                        for _ in range(maxecyc - len(conf.job["ecyc"])):
                            conf.job["ecyc"].append(conf.job["ecyc"][-1])

                # calculate min of each cycle:
                minecyc = []
                if config.evaluate_rrho:
                    rrho_energy = "energy_rrho"
                else:
                    rrho_energy = "axqzv"  # to get 0.0 contribution
                for i in range(maxecyc):
                    try:
                        minecyc.append(
                            min(
                                [
                                    conf.job["ecyc"][i]
                                    + getattr(conf, "optimization_info").get(
                                        rrho_energy, 0.0
                                    )
                                    for conf in calculate + prev_calculated
                                    if conf.job["ecyc"][i] is not None
                                ]
                            )
                        )
                    except (ValueError) as e:
                        minecyc.append(0.0)
                        print(e)
                # evalulate ΔE
                for conf in sorted(
                    calculate + prev_calculated, key=lambda x: int(x.id)
                ):
                    conf.job["decyc"] = []
                    for i in range(maxecyc):
                        conf.job["decyc"].append(
                            (
                                conf.job["ecyc"][i]
                                + getattr(conf, "optimization_info").get(
                                    rrho_energy, 0.0
                                )
                                - minecyc[i]
                            )
                            * AU2KCAL
                        )
                if run == 1:
                    print("")
                    for conf in sorted(
                        calculate + prev_calculated, key=lambda x: int(x.id)
                    ):
                        print(
                            f"CONF{conf.id :<{config.lenconfx}} initial ΔG =  "
                            f"{conf.job['decyc'][0]:^5.2f} kcal/mol and "
                            f"current ΔG =  {conf.job['decyc'][-1]:^5.2f} kcal/mol."
                            f" ({conf.optimization_info['convergence']})"
                        )
                    previouscycle = maxecyc
                    print("")
                else:
                    print("")
                    for conf in sorted(
                        calculate + prev_calculated, key=lambda x: int(x.id)
                    ):
                        print(
                            f"CONF{conf.id :<{config.lenconfx}} previous ΔG =  "
                            f"{conf.job['decyc'][previouscycle-1]:^5.2f} kcal/mol and "
                            f"current ΔG =  {conf.job['decyc'][-1]:^5.2f} kcal/mol."
                            f" ({conf.optimization_info['convergence']})"
                        )
                    previouscycle = maxecyc
                    print("")
                if run_spearman:
                    num_eval = 3
                    try:
                        toevaluate = []
                        for i in range(maxecyc_prev, maxecyc):
                            if i + num_eval <= maxecyc:
                                toevaluate.append(i)
                        _ = max(toevaluate)
                        digits1 = 4
                    except ValueError:
                        # need to do another optimization cycle
                        run += 1
                        toc = time.perf_counter()
                        timings.append(toc - tic)
                        cycle_spearman.append("")
                        nconf_cycle.append(len(calculate) + len(prev_calculated))
                        print(f"CYCLE {run} performed in {toc -tic:0.4f} seconds")
                        continue

                    evalspearman = []
                    for i in toevaluate:
                        deprevious = [
                            conf.job["decyc"][i - 1]
                            for conf in sorted(
                                calculate + prev_calculated, key=lambda x: int(x.id)
                            )
                        ]
                        decurrent = [
                            conf.job["decyc"][i - 1 + num_eval]
                            for conf in sorted(
                                calculate + prev_calculated, key=lambda x: int(x.id)
                            )
                        ]

                        spearman_v = spearman(deprevious, decurrent)
                        if i in toevaluate[-2:]:
                            print(
                                f"Evaluating Spearman coeff. from {i:>{digits1}} --> "
                                f"{i+num_eval:>{digits1}}"
                                f" =  {spearman_v:>.4f}"
                            )
                            evalspearman.append(spearman_v)
                        else:
                            print(
                                f"{'':>10} Spearman coeff. from {i:>{digits1}} --> "
                                f"{i+num_eval:>{digits1}}"
                                f" =  {spearman_v:>.4f}"
                            )
                    print(
                        f"Final averaged Spearman correlation coefficient: "
                        f"{(sum(evalspearman)/2):>.4f}"
                    )

                    if (
                        len(evalspearman) >= 2
                        and sum(evalspearman) / 2 >= config.spearmanthr
                    ):
                        print("\nPES is assumed to be parallel")
                        # adjust threshold Ewin:
                        if ewin > lower_limit:
                            if (ewin - (ewin_increase / 3)) < lower_limit:
                                ewin = lower_limit
                            else:
                                ewin += -(ewin_increase / 3)
                            print(
                                f"Updated optimization threshold to: {ewin:.2f} kcal/mol"
                            )
                        else:
                            print(
                                f"Current optimization threshold: {ewin:.2f} kcal/mol"
                            )
                        cycle_spearman.append(f"{sum(evalspearman)/2:.3f}")

                for conf in list(calculate):
                    if conf.job["decyc"][-1] > ewin and conf.job["grad_norm"] < 0.01:
                        print(
                            f"CONF{conf.id} is above {ewin} kcal/mol and gradient "
                            f"norm ({conf.job['grad_norm']}) is below {0.01}."
                        )
                        if conf.job["decyc"][-1] < ewin_initial:
                            print(
                                f"CONF{conf.id} is removed because of the "
                                "lowered threshold!"
                            )
                        # transfer energies and remove conf
                        conf.optimization_info["energy"] = conf.job["energy"]
                        conf.optimization_info["info"] = "calculated"
                        conf.optimization_info["cycles"] = conf.job["cycles"]
                        conf.optimization_info["ecyc"] = conf.job["ecyc"]
                        conf.optimization_info["decyc"] = conf.job["decyc"]
                        conf.optimization_info["method"] = instruction_opt["method"]
                        conf.optimization_info[
                            "convergence"
                        ] = "stopped_before_converged"
                        print(
                            f"CONF{conf.id} is above threshold, dont optimize "
                            f"further and remove conformer."
                        )
                        store_confs.append(calculate.pop(calculate.index(conf)))
                    elif conf.job["decyc"][-1] > ewin and conf.job["grad_norm"] > 0.01:
                        print(
                            f"CONF{conf.id} is above {ewin} kcal/mol but "
                            f"gradient norm ({conf.job['grad_norm']}) is "
                            f"above {0.01} --> not sorted out!"
                        )
                toc = time.perf_counter()
                timings.append(toc - tic)
                nconf_cycle.append(len(calculate) + len(prev_calculated))
                print(f"\nCYCLE {run} performed in { toc - tic:0.4f} seconds")
                #
                if maxecyc >= stopcycle:
                    print("")
                    for conf in list(calculate):
                        # don't optimize further:
                        print(
                            f"!!! Geometry optimization STOPPED because of "
                            f"optcycle limit of {stopcycle} cycles reached for: "
                            f"CONF{conf.id} within {conf.job['cycles']:>3} cycles"
                        )
                        conf.optimization_info["info"] = "calculated"
                        conf.optimization_info["energy"] = conf.job["energy"]
                        conf.optimization_info["cycles"] = conf.job["cycles"]
                        conf.optimization_info["ecyc"] = conf.job["ecyc"]
                        conf.optimization_info["decyc"] = conf.job["decyc"]
                        conf.optimization_info[
                            "convergence"
                        ] = "converged"  #### THIS IS NOT CORRECT /or can lead to errors!
                        conf.optimization_info["method"] = instruction_opt["method"]
                        prev_calculated.append(calculate.pop(calculate.index(conf)))
                #
                maxecyc_prev = maxecyc
                run += 1
            # END while loop
        else:
            # use standard optimization!
            # update instruct_opt
            tic = time.perf_counter()
            del instruction_opt["optcycles"]
            # calculate first round
            calculate = run_in_parallel(
                config,
                q,
                resultq,
                job,
                config.maxthreads,
                config.omp,
                calculate,
                instruction_opt,
                config.balance,
                config.func,
            )
            # check if optimization crashed
            for conf in list(calculate):
                if not conf.job["success"]:
                    print(f"removing CONF{conf.id} because optimization crashed.")
                    conf.optimization_info["info"] = "failed"
                    conf.optimization_info["convergence"] = "not_converged"
                    conf.optimization_info["method"] = instruction_opt["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                elif conf.job["success"]:
                    if conf.job["converged"]:
                        # don't optimize further:
                        conf.optimization_info["info"] = "calculated"
                        conf.optimization_info["energy"] = conf.job["energy"]
                        conf.optimization_info["cycles"] = conf.job["cycles"]
                        conf.optimization_info["ecyc"] = conf.job["ecyc"]
                        conf.optimization_info["decyc"] = conf.job["decyc"]
                        conf.optimization_info["convergence"] = "converged"
                        conf.optimization_info["method"] = instruction_opt["method"]
                        # prev_calculated to keep it consistent with new ensemble optimizer
                        prev_calculated.append(calculate.pop(calculate.index(conf)))
                    else:
                        print(f"{'ERROR:':{WARNLEN}}CONF{conf.id} fell through sorting")
            toc = time.perf_counter()
            timings.append(toc - tic)
            # ********************end standard optimization *********************
        print("Finished optimizations!".center(70, "*"))
        if config.opt_spearman and (
            len(timings) == len(cycle_spearman) == len(nconf_cycle)
        ):
            try:
                tl = max([len(f"{i: .2f}") for i in timings])
                if tl > 7:
                    tmp1 = tl
                else:
                    tmp1 = 7
                print("Timings:")
                print(f"Cycle:  [s]  {'#nconfs':^{tmp1}}  Spearman coeff.")
                for i in range(len(timings)):
                    print(
                        f"{i+1:>4}  {timings[i]:> {tl}.2f}  {nconf_cycle[i]:^{tmp1}}  {cycle_spearman[i]}"
                    )
                print("sum:  {:> .2f}".format(sum(timings)))
            except Exception as e:
                print(e)
        else:
            print("Timings:")
            print("Cycle:  [s]")
            for i in timings:
                print("{:4}  {:>.2f}".format(timings.index(i) + 1, i))
            print("sum:  {:>.2f}".format(sum(timings)))
        print("\nCONVERGED optimizations for the following remaining conformers:")
        prev_calculated.sort(key=lambda x: int(x.id))
    # end if calculate--
    for conf in list(prev_calculated):
        print(
            f"Converged optimization for {'CONF' + str(conf.id):{config.lenconfx+4}} "
            f"after {conf.optimization_info['cycles'] :>3} cycles: "
            f"{conf.optimization_info['energy']:>.7f}"
        )
        calculate.append(prev_calculated.pop(prev_calculated.index(conf)))

    ensembledata.nconfs_per_part["part2_opt"] = len(calculate)
    # save current data to jsonfile
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )
    if config.crestcheck:
        calculate, prev_calculated, store_confs = crest_routine(
            config, calculate, config.func, store_confs
        )

    # reset
    for conf in calculate:
        conf.reset_job_info()
    if not calculate and not prev_calculated:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)
    # ******************************Optimization done****************************
    # Start Gsolv (COSMO-RS) calculation (or only gas phase single-point)
    instruction_gsolv = {
        "func": config.func,
        "prepinfo": ["low+"],  # TM m4 scfconv6
        "basis": getattr(config, "basis", config.func_basis_default[config.func]),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": config.smgsolv2,
        "temperature": config.temperature,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
        "gfn_version": config.part2_gfnv,
    }
    if config.multitemp:
        instruction_gsolv["trange"] = [
            i for i in frange(config.trange[0], config.trange[1], config.trange[2])
        ]
    else:
        instruction_gsolv["trange"] = []
    if config.solvent == "gas":
        print("\nCalculating single-point energies!")
        instruction_gsolv["jobtype"] = "sp"
        instruction_gsolv["method"], _ = config.get_method_name(
            instruction_gsolv["jobtype"],
            func=instruction_gsolv["func"],
            basis=instruction_gsolv["basis"],
        )
        folder = instruction_gsolv["func"]
        name = "lowlevel single-point"
    else:
        print(
            "\nCalculating single-point energies and solvation contribution (G_solv)!"
        )
        if config.smgsolv2 in config.smgsolv_2:
            # additive Gsolv
            # COSMO-RS
            if "cosmors" in config.smgsolv2 and config.smgsolv2 != "dcosmors":
                job = TmJob
                instruction_gsolv["prepinfo"] = ["low+"]
                exc_fine = {"cosmors": "normal", "cosmors-fine": "fine"}
                tmp = {
                    "jobtype": "cosmors",
                    "cosmorssetup": config.external_paths["cosmorssetup"],
                    "cosmorsparam": exc_fine.get(config.smgsolv2, "normal"),
                    "cosmothermversion": config.external_paths["cosmothermversion"],
                    "ctd-param": config.cosmorsparam,
                    "copymos": str(instruction_gsolv["func"]),
                }
                instruction_gsolv.update(tmp)
                instruction_gsolv["method"], instruction_gsolv[
                    "method2"
                ] = config.get_method_name(
                    "cosmors",
                    func=instruction_gsolv["func"],
                    basis=instruction_gsolv["basis"],
                    sm=instruction_gsolv["sm"],
                )
                folder = str(instruction_gsolv["func"]) + "/COSMO"
                name = "lowlevel COSMO-RS"
            # GBSA-Gsolv / ALPB-Gsolv
            elif config.smgsolv2 in ("gbsa_gsolv", "alpb_gsolv"):
                instruction_gsolv["jobtype"] = instruction_gsolv["sm"]
                if config.prog == "orca":
                    instruction_gsolv["progpath"] = config.external_paths["orcapath"]
                instruction_gsolv["xtb_driver_path"] = config.external_paths["xtbpath"]
                instruction_gsolv["method"], instruction_gsolv[
                    "method2"
                ] = config.get_method_name(
                    instruction_gsolv["jobtype"],
                    func=instruction_gsolv["func"],
                    basis=instruction_gsolv["basis"],
                    sm=instruction_gsolv["sm"],
                    gfn_version=instruction_gsolv["gfn_version"],
                )
                if (
                    conf.lowlevel_sp_info["info"] == "calculated"
                    and conf.lowlevel_sp_info["method"] == instruction_gsolv["method"]
                ):
                    # do not calculate gas phase sp again!
                    instruction_gsolv["energy"] = conf.lowlevel_sp_info["energy"]
                    instruction_gsolv["prepinfo"] = []
                # else:
                name = "lowlevel additive solvation"
                folder = str(instruction_gsolv["func"]) + "/Gsolv2"
            # SMD_Gsolv
            elif config.smgsolv2 == "smd_gsolv":
                job = OrcaJob
                instruction_gsolv["jobtype"] = "smd_gsolv"
                instruction_gsolv["progpath"] = config.external_paths["orcapath"]
                instruction_gsolv["method"], instruction_gsolv[
                    "method2"
                ] = config.get_method_name(
                    "smd_gsolv",
                    func=instruction_gsolv["func"],
                    basis=instruction_gsolv["basis"],
                    sm=instruction_gsolv["sm"],
                )
                name = "lowlevel SMD_Gsolv"
                folder = str(instruction_gsolv["func"]) + "/Gsolv2"
        else:
            # with implicit solvation
            instruction_gsolv["jobtype"] = "sp_implicit"
            if config.prog == "orca":
                instruction_gsolv["progpath"] = config.external_paths["orcapath"]
            instruction_gsolv["method"], instruction_gsolv[
                "method2"
            ] = config.get_method_name(
                "sp_implicit",
                func=instruction_gsolv["func"],
                basis=instruction_gsolv["basis"],
                sm=instruction_gsolv["sm"],
            )
            name = "lowlevel single-point"
            folder = instruction_gsolv["func"]

    for conf in list(calculate):
        if conf.removed:
            store_confs.append(calculate.pop(calculate.index(conf)))
            print(f"CONF{conf.id} is removed as requested by the user!")
            continue
        if conf.lowlevel_sp_info["info"] == "failed":
            conf = calculate.pop(calculate.index(conf))
            store_confs.append(conf)
            print(f"Calculation of CONF{conf.id} failed in the previous run!")
        elif conf.lowlevel_sp_info["info"] == "not_calculated":
            # has to be calculated now
            #  take opt sp as lowlevel sp !
            if conf.optimization_info["method"] == instruction_gsolv["method"]:
                if conf.optimization_info.get("method", "not_found") == "calculated":
                    conf = calculate.pop(calculate.index(conf))
                    conf.lowlevel_sp_info["info"] = "calculated"
                    conf.lowlevel_sp_info["method"] = instruction_gsolv["method"]
                    conf.lowlevel_sp_info["energy"] = conf.optimization_info["energy"]
                    conf.job["success"] = True
                    prev_calculated.append(conf)
                    continue
        elif conf.lowlevel_sp_info["info"] == "prep-failed":
            print(
                f"Preparation step for CONF{conf.id} failed in the previous "
                "run and is retried now!"
            )
            # is retried now!
        elif conf.lowlevel_sp_info["info"] == "calculated":
            conf = calculate.pop(calculate.index(conf))
            if config.solvent != "gas":
                # check if solvation calculation is calculated as well
                if conf.lowlevel_gsolv_info["info"] == "failed":
                    store_confs.append(conf)
                    print(
                        f"Calculation of the solvation contribution for CONF"
                        f"{conf.id} failed in the previous run!"
                    )
                elif conf.lowlevel_gsolv_info["info"] == "not_calculated":
                    calculate.append(conf)
                elif conf.lowlevel_gsolv_info["info"] == "calculated":
                    conf.job["success"] = True
                    prev_calculated.append(conf)
                else:
                    print("UNEXPECTED BEHAVIOUR")
            elif config.solvent == "gas":
                conf.job["success"] = True
                prev_calculated.append(conf)
        else:
            print("\nMISSING STUFF!\n")

    if not calculate and not prev_calculated:
        print(f"{'ERROR:':{WARNLEN}}No conformers left!")
        print("Going to exit!")
        sys.exit(1)
    if prev_calculated:
        check_for_folder(config.cwd, [i.id for i in prev_calculated], folder)
        if config.solvent == "gas":
            print("The low level single-point was calculated before for:")
        else:
            print("The low level gsolv calculation was calculated before for:")
        print_block(["CONF" + str(i.id) for i in prev_calculated])
    pl = config.lenconfx + 4 + len(str("/" + folder))

    check = {True: "was successful", False: "FAILED"}
    if calculate:
        if config.solvent == "gas":
            print("The low level single-point is now calculated for:")
        if config.solvent != "gas" and config.smgsolv2 in config.smgsolv_2:
            print("The low level gsolv calculation is now calculated for:")
            # need to create folders
            save_errors, store_confs, calculate = new_folders(
                config.cwd, calculate, folder, save_errors, store_confs
            )
            # need to copy optimized coord to COSMO/GSOLV2 folder
            for conf in calculate:
                tmp1 = os.path.join(
                    config.cwd,
                    "CONF" + str(conf.id),
                    instruction_gsolv["func"],
                    "coord",
                )
                tmp2 = os.path.join("CONF" + str(conf.id), folder, "coord")
                try:
                    shutil.copy(tmp1, tmp2)
                except FileNotFoundError:
                    print(f"{'ERROR:':{WARNLEN}}can't copy optimized geometry!")
        print_block(["CONF" + str(i.id) for i in calculate])
        # parallel execution:
        calculate = run_in_parallel(
            config,
            q,
            resultq,
            job,
            config.maxthreads,
            config.omp,
            calculate,
            instruction_gsolv,
            config.balance,
            folder,
        )

        for conf in list(calculate):
            if instruction_gsolv["jobtype"] in("sp", "sp_implicit"):
                line = (
                    f"{name} calculation {check[conf.job['success']]}"
                    f" for {last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f}"
                )
                print(line)
                if not conf.job["success"]:
                    save_errors.append(line)
                    conf.lowlevel_sp_info["info"] = "failed"
                    conf.lowlevel_sp_info["method"] = conf.job["method"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.lowlevel_sp_info["energy"] = conf.job["energy"]
                    conf.lowlevel_sp_info["info"] = "calculated"
                    conf.lowlevel_sp_info["method"] = conf.job["method"]
                    if instruction_gsolv["jobtype"] == "sp_implicit":
                        conf.lowlevel_gsolv_info["energy"] = 0.0
                        conf.lowlevel_gsolv_info["info"] = "calculated"
                        conf.lowlevel_gsolv_info["method"] = conf.job["method2"]
            elif instruction_gsolv["jobtype"] in (
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
                    conf.lowlevel_sp_info["info"] = "failed"
                    conf.lowlevel_sp_info["method"] = conf.job["method"]
                    conf.lowlevel_gsolv_info["info"] = "failed"
                    conf.lowlevel_gsolv_info["method"] = conf.job["method2"]
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.lowlevel_sp_info["energy"] = conf.job["energy"]
                    conf.lowlevel_sp_info["info"] = "calculated"
                    conf.lowlevel_sp_info["method"] = instruction_gsolv["method"]
                    conf.lowlevel_gsolv_info["energy"] = conf.job["energy2"]
                    conf.lowlevel_gsolv_info["gas-energy"] = conf.job["energy"]
                    conf.lowlevel_gsolv_info["info"] = "calculated"
                    conf.lowlevel_gsolv_info["method"] = instruction_gsolv["method2"]
                    conf.lowlevel_gsolv_info["range"] = conf.job["erange1"]
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
        for conf in list(prev_calculated):
            conf.job["workdir"] = os.path.normpath(
                os.path.join(config.cwd, "CONF" + str(conf.id), folder)
            )
            if instruction_gsolv["jobtype"] in ("sp", "sp_implicit"):
                print(
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.lowlevel_sp_info['energy']:>.7f}"
                )
            elif instruction_gsolv["jobtype"] in (
                "cosmors",
                "smd_gsolv",
                "gbsa_gsolv",
                "alpb_gsolv",
            ):
                print(
                    f"{name} calculation {check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 3):>{pl}}: "
                    f"{conf.lowlevel_gsolv_info['energy']:>.7f}"
                )
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
    # reset
    for conf in calculate:
        conf.reset_job_info()
    if not calculate:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)

    # ***************************************************************************
    # Starting grrho calculation on DFT geometry (bhess)
    if config.evaluate_rrho:
        if config.solvent == "gas":
            print("\nCalculating lowlevel G_mRRHO on DFT geometry!")
        else:
            print(
                "\nCalculating lowlevel G_mRRHO with implicit solvation "
                "on DFT geometry!"
            )
        for conf in list(calculate):
            if conf.lowlevel_grrho_info["info"] == "not_calculated":
                pass
            if conf.lowlevel_grrho_info["info"] == "prep-failed":
                # try again
                pass
            elif conf.lowlevel_grrho_info["info"] == "failed":
                conf = calculate.pop(calculate.index(conf))
                conf.__class__ = job
                store_confs.append(conf)
                print(f"Calculation of CONF{conf.id} failed in the previous run!")
            elif conf.lowlevel_grrho_info["info"] == "calculated":
                conf = calculate.pop(calculate.index(conf))
                conf.__class__ = job
                conf.job["success"] = True
                prev_calculated.append(conf)

        if not calculate and not prev_calculated:
            print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
            print("Going to exit!")
            sys.exit(1)
        folderrho = "rrho_part2"
        if prev_calculated:
            check_for_folder(config.cwd, [i.id for i in prev_calculated], folderrho)
            print("The G_mRRHO calculation was performed before for:")
            print_block(["CONF" + str(i.id) for i in prev_calculated])
        pl = config.lenconfx + 4 + len(str("/" + folderrho))
        instruction_rrho = {
            "jobtype": "rrhoxtb",
            "func": getattr(config, "part2_gfnv"),
            "gfn_version": getattr(config, "part2_gfnv"),
            "temperature": config.temperature,
            "charge": config.charge,
            "unpaired": config.unpaired,
            "solvent": config.solvent,
            "progpath": config.external_paths["xtbpath"],
            "bhess": config.bhess,
            "sm_rrho": config.sm_rrho,
            "rmsdbias": config.rmsdbias,
            "cwd": config.cwd,
            "consider_sym": config.consider_sym,
            "energy": 0.0,
            "energy2": 0.0,
            "success": False,
            "imagthr": config.imagthr,
            "sthr": config.sthr,
            "scale":config.scale,
        }
        instruction_rrho["method"], _ = config.get_method_name(
            "rrhoxtb",
            bhess=config.bhess,
            gfn_version=instruction_rrho["gfn_version"],
            sm=instruction_rrho["sm_rrho"],
            solvent=instruction_rrho["solvent"],
        )
        if config.multitemp:
            instruction_rrho["trange"] = [
                i for i in frange(config.trange[0], config.trange[1], config.trange[2])
            ]
        else:
            instruction_rrho["trange"] = []

        if calculate:
            print("The lowlevel G_mRRHO calculation is now performed for:")
            print_block(["CONF" + str(i.id) for i in calculate])
            # create folders:
            save_errors, store_confs, calculate = new_folders(
                config.cwd, calculate, folderrho, save_errors, store_confs
            )
            # copy optimized geoms to folder
            for conf in list(calculate):
                try:
                    tmp_from = os.path.join(
                        config.cwd, "CONF" + str(conf.id), config.func
                    )
                    tmp_to = os.path.join(config.cwd, "CONF" + str(conf.id), folderrho)
                    shutil.copy(
                        os.path.join(tmp_from, "coord"), os.path.join(tmp_to, "coord")
                    )
                except shutil.SameFileError:
                    pass
                except FileNotFoundError:
                    if not os.path.isfile(os.path.join(tmp_from, "coord")):
                        print_errors(
                            f"{'ERROR:':{WARNLEN}}while copying the coord file from {tmp_from}! "
                            "The corresponding file does not exist.",
                            save_errors
                        )
                    elif not os.path.isdir(tmp_to):
                        print_errors(f"{'ERROR:':{WARNLEN}}Could not create folder {tmp_to}!", save_errors)
                    print_errors(f"{'ERROR:':{WARNLEN}}Removing conformer {conf.name}!", save_errors)
                    conf.lowlevel_grrho_info["info"] = "prep-failed"
                    store_confs.append(calculate.pop(calculate.index(conf)))
            # parallel execution:
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
                folderrho,
            )
            check = {True: "was successful", False: "FAILED"}
            # check if too many calculations failed

            ###
            for conf in list(calculate):
                print(
                    f"The lowlevel G_mRRHO calculation @ {conf.job['symmetry']} "
                    f"{check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.job['energy']:>.8f}"
                )
                if not conf.job["success"]:
                    conf.lowlevel_grrho_info["info"] = "failed"
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.sym = conf.job["symmetry"]
                    conf.lowlevel_grrho_info["rmsd"] = conf.job["rmsd"]
                    conf.lowlevel_grrho_info["energy"] = conf.job["energy"]
                    conf.lowlevel_grrho_info["info"] = "calculated"
                    conf.lowlevel_grrho_info["method"] = instruction_rrho["method"]
                    conf.lowlevel_grrho_info["range"] = conf.job["erange1"]
                    conf.lowlevel_hrrho_info["range"] = conf.job["erange2"]
                    conf.lowlevel_hrrho_info["info"] = "calculated"
                    conf.lowlevel_hrrho_info["method"] = instruction_rrho["method"]
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
                    f"The lowlevel G_mRRHO calculation @ {conf.sym} "
                    f"{check[conf.job['success']]} for "
                    f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                    f"{conf.lowlevel_grrho_info['energy']:>.8f}"
                )
                calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
        if not calculate:
            print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
            print("Going to exit!")
            sys.exit(1)

    # printout for part2 -------------------------------------------------------
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("* Gibbs free energies of part2 *".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "xtb_energy"),
        lambda conf: getattr(conf, "rel_xtb_energy"),
        lambda conf: getattr(conf, "lowlevel_sp_info")["energy"],
        lambda conf: getattr(conf, "lowlevel_gsolv_info")["energy"],
        lambda conf: getattr(conf, "lowlevel_grrho_info")["energy"],
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
        "",  # Gsolv
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
        (5, 2),
        (5, 2),
    ]
    if config.solvent == "gas":
        # Energy
        columndescription[3] = instruction_gsolv["method"]
    elif config.solvent != "gas":
        # Energy
        columndescription[3] = instruction_gsolv["method"]
        # Gsolv
        columndescription[4] = instruction_gsolv["method2"]
    if config.evaluate_rrho:
        # Grrho
        columndescription[5] = instruction_rrho["method"]
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho:
            # ignore rrho in printout
            columncall.pop(5)
            columnheader.pop(5)
            columndescription.pop(5)
            columnformat.pop(5)
        if config.solvent == "gas":
            # ignore Gsolv
            columncall.pop(4)
            columnheader.pop(4)
            columndescription.pop(4)
            columnformat.pop(4)

    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "lowlevel_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "lowlevel_gsolv_info"
        e = "lowlevel_sp_info"
        conf.calc_free_energy(e=e, solv=solv, rrho=rrho)

    calculate = calc_boltzmannweights(calculate, "free_energy", config.temperature)
    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
    except ValueError:
        raise ValueError
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
        if conf.free_energy == minfree:
            lowestconf = conf.id
            ensembledata.bestconf["part2"] = conf.id

    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part2.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
        columndescription2=columndescription2,
    )
    # printout for part2 -------------------------------------------------------

    # ***************************************************************************
    # SD on solvation:
    # DCOSMO-RS_GSOLV
    instruction_gsolv_compare = {
        "func": config.func,
        "basis": getattr(config, "basis", config.func_basis_default[config.func]),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": "alpb_gsolv",
        "temperature": config.temperature,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
        "gfn_version": config.part2_gfnv,
    }
    instruction_gsolv_compare["trange"] = []
    instruction_gsolv_compare["prepinfo"] = []
    instruction_gsolv_compare["xtb_driver_path"] = config.external_paths["xtbpath"]
    instruction_gsolv_compare["jobtype"] = instruction_gsolv_compare["sm"]
    _, instruction_gsolv_compare["method"] = config.get_method_name(
        instruction_gsolv_compare["jobtype"],
        func=instruction_gsolv_compare["func"],
        basis=instruction_gsolv_compare["basis"],
        sm=instruction_gsolv_compare["sm"],
        gfn_version=instruction_gsolv_compare["gfn_version"],
    )
    folder_compare = "alpb_gsolv"
    name = "alpb_gsolv".upper()
    pl = config.lenconfx + 4 + len(str("/" + folder_compare))
    if (
        config.solvent != "gas"
        and config.sm2 == "dcosmors"
        and config.smgsolv2 in ("cosmors", "cosmors-fine")
    ):
        dorun = True
        while dorun:
            print(
                "\nCalculating ALPB_Gsolv values for evaluation of the std. dev. of Gsolv."
            )
            for conf in calculate:
                if conf.id == ensembledata.bestconf["part2"]:
                    gsolv_min = conf.lowlevel_gsolv_info["energy"]
                    gsolv_min_id = conf.id
                    try:
                        dcosmors_gsolv_min = (
                            conf.optimization_info["energy"]
                            - conf.lowlevel_gsolv_info["gas-energy"]
                        )
                    except (TypeError, KeyError):
                        print_errors(
                            f"{'ERROR:':{WARNLEN}}Can't calculate DCOSMO-RS_gsolv. Skipping SD of Gsolv!", save_errors
                        )
                        dorun = False
                        break

            for conf in list(calculate):
                if conf.lowlevel_gsolv_compare_info["info"] == "not_calculated":
                    pass
                if conf.lowlevel_gsolv_compare_info["info"] == "prep-failed":
                    # try again
                    pass
                elif conf.lowlevel_gsolv_compare_info["info"] == "failed":
                    # dont remove conformer this is only to calculate SD
                    conf = calculate.pop(calculate.index(conf))
                    conf.job["success"] = True
                    conf.lowlevel_gsolv_compare_info["energy"] = 0.0
                    prev_calculated.append(conf)
                    print(f"Calculation of CONF{conf.id} failed in the previous run!")
                elif conf.lowlevel_gsolv_compare_info["info"] == "calculated":
                    conf = calculate.pop(calculate.index(conf))
                    conf.job["success"] = True
                    prev_calculated.append(conf)
            # need to create folders
            save_errors, store_confs, calculate = new_folders(
                config.cwd, calculate, folder_compare, save_errors, store_confs
            )
            # need to copy optimized coord to COSMO/GSOLV2 folder
            for conf in calculate:
                tmp1 = os.path.join(
                    config.cwd,
                    "CONF" + str(conf.id),
                    instruction_gsolv_compare["func"],
                    "coord",
                )
                tmp2 = os.path.join("CONF" + str(conf.id), folder_compare, "coord")
                try:
                    shutil.copy(tmp1, tmp2)
                except FileNotFoundError:
                    print_errors(f"{'ERROR:':{WARNLEN}}can't copy optimized geometry!", save_errors)

            if calculate:
                calculate = run_in_parallel(
                    config,
                    q,
                    resultq,
                    job,
                    config.maxthreads,
                    config.omp,
                    calculate,
                    instruction_gsolv_compare,
                    config.balance,
                    folder_compare,
                )
                for conf in calculate:
                    line = (
                        f"{name} calculation {check[conf.job['success']]} for "
                        f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                        f"{conf.job['energy2']:>.8f}"
                    )
                    print(line)
                    if not conf.job["success"]:
                        save_errors.append(line)
                        conf.lowlevel_gsolv_compare_info["info"] = "failed"
                        conf.lowlevel_gsolv_compare_info["method"] = conf.job["method"]
                        # store_confs.append(calculate.pop(calculate.index(conf)))
                        print("ERROR")
                    else:
                        conf.lowlevel_gsolv_compare_info["energy"] = conf.job["energy2"]
                        conf.lowlevel_gsolv_compare_info["info"] = "calculated"
                        conf.lowlevel_gsolv_compare_info["method"] = instruction_gsolv[
                            "method"
                        ]
            if prev_calculated:
                for conf in list(prev_calculated):
                    conf.job["workdir"] = os.path.normpath(
                        os.path.join(config.cwd, "CONF" + str(conf.id), folder_compare)
                    )
                    line = (
                        f"{name} calculation {check[conf.job['success']]} for "
                        f"{last_folders(conf.job['workdir'], 2):>{pl}}: "
                        f"{conf.lowlevel_gsolv_compare_info['energy']:>.8f}"
                    )
                    print(line)
                    calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
            for conf in calculate:
                if conf.id == gsolv_min_id:
                    alpb_gsolv_min = conf.lowlevel_gsolv_compare_info["energy"]

            print("\nSD of solvation models (all units in kcal/mol):")
            print(
                f"CONFX ΔG(COSMO-RS) ΔG(DCOSMO-RS_gsolv) ΔG(ALPB_gsolv)  "
                f"SD(COSMO-RS 40%, DCOSMO-RS_gsolv 40%, ALPB_gsolv 20%)"
            )
            print("".ljust(PLENGTH, "-"))
            pl = max([len(str(conf.id)) for conf in calculate])
            for conf in calculate:
                dgsolv = -gsolv_min + conf.lowlevel_gsolv_info["energy"]
                dgdcosmors = -dcosmors_gsolv_min + (
                    conf.optimization_info["energy"]
                    - conf.lowlevel_gsolv_info["gas-energy"]
                )
                dgalpb = -alpb_gsolv_min + conf.lowlevel_gsolv_compare_info["energy"]
                conf.lowlevel_gsolv_compare_info["std_dev"] = calc_weighted_std_dev(
                    [dgsolv, dgdcosmors, dgalpb], weights=[0.4, 0.4, 0.2]
                )
                print(
                    f"CONF{conf.id:<{pl}} {dgsolv*AU2KCAL:^ 12.2f} "
                    f"{dgdcosmors*AU2KCAL:^ 19.2f} {dgalpb*AU2KCAL:^ 14.2f}"
                    f" {conf.lowlevel_gsolv_compare_info['std_dev']*AU2KCAL:^ 41.2f}"
                )
            print("".ljust(PLENGTH, "-"))
            break
    # END SD Gsolv

    # calculate average G correction
    print("\nCalculating Boltzmann averaged free energy of ensemble!\n")
    avGcorrection = {
        "avGcorrection": {},
        "avG": {},
        "avE": {},
        "avGsolv": {},
        "avGRRHO": {},
    }
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
                # f"{'avGcorrection(T) /a.u.':>22}"
            )
        elif not config.evaluate_rrho:
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGsolv(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
                # f"{'avGcorrection(T) /a.u.':>22}"
            )
        elif config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGmRRHO(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
                # f"{'avGcorrection(T) /a.u.':>22}"
            )
    else:
        line = (
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
            f"{'avGmRRHO(T) /a.u.':>16} {'avGsolv(T) /a.u.':>16} "
            f"{'avG(T) /a.u.':>14}"
            # f" {'avGcorrection(T) /a.u.':>22}"
        )
    print(line)
    print("".ljust(int(PLENGTH), "-"))
    for temperature in trange:
        # get free energy at (T)
        for conf in calculate:
            if not config.evaluate_rrho:
                rrho = None
            else:
                rrho = "lowlevel_grrho_info"
            if config.solvent == "gas":
                solv = None
            else:
                solv = "lowlevel_gsolv_info"
            e = "lowlevel_sp_info"
            conf.calc_free_energy(e=e, solv=solv, rrho=rrho, t=temperature)
        try:
            minfreeT = min(
                [conf.free_energy for conf in calculate if conf.free_energy is not None]
            )
        except ValueError:
            raise ValueError

        calculate = calc_boltzmannweights(calculate, "free_energy", temperature)
        avG = 0.0
        avE = 0.0
        avGRRHO = 0.0
        # avHRRHO but with new boltzmann weights?
        avGsolv = 0.0
        for conf in calculate:
            avG += conf.bm_weight * conf.free_energy
            avE += conf.bm_weight * conf.lowlevel_sp_info["energy"]
            avGRRHO += conf.bm_weight * conf.lowlevel_grrho_info["range"].get(
                temperature, 0.0
            )
            avGsolv += conf.bm_weight * conf.lowlevel_gsolv_info["range"].get(
                temperature, 0.0
            )

        avGcorrection["avG"][temperature] = avG
        avGcorrection["avE"][temperature] = avE
        avGcorrection["avGRRHO"][temperature] = avGRRHO
        avGcorrection["avGsolv"][temperature] = avGsolv
        for conf in calculate:
            if conf.free_energy == minfreeT:
                avGcorrection["avGcorrection"][temperature] = avG - conf.free_energy
        # printout:
        if not config.evaluate_rrho or config.solvent == "gas":
            if not config.evaluate_rrho and config.solvent == "gas":
                line = (
                    f"{temperature:^15} {avE:>14.7f}  {avG:>14.7f} "
                    # f"{avGcorrection['avGcorrection'][temperature]:>22.7f}"
                )
            elif not config.evaluate_rrho:
                line = (
                    f"{temperature:^15} {avE:>14.7f} {avGsolv:>16.7f} "
                    f"{avG:>14.7f} "
                    # f"{ avGcorrection['avGcorrection'][temperature]:>22.7f}"
                )
            elif config.solvent == "gas":
                line = (
                    f"{temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                    f"{avG:>14.7f} "
                    # f"{ avGcorrection['avGcorrection'][temperature]:>22.7f}"
                )
        else:
            line = (
                f"{temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                f"{avGsolv:>16.7f} {avG:>14.7f} "
                # f"{ avGcorrection['avGcorrection'][temperature]:>22.7f}"
            )
        if temperature == config.temperature:
            print(line, "    <<==part2==")
        else:
            print(line)
    print("".ljust(int(PLENGTH), "-"))
    print("")

    # reset boltzmannweights to correct temperature
    # get free energy at (T)
    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "lowlevel_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "lowlevel_gsolv_info"
        e = "lowlevel_sp_info"
        conf.calc_free_energy(e=e, solv=solv, rrho=rrho)
    # ensembledata is used to store avGcorrection
    ensembledata.comment = [
        lowestconf,
        "storage for avGcorrection of ensemble",
        f"corresponding to CONF{lowestconf}",
    ]
    ensembledata.avGcorrection = avGcorrection

    if not calculate:
        print_errors(f"{'ERROR:':{WARNLEN}}No conformers left!", save_errors)
        print("Going to exit!")
        sys.exit(1)

    calculate = calc_boltzmannweights(calculate, "free_energy", config.temperature)

    try:
        minfree = min(
            [conf.free_energy for conf in calculate if conf.free_energy is not None]
        )
    except ValueError:
        raise ValueError
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL

    for conf in calculate:
        if conf.free_energy == minfree:
            ensembledata.bestconf["part2"] = conf.id

    #
    print("")
    onlyprintout = deepcopy(calculate)
    onlyprintout.sort(reverse=True, key=lambda x: float(x.bm_weight))
    for i in (100, 95, 90, 80, 70):
        sumup = 0.0
        for conf in list(onlyprintout):
            sumup += conf.bm_weight
            if sumup > (i / 100):
                if conf.bm_weight <= (1 - (i / 100)):
                    onlyprintout.pop(onlyprintout.index(conf))

        # write ensemble
        outfile = f"enso_ensemble_part2_p_{i}.xyz"
        # move_recursively(config.cwd, outfile)
        kwargs = {"energy": "xtb_energy", "rrho": "lowlevel_grrho_info"}
        write_trj(
            sorted(onlyprintout, key=lambda x: float(x.free_energy)),
            config.cwd,
            outfile,
            config.func,
            config.nat,
            "free_energy",
            overwrite=True,
            **kwargs,
        )

    # SORTING for the next part:
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("Conformers considered further".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    # evaluate conformer consideration based on Boltzmann-population
    calculate.sort(reverse=True, key=lambda x: float(x.bm_weight))
    sumup = 0.0
    for conf in list(calculate):
        sumup += conf.bm_weight
        if sumup >= (config.part2_threshold / 100):
            if conf.bm_weight < (1 - (config.part2_threshold / 100)):
                mol = calculate.pop(calculate.index(conf))
                mol.part_info["part2"] = "refused"
                store_confs.append(mol)
            else:
                conf.part_info["part2"] = "passed"
        else:
            conf.part_info["part2"] = "passed"

    ensembledata.nconfs_per_part["part2"] = len(calculate)
    if calculate:
        print(
            f"\nConformers that are below the Boltzmann threshold G_thr(2) "
            f"of {config.part2_threshold}%:"
        )
        print_block(["CONF" + str(i.id) for i in calculate])
    else:
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
    # write ensemble
    move_recursively(config.cwd, "enso_ensemble_part2.xyz")
    if config.evaluate_rrho:
        kwargs = {"energy": "xtb_energy_unbiased", "rrho": "lowlevel_grrho_info"}
    else:
        kwargs = {"energy": "xtb_energy_unbiased"}
    write_trj(
        sorted(calculate, key=lambda x: float(x.free_energy)),
        config.cwd,
        "enso_ensemble_part2.xyz",
        config.func,
        config.nat,
        "free_energy",
        **kwargs,
    )

    # write coord.enso_best
    for conf in calculate:
        if conf.id == ensembledata.bestconf["part2"]:
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
                        conf.free_energy, conf.lowlevel_grrho_info["energy"], conf.id
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
        print(
            "\n***---------------------------------------------------------***"
        )
        print(
            "Printing most relevant errors again, just for user convenience:"
        )
        for _ in list(save_errors):
            print(save_errors.pop())
        print(
            "***---------------------------------------------------------***"
        )

    tmp = int((PLENGTH - len("END of Part2")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part2" + "".rjust(tmp, "<"))
    if config.progress:
        print("#>>># CENSO: Finished part2", file=sys.stderr)
    return config, calculate, store_confs, ensembledata
