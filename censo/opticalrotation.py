"""
module for the calculation of optical rotation 
"""

import os
import shutil
import sys
from random import normalvariate
from multiprocessing import JoinableQueue as Queue
from .cfg import PLENGTH, DIGILEN, AU2KCAL, WARNLEN, dfa_settings
from .parallel import run_in_parallel
from .orca_job import OrcaJob
from .tm_job import TmJob
from .utilities import (
    calc_boltzmannweights,
    printout,
    print_block,
    new_folders,
    last_folders,
    print,
    print_errors,
    calc_std_dev,
    ensemble2coord,
)


def part5(config, conformers, store_confs, ensembledata):
    """
    Calculate optical rotation on the populated
    conformers (either directly from part2 OPTIMIZATION or after REFINEMENT
    (part3))
    """
    save_errors = []
    print("\n" + "".ljust(PLENGTH, "-"))
    print("OPTICAL ROTATION MODE - PART5".center(PLENGTH, " "))
    if config.progress:
        print("#>>># CENSO: Starting part5", file=sys.stderr)
    print("".ljust(PLENGTH, "-") + "\n")
    # print flags for part5
    info = []
    info.append(["optical_rotation", "Part5"])
    info.append(["freq_or", "frequency in [nm]"])
    info.append(["func_or_scf", "functional for SCF"])
    info.append(["func_or", "functional for optical rotation"])
    info.append(["basis_or", "basis set for optical rotation"])
    if config.part3:
        info.append(["part3_threshold", "Boltzmann sum threshold employed"])
    elif config.part2:
        info.append(["part2_P_threshold", "Boltzmann sum threshold employed"])
    elif config.part1:
        info.append(["part2_P_threshold", "Boltzmann sum threshold employed"])
    if config.solvent != "gas":
        info.append(["solvent", "solvent"])
        if config.prog == "tm":
            info.append(["printoption", "solvation model", "cosmo"])

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
    # # end print

    calculate = []  # has to be calculated in this run
    prev_calculated = []  # was already calculated in a previous run
    try:
        store_confs
    except NameError:
        store_confs = []  # stores all confs which are sorted out!

    # setup queues
    q = Queue()
    resultq = Queue()

    unoptimized_warning = False
    # sort conformers:
    for conf in list(conformers):
        if conf.removed:
            store_confs.append(conformers.pop(conformers.index(conf)))
            print(f"CONF{conf.id} is removed as requested by the user!")
            continue
        if (
            conf.part_info["part2"] != "passed"
            and conf.optimization_info["info"] != "calculated"
        ):
            unoptimized_warning = True

        if config.part3:
            # calc boltzmann weights from part3
            energy = "highlevel_sp_info"
            rrho = "highlevel_grrho_info"
            gsolv = "highlevel_gsolv_info"
            boltzmannthr = config.part3_threshold
        elif config.part2:
            # part3 is not calculated use boltzmann weights directly from part2
            energy = "lowlevel_sp_info"
            rrho = "lowlevel_grrho_info"
            gsolv = "lowlevel_gsolv_info"
            boltzmannthr = config.part2_P_threshold
        elif config.part1:
            # part2 is not calculated use boltzmann weights directly from part1
            # --> misappropriate config.part2_P_threshold
            # This means starting from not DFT optimized geometries!
            energy = "prescreening_sp_info"
            rrho = "prescreening_grrho_info"
            gsolv = "prescreening_gsolv_info"
            boltzmannthr = config.part2_P_threshold
        else:
            print("UNEXPECTED BEHAVIOUR")
        mol = conformers.pop(conformers.index(conf))
        if getattr(conf, energy)["info"] != "calculated":
            store_confs.append(mol)
            continue
        elif getattr(conf, rrho)["info"] != "calculated" and config.evaluate_rrho:
            store_confs.append(mol)
            continue
        elif getattr(conf, gsolv)["info"] != "calculated" and config.solvent != "gas":
            store_confs.append(mol)
            continue
        else:
            calculate.append(mol)

    if unoptimized_warning:
        print_errors(
            f"{'INFORMATION:':{WARNLEN}}Conformers have not been optimized at DFT level!!!\n"
            f"{'':{WARNLEN}}Use results with care!\n",
            save_errors,
        )
    # SI geometry:
    if unoptimized_warning:
        ensembledata.si["part5"]["Geometry"] = "GFNn-xTB (input geometry)"
    else:
        ensembledata.si["part5"]["Geometry"], _ = config.get_method_name(
            "xtbopt",
            func=config.func,
            basis=config.basis,
            solvent=config.solvent,
            sm=config.sm2,
        )
        ensembledata.si["part5"]["Geometry"] += f" @optlevel: {config.optlevel2}"
    #
    if not calculate and not prev_calculated:
        print(f"{'ERROR:':{WARNLEN}}No conformers left!")
        print("Going to exit!")
        sys.exit(1)

    calculate.sort(key=lambda x: int(x.id))
    print(f"Considering the following {len(calculate)} conformers:")
    print_block(["CONF" + str(i.id) for i in calculate])

    # Calculate Boltzmann weight for confs:
    if config.part3:
        using_part = "part3 - refinement"
        if not config.evaluate_rrho:
            rrho = None
            rrho_method = None
        else:
            rrho_method, _ = config.get_method_name(
                "rrhoxtb",
                bhess=config.bhess,
                gfn_version=config.part3_gfnv,
                sm=config.sm_rrho,
                solvent=config.solvent,
            )
        if config.solvent == "gas":
            gsolv = None
            energy_method, _ = config.get_method_name(
                "xtbopt",
                func=config.func3,
                basis=config.basis3,
                sm=config.smgsolv3,
                gfn_version=config.part3_gfnv,
                solvent=config.solvent,
            )
        else:
            if config.smgsolv3 in ("cosmors", "cosmors-fine"):
                tmp_name = "cosmors"
            elif config.smgsolv3 in ("alpb_gsolv", "gbsa_gsolv", "smd_gsolv"):
                tmp_name = config.smgsolv3
            else:
                tmp_name = "sp_implicit"
            energy_method, solv_method = config.get_method_name(
                tmp_name,
                func=config.func3,
                basis=config.basis3,
                sm=config.smgsolv3,
                gfn_version=config.part3_gfnv,
                solvent=config.solvent,
            )
    elif config.part2:
        using_part = "part2 - optimization"
        if not config.evaluate_rrho:
            rrho = None
            rrho_method = None
        else:
            rrho_method, _ = config.get_method_name(
                "rrhoxtb",
                bhess=config.bhess,
                gfn_version=config.part2_gfnv,
                sm=config.sm_rrho,
                solvent=config.solvent,
            )
        if config.solvent == "gas":
            gsolv = None
            energy_method, _ = config.get_method_name(
                "xtbopt",
                func=config.func,
                basis=config.basis,
                sm=config.smgsolv2,
                gfn_version=config.part2_gfnv,
                solvent=config.solvent,
            )
        else:
            if config.smgsolv2 in ("cosmors", "cosmors-fine"):
                tmp_name = "cosmors"
            elif config.smgsolv2 in ("alpb_gsolv", "gbsa_gsolv", "smd_gsolv"):
                tmp_name = config.smgsolv2
            else:
                tmp_name = "sp_implicit"
            energy_method, solv_method = config.get_method_name(
                tmp_name,
                func=config.func,
                basis=config.basis,
                sm=config.smgsolv2,
                gfn_version=config.part2_gfnv,
                solvent=config.solvent,
            )
    elif config.part1:
        using_part = "part1 - prescreening"
        # on DFT unoptimized geometries!
        if not config.evaluate_rrho:
            rrho = None
            rrho_method = None
        else:
            rrho_method, _ = config.get_method_name(
                "rrhoxtb",
                bhess=config.bhess,
                gfn_version=config.part1_gfnv,
                sm=config.sm_rrho,
                solvent=config.solvent,
            )
        if config.solvent == "gas":
            gsolv = None
            energy_method, _ = config.get_method_name(
                "xtbopt",
                func=config.func,
                basis=config.basis,
                sm=config.smgsolv1,
                gfn_version=config.part1_gfnv,
                solvent=config.solvent,
            )
        else:
            if config.smgsolv1 in ("cosmors", "cosmors-fine"):
                tmp_name = "cosmors"
            elif config.smgsolv2 in ("alpb_gsolv", "gbsa_gsolv", "smd_gsolv"):
                tmp_name = config.smgsolv1
            else:
                tmp_name = "sp_implicit"
            energy_method, solv_method = config.get_method_name(
                tmp_name,
                func=config.func,
                basis=config.basis,
                sm=config.smgsolv1,
                gfn_version=config.part1_gfnv,
                solvent=config.solvent,
            )

    # SI information:-------------------------------------------------------
    ensembledata.si["part5"]["Energy"] = f"using Energy from {using_part}"
    ensembledata.si["part5"]["Energy_settings"] = f"see {using_part}"
    ensembledata.si["part5"]["G_mRRHO"] = f"using G_mRRHO from {using_part}"
    ensembledata.si["part5"]["G_solv"] = f"using G_solv from {using_part}"
    ensembledata.si["part5"]["Threshold"] = f"Boltzmann sum threshold: {boltzmannthr} %"
    ensembledata.si["part5"]["main QM code"] = str("tm").upper()
    # END SI information--------------------------------------------------------

    for conf in calculate:
        conf.calc_free_energy(
            e=energy,
            solv=gsolv,
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

    # printout for part5 -------------------------------------------------------
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("* Gibbs free energies used in part5 *".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, energy)["energy"],
        # lambda conf: getattr(conf, gsolv)["energy"],
        lambda conf: getattr(conf, gsolv)
        .get("range", {})
        .get(config.temperature, getattr(conf, gsolv, {"energy": 0.0})["energy"]),
        # lambda conf: getattr(conf, rrho)["energy"],
        lambda conf: conf.get_mrrho(config.temperature, rrho, config.consider_sym),
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "rel_free_energy"),
        lambda conf: getattr(conf, "bm_weight") * 100,
    ]
    columnheader = [
        "CONF#",
        "E [Eh]",
        "Gsolv [Eh]",
        "GmRRHO [Eh]",
        "Gtot",
        "ΔGtot",
        "Boltzmannweight",
    ]
    columndescription = [
        "",
        "",
        "",
        "",
        "[Eh]",
        "[kcal/mol]",
        f"  % at {config.temperature:.2f} K",
    ]
    columnformat = ["", (12, 7), (12, 7), (12, 7), (12, 7), (5, 2), (5, 2)]
    columndescription[1] = energy_method
    if config.solvent != "gas":
        columndescription[2] = solv_method
    columndescription[3] = rrho_method

    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho:
            # ignore rrho in printout
            columncall.pop(3)
            columnheader.pop(3)
            columndescription.pop(3)
            columnformat.pop(3)
        if config.solvent == "gas":
            columncall.pop(2)
            columnheader.pop(2)
            columndescription.pop(2)
            columnformat.pop(2)

    printout(
        os.path.join(config.cwd, "part5.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
    )
    calculate.sort(reverse=True, key=lambda x: float(x.bm_weight))
    sumup = 0.0
    for conf in list(calculate):
        sumup += conf.bm_weight
        if sumup >= boltzmannthr:
            if conf.bm_weight < (1 - boltzmannthr):
                store_confs.append(calculate.pop(calculate.index(conf)))
    print(f"\nConformers that are below the Boltzmann-thr of {boltzmannthr}:")
    ensembledata.nconfs_per_part["part5"] = len(calculate)
    print_block(["CONF" + str(i.id) for i in calculate])

    # create folder
    folder = "OR"
    save_errors, store_confs, calculate = new_folders(
        config.cwd, calculate, folder, save_errors, store_confs
    )
    if config.part3 and config.part2 or config.part2:
        # need to copy optimized coord to folder
        for conf in list(calculate):
            tmp1 = os.path.join(config.cwd, "CONF" + str(conf.id), config.func, "coord")
            tmp2 = os.path.join("CONF" + str(conf.id), folder, "coord")
            try:
                shutil.copy(tmp1, tmp2)
            except FileNotFoundError:
                print(f"{'ERROR:':{WARNLEN}}can't copy optimized geometry!")
                store_confs.append(calculate.pop(calculate.index(conf)))
    elif config.part3:
        # structures can be DFT optimized or not (part2 might not have been run)
        for conf in list(calculate):
            tmp1 = os.path.join(config.cwd, "CONF" + str(conf.id), "part3", "coord")
            tmp2 = os.path.join("CONF" + str(conf.id), folder, "coord")
            try:
                shutil.copy(tmp1, tmp2)
            except FileNotFoundError:
                print(f"{'ERROR:':{WARNLEN}}can't copy geometry!")
                store_confs.append(calculate.pop(calculate.index(conf)))
    elif config.part1:
        # do not use coord from folder config.func it could be optimized if
        # part2 has ever been run, take coord from ensemble file
        # write coord to folder
        calculate, store_confs, save_errors = ensemble2coord(
            config, folder, calculate, store_confs, save_errors
        )

    if not calculate:
        print(f"{'ERROR:':{WARNLEN}}No conformers left!")
        print("Going to exit!")
        sys.exit(1)

    # check if OR calculated before!
    for conf in list(calculate):
        if getattr(conf, "optical_rotation_info")["info"] == "calculated":
            prev_calculated.append(calculate.pop(calculate.index(conf)))
        elif getattr(conf, "optical_rotation_info")["info"] == "failed":
            print(
                f"{'INFORMATION:':{WARNLEN}}The calculation failed for CONF{conf.id} in the previous run."
            )
            store_confs.append(calculate.pop(calculate.index(conf)))

    instruction_or = {
        "jobtype": "opt-rot_sp",  # opt-rot only escf ; opt-rot_sp SP then escf
        "func": config.func_or_scf,
        "func2": config.func_or,
        "basis": getattr(
            config,
            "basis_or",
            dfa_settings.composite_method_basis.get(config.func_or_scf, "def2-mTZVPP"),
        ),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": "gas",
        "sm": "cosmo",
        "success": False,
        "freq_or": config.freq_or,
    }
    if config.prog == "orca":
        print(
            f"{'ERROR:':{WARNLEN}}Can't calculate OR with ORCA! You can use TM instead."
        )
        print("Going to exit!")
        sys.exit(1)
        # ORCA can't calculate optical rotation!!! -->
        job = OrcaJob
        if config.solvent != "gas":
            instruction_or["solvent"] = config.solvent
            instruction_or["sm"] = "cpcm"
    if config.prog == "tm":
        job = TmJob
        instruction_or["prepinfo"] = ["clear", "-grid", "2", "-scfconv", "6"]
        if config.basis == config.basis_or:
            instruction_or["copymos"] = config.func
            instruction_or["jobtype"] = "opt-rot"
        if config.solvent != "gas":
            instruction_or["solvent"] = config.solvent
            instruction_or["sm"] = "cosmo"

    instruction_or["method"], _ = config.get_method_name(
        instruction_or["jobtype"],
        func=instruction_or["func2"],
        basis=instruction_or["basis"],
        solvent=instruction_or["solvent"],
        prog=config.prog,
        func2=instruction_or["func"],
        sm=instruction_or["sm"],
    )
    # SI update-----------------------------------------------------------------
    ensembledata.si["part5"]["Optical rotation"] = (
        instruction_or["method"]
        + f" using {' '.join(instruction_or['prepinfo']).replace('-', '').replace('clear', '')}"
    )
    # END SI update--------------------------------------------------------------

    print(f"\nOptical-rotation is calculated at {instruction_or['method']} level.\n")
    check = {True: "was successful", False: "FAILED"}
    pl = config.lenconfx + 4 + len(str("/" + folder))

    if calculate:
        calculate = run_in_parallel(
            config,
            q,
            resultq,
            job,
            config.maxthreads,
            config.omp,
            calculate,
            instruction_or,
            config.balance,
            folder,
        )
        calculate.sort(key=lambda x: int(x.id))
        try:
            max_fmt = max(
                [
                    len(str(item.job["erange1"].get(config.freq_or[0])).split(".")[0])
                    for item in calculate
                ]
            )
            max_fmt += 9
        except Exception:
            max_fmt = 16
        for conf in list(calculate):
            line = (
                f"Optical-rotation calculation {check[conf.job['success']]}"
                f" for {last_folders(conf.job['workdir'], 2):>{pl}} at {config.freq_or[0]} nm: "
                f"{conf.job['erange1'].get(config.freq_or[0], 0.00):> {max_fmt}.7f}"
                f"   populated to {conf.bm_weight*100:.2f} %"
            )
            print(line)
            if not conf.job["success"]:
                save_errors.append(line)
                conf.optical_rotation_info["info"] = "failed"
                conf.optical_rotation_info["method"] = instruction_or["method"]
                conf.part_info["part5"] = "refused"
                store_confs.append(calculate.pop(calculate.index(conf)))
            else:
                conf.optical_rotation_info["range"] = conf.job["erange1"]
                conf.optical_rotation_info["info"] = "calculated"
                conf.optical_rotation_info["method"] = instruction_or["method"]
                conf.part_info["part5"] = "passed"

    if prev_calculated:
        try:
            max_fmt = max(
                [
                    len(
                        str(
                            item.optical_rotation_info["range"].get(config.freq_or[0])
                        ).split(".")[0]
                    )
                    for item in prev_calculated
                ]
            )
            max_fmt += 9
        except Exception as e:
            print(e)
            max_fmt = 16
        prev_calculated.sort(key=lambda x: int(x.id))
        for conf in list(prev_calculated):
            conf.job["workdir"] = os.path.normpath(
                os.path.join(config.cwd, "CONF" + str(conf.id), folder)
            )
            print(
                f"Optical-rotation calculation {check[True]} for "
                f"{last_folders(conf.job['workdir'], 2):>{pl}} at {config.freq_or[0]} nm: "
                f"{conf.optical_rotation_info['range'].get(config.freq_or[0], 0.0000):> {max_fmt}.7f}"
                f"  populated to {conf.bm_weight*100:.2f} % "
            )
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))

    outlen = config.lenconfx + 1
    if len("ENSEMBLE") >= config.lenconfx:
        outlen = len("ENSEMBLE") + 1

    for freq in config.freq_or:
        averaged_or = 0.0
        with open(
            os.path.join(config.cwd, f"OR_{int(freq)}.dat"), "w", newline=None
        ) as outdata:
            outdata.write(
                f"{'#label':{outlen}} {'#unmod_alpha':>{max_fmt}} {'#%pop':^7} {'#pop_alpha':{max_fmt}} \n"
            )
            for conf in calculate:
                averaged_or += conf.bm_weight * conf.optical_rotation_info["range"].get(
                    freq, 0.0
                )
                # CONFX
                outdata.write(
                    f"{'CONF'+str(conf.id):{outlen}} "
                    + f"{conf.optical_rotation_info['range'].get(freq, 0.0):> {max_fmt}.7f} "
                    + f"{conf.bm_weight*100:> 6.2f} "
                    + f"{conf.optical_rotation_info['range'].get(freq, 0.0)*conf.bm_weight:> {max_fmt}.7f}\n"
                )
            # ENSEMBLE
            outdata.write(
                f"{'ENSEMBLE':{outlen}} "
                + f"{'-':^{max_fmt-1}} "
                + f"{100.00:> 6.2f} "
                + f"{averaged_or:> {max_fmt}.7f}\n"
            )
        print(
            f"\nAveraged specific rotation at {freq} nm : "
            f"{averaged_or:> {max_fmt}.3f}   in deg*[dm(g/cc)]^(-1)"
        )

    if (
        all(
            [
                conf.lowlevel_gsolv_compare_info["std_dev"] is not None
                for conf in calculate
            ]
        )
        and calculate
    ):
        for freq in config.freq_or:
            all_or = []
            for _ in range(1000):
                averaged_or = 0.0
                for conf in calculate:
                    conf.calc_free_energy(
                        e=energy,
                        solv=gsolv,
                        rrho=rrho,
                        t=config.temperature,
                        consider_sym=config.consider_sym,
                    )
                    conf.free_energy += normalvariate(
                        0.0, conf.lowlevel_gsolv_compare_info["std_dev"]
                    )
                calculate = calc_boltzmannweights(
                    calculate, "free_energy", config.temperature
                )
                for conf in calculate:
                    averaged_or += conf.bm_weight * conf.optical_rotation_info[
                        "range"
                    ].get(freq)
                all_or.append(averaged_or)
            try:
                max_fmt = max(
                    [
                        len(
                            str(
                                item.optical_rotation_info["range"].get(
                                    config.freq_or[0]
                                )
                            ).split(".")[0]
                        )
                        for item in calculate
                    ]
                )
                max_fmt += 9
            except Exception as e:
                print(e)
                max_fmt = 16
            print(
                f"    SD based on SD of Gsolv (part2)    "
                f": {calc_std_dev(all_or):> {max_fmt}.3f}   in deg*[dm(g/cc)]^(-1)"
            )

    if calculate:
        for freq in config.freq_or:
            all_or = []
            for _ in range(1000):
                averaged_or = 0.0
                for conf in calculate:
                    conf.calc_free_energy(
                        e=energy,
                        solv=gsolv,
                        rrho=rrho,
                        t=config.consider_sym,
                        consider_sym=config.consider_sym,
                    )
                    conf.free_energy += normalvariate(0.0, (0.4 / AU2KCAL))
                calculate = calc_boltzmannweights(
                    calculate, "free_energy", config.temperature
                )
                for conf in calculate:
                    averaged_or += conf.bm_weight * conf.optical_rotation_info[
                        "range"
                    ].get(freq)
                all_or.append(averaged_or)
            try:
                max_fmt = max(
                    [
                        len(
                            str(
                                item.optical_rotation_info["range"].get(
                                    config.freq_or[0]
                                )
                            ).split(".")[0]
                        )
                        for item in calculate
                    ]
                )
                max_fmt += 9
            except Exception as e:
                print(e)
                max_fmt = 16
            print(
                f"    SD based on SD in G of 0.4 kcal/mol"
                f": {calc_std_dev(all_or):> {max_fmt}.3f}   in deg*[dm(g/cc)]^(-1)"
            )
    if calculate:
        # print min and max OR values and corresponding conformers, but no Boltzmann weighting
        print("")
        output = []
        output.append(f"{'frequency'}     {'min(OR)'} {'CONF'}    {'max(OR)'} {'CONF'}")
        output.append("".ljust(int(50), "-"))
        for freq in config.freq_or:
            tmp_or = {}
            for conf in calculate:
                tmp_or[conf.id] = conf.optical_rotation_info["range"].get(freq)
            minx = min(list(tmp_or.values()))
            minid = next(key for key, value in tmp_or.items() if value == minx)
            maxx = max(list(tmp_or.values()))
            maxid = next(key for key, value in tmp_or.items() if value == maxx)

            output.append(
                f"{float(freq): ^12}  {minx: .2f} {'CONF'+str(minid)}    {maxx: .2f} {'CONF'+str(maxid)}"
            )
        output.append("".ljust(int(50), "-"))
        output.append("* min and max values are not Boltzmann weighted.")
        for line in output:
            print(line)

    for conf in calculate:
        conf.reset_job_info()

    # save current data to jsonfile
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
    if not calculate:
        print("ERROR: No conformers left!")
        print("Going to exit!")
        sys.exit(1)

    # end printout for part5
    tmp = int((PLENGTH - len("END of Part5")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part5" + "".rjust(tmp, "<"))
    if config.progress:
        print("#>>># CENSO: Finished part5", file=sys.stderr)
    return config, calculate, store_confs, ensembledata
