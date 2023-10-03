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
from .cfg import PLENGTH, CODING, AU2KCAL, DIGILEN, WARNLEN, qm_prepinfo, dfa_settings
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
    conf_in_interval,
)
from .orca_processor import OrcaProc
from .tm_processor import TmProc
from .parallel import run_in_parallel

from censo.core import CensoCore
from censo.settings import CensoSettings
from censo.parallel import ProcessHandler


class Optimization(Part):
    
    alt_name = "part2"

    def __init__(self, core: CensoCore, settings: CensoSettings):
        super().__init__(core, settings, "optimization")


    @timeit
    def run(self) -> None:
        """
        Optimization of the ensemble at DFT level (possibly with implicit solvation)

        Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
        by calculating the thermodynamics of the optimized geometries

        Alternatively just run the complete optimization for every conformer with xtb as driver (decide with 'opt_spearman')

        TODO - implement regular optimization (no xtb driver)
        """

        """
        don't use this for now, use new keyword 'maxcyc' instead
        if config.nat > 200:
            # stopcycle = don't optimize more than stopcycle cycles
            stopcycle = config.nat * 2
        else:
            stopcycle = 200
        """

        handler = ProcessHandler(self._instructions)

        # set folder to do the calculations in
        folder = os.path.join(self.core.workdir, self.__class__.__name__.lower())
        if os.path.isdir(folder):
            # TODO - warning? stderr?
            print(f"Folder {folder} already exists. Potentially overwriting files.")
        elif os.system(f"mkdir {folder}") != 0 and not os.path.isdir(folder):
            # TODO - error handling stderr case? is this reasonable?
            raise RuntimeError(f"Could not create directory for {self.__class__.__name__.lower()}")

        if self._instructions["opt_spearman"]:
            """
            optimization using spearman threshold, updated every 'optcycles' steps
            """
            # TODO - run single-point for every conformer first in order to set preliminary ordering and fuzzy threshold? or take from screening

            # make a separate list of conformers that only includes (considered) conformers that are not converged
            # at this point it's just self.core.conformers
            confs_nc = self.core.conformers

            stopcond_converged = False
            ncyc = 0
            while (
                not stopcond_converged 
                and ncyc < self._instructions["maxcyc"]
            ):
                # basically TL;DR:
                # while the number of cycles is less than stopcycle, optimize for optcycles steps
                # get resulting geometries and energies and update conformers with this
                # update ensemble

                handler.conformers = confs_nc

                # update number of cycles
                ncyc += self._instructions["optcycles"]

                # run optimizations for 'optcycles' steps
                results_opt = handler.execute(["xtb_opt"], os.path.join(folder, "opt"))

                # TODO - crestcheck if ncyc >= 2

                # run xtb_rrho for finite temperature contributions
                # for now only after the first 'optcycles' steps or after at least 6 cycles are done
                # TODO - make this better
                if ncyc >= 6:
                    """run xtb_rrho"""
                    results_rrho = handler.execute(["xtb_rrho"], os.path.join(folder, "rrho"))

                # TODO - recalculate fuzzy threshold

                # kick out conformers above threshold
                # TODO - do this for confs_nc and self.core.conformers
                threshold = self._instructions.get("threshold", None)
                if not threshold is None:
                    # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
                    limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
                    
                    # filter out all conformers below threshold
                    # so that 'filtered' contains all conformers that should not be considered any further
                    filtered = [
                        conf for conf in filter(
                            lambda x: self.key(x) > limit + threshold, 
                            self.core.conformers
                        )
                    ]
                    
                    # update the conformer list in core (remove conf if below threshold)
                    self.core.update_conformers(filtered)  

                    # also remove conformers from confs_nc
                    for conf in filtered:
                        confs_nc.remove(conf)
                else:
                    """
                    TODO
                    print warning that no threshold could be determined
                    (should not happen but doesn't necessarily brake the program)
                    """
                    print("...")

                # update list of converged conformers
                for conf in confs_nc:
                    if conf.results[self.__class__.__name__.lower()]["xtb_opt"]["converged"]:
                        confs_nc.remove(conf)

                # check if all (considered) conformers converged - End of While-loop
                stopcond_converged = len(confs_nc) == 0


    def key(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy and Gmrrho
        """


    def write_results(self) -> None:
        """
        formatted write of part results (optional)
        """
        pass


def part2(config, conformers, store_confs, ensembledata):
                maxecyc = max([len(conf.job["ecyc"]) for conf in calculate])
                print(f"Max number of performed iterations: {maxecyc}")
                if len(calculate + prev_calculated) == 1:
                    # can't do spearman with only one conf
                    run_spearman = False
                elif len(calculate) > 1 and run > 1:
                    run_spearman = True
                else:
                    run_spearman = True
                if run == 1:
                    # only evaluate spearman starting from second cycle
                    print("Spearman rank evaluation is performed in the next cycle.")
                    cycle_spearman.append("")
                    run_spearman = False

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
                    # rrho_energy = "energy_rrho"
                    rrho_energy = "rrho_optimization"
                else:
                    rrho_energy = None  # to get 0.0 contribution
                for i in range(maxecyc):
                    try:
                        minecyc.append(
                            min(
                                [
                                    conf.job["ecyc"][i]
                                    + conf.get_mrrho(
                                        config.fixed_temperature,
                                        rrho=rrho_energy,
                                        consider_sym=config.consider_sym,
                                    )
                                    # + getattr(conf, "optimization_info").get(
                                    #    rrho_energy, 0.0
                                    # )
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
                                + conf.get_mrrho(
                                    config.fixed_temperature,
                                    rrho=rrho_energy,
                                    consider_sym=config.consider_sym,
                                )
                                # + getattr(conf, "optimization_info").get(
                                #    rrho_energy, 0.0
                                # )
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
                    cycle_spearman.append(f"{sum(evalspearman)/2:.3f}")
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
                            f"CONF{conf.id} is above threshold, don't optimize "
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
                            f"{'CONF' + str(conf.id):{config.lenconfx+4}} within {conf.job['cycles']:>3} cycles"
                        )
                        conf.optimization_info["info"] = "calculated"
                        conf.optimization_info["energy"] = conf.job["energy"]
                        conf.optimization_info["cycles"] = conf.job["cycles"]
                        conf.optimization_info["ecyc"] = conf.job["ecyc"]
                        conf.optimization_info["decyc"] = conf.job["decyc"]
                        conf.optimization_info[
                            "convergence"
                        ] = (
                            "converged"
                        )  #### THIS IS NOT CORRECT /or can lead to errors!
                        conf.optimization_info["method"] = instruction_opt["method"]
                        prev_calculated.append(calculate.pop(calculate.index(conf)))
                #
                maxecyc_prev = maxecyc
                run += 1
                if (run >= 3 or
                    ( (nconf_start -len(calculate)) >= 5 )):
                    print(f"{'INFORMATION:':{WARNLEN}}Saving enso.json to update optimization information.")
                    # save current data to jsonfile
                    # safety measure for runs with thousands of conformers to keep
                    # information on optimized structures and make restart easier
                    # in case of a crash or cluster stop in part2
                    config.write_json(
                        config.cwd,
                        [i.provide_runinfo() for i in calculate]
                        + [i.provide_runinfo() for i in prev_calculated]
                        + [i.provide_runinfo() for i in store_confs]
                        + [ensembledata],
                        config.provide_runinfo(),
                    )
            # END while loop
        else:
            # use standard optimization!
            # update instruct_opt
            tic = time.perf_counter()
            del instruction_opt["optcycles"]
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
        if config.opt_spearman and (len(timings) == len(nconf_cycle)):
            try:
                tl = max([len(f"{i: .2f}") for i in timings])
                if tl > 7:
                    tmp1 = tl
                else:
                    tmp1 = 7
                print("Timings:")
                print(f"Cycle:  [s]  {'#nconfs':^{tmp1}}  Spearman coeff.")
                for i in range(len(timings)):
                    try:
                        print(
                            f"{i+1:>4}  {timings[i]:> {tl}.2f}  {nconf_cycle[i]:^{tmp1}}  {cycle_spearman[i]}"
                        )
                    except Exception as e:
                        print(i, len(nconf_cycle), len(cycle_spearman))
                        print(f"{i+1:>4}  {timings[i]:> {tl}.2f}  {'-':^{tmp1}}  {'-'}")
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
    end_opt_time = time.perf_counter()
    ensembledata.part_info["part2_opt"] = end_opt_time - start_opt_time
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

    if config.prog == "tm":
        job = TmProc
    elif config.prog == "orca":
        job = OrcaProc
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
        "basis": getattr(
            config,
            "basis",
            dfa_settings.composite_method_basis.get(config.func, "def2-TZVP"),
        ),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": config.smgsolv2,
        "temperature": config.temperature,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
        "gfn_version": config.part2_gfnv,
        "onlyread": config.onlyread,
    }
    tmp_SI = None
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
                job = TmProc
                exc_fine = {"cosmors": "normal", "cosmors-fine": "fine"}
                tmp = {
                    "jobtype": "cosmors",
                    "cosmorssetup": config.external_paths["cosmorssetup"],
                    "cosmorsparam": exc_fine.get(config.smgsolv2, "normal"),
                    "cosmothermversion": config.external_paths["cosmothermversion"],
                    "ctd-param": config.cosmorsparam,
                    "copymos": str(instruction_gsolv["func"]),
                }
                if config.vapor_pressure:
                    tmp["vapor_pressure"] = True
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
                    tmp_SI = instruction_gsolv["prepinfo"]
                    instruction_gsolv["prepinfo"] = []
                name = "lowlevel additive solvation"
                folder = str(instruction_gsolv["func"]) + "/Gsolv2"
            # SMD_Gsolv
            elif config.smgsolv2 == "smd_gsolv":
                job = OrcaProc
                instruction_gsolv["jobtype"] = "smd_gsolv"
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
            if instruction_gsolv["jobtype"] in ("sp", "sp_implicit"):
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
                        conf.lowlevel_gsolv_info["range"] = {
                            conf.job["temperature"]: 0.0
                        }
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
                    if (
                        instruction_gsolv["jobtype"] in ("gbsa_gsolv", "alpb_gsolv")
                        and conf.lowlevel_sp_info["info"] == "calculated"
                        and conf.lowlevel_sp_info["method"]
                        == instruction_gsolv["method"]
                    ):
                        conf.job["energy"] = conf.lowlevel_sp_info["energy"]
                    conf.lowlevel_sp_info["energy"] = conf.job["energy"]
                    conf.lowlevel_sp_info["info"] = "calculated"
                    conf.lowlevel_sp_info["method"] = instruction_gsolv["method"]
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
                    f"{conf.lowlevel_gsolv_info['range'].get(config.temperature, 0.0):>.7f}"
                )
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))

    # create data for possible SI generation:-----------------------------------
    # Energy:
    ensembledata.si["part2"]["Energy"] = instruction_gsolv["method"]
    if "DOGCP" in instruction_gsolv["prepinfo"]:
        ensembledata.si["part2"]["Energy"] += " + GCP"
    # Energy_settings:
    try:
        if job == TmProc:
            if tmp_SI is not None:
                tmp = " ".join(qm_prepinfo["tm"][tmp_SI[0]]).replace("-", "")
            else:
                tmp = " ".join(
                    qm_prepinfo["tm"][instruction_gsolv["prepinfo"][0]]
                ).replace("-", "")
        elif job == OrcaProc:
            if tmp_SI is not None:
                tmp = " ".join(qm_prepinfo["orca"][tmp_SI[0]])
            else:
                tmp = " ".join(qm_prepinfo["orca"][instruction_gsolv["prepinfo"][0]])
        ensembledata.si["part2"]["Energy_settings"] = tmp
    except (TypeError, IndexError):
        pass
    # GmRRHO is done below!
    # Solvation:
    if config.solvent == "gas":
        ensembledata.si["part2"]["G_solv"] = "gas-phase"
    else:
        if instruction_gsolv["jobtype"] in (
            "cosmors",
            "cosmors-fine",
            "cosmors-normal",
        ):
            ensembledata.si["part2"]["G_solv"] = (
                instruction_gsolv["method2"]
                + " param= "
                + instruction_gsolv["ctd-param"]
            )
        else:
            ensembledata.si["part2"]["G_solv"] = instruction_gsolv["method2"]
    # Geometry is done above:
    # QM-CODE:
    ensembledata.si["part2"]["main QM code"] = str(config.prog).upper()
    # Threshold:
    ensembledata.si["part2"][
        "Threshold"
    ] = f"part2_threshold: {config.opt_limit} kcal/mol, Boltzmann sum threshold: {config.part2_P_threshold} %"
    # END SI generation --------------------------------------------------------

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
                conf.sym = conf.lowlevel_grrho_info.get("sym", conf.sym)
                conf.linear = conf.lowlevel_grrho_info.get("linear", conf.linear)
                conf.symnum = conf.lowlevel_grrho_info.get("symnum", conf.symnum)
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
            "scale": config.scale,
            "onlyread": config.onlyread,
        }
        instruction_rrho["method"], _ = config.get_method_name(
            "rrhoxtb",
            bhess=config.bhess,
            gfn_version=instruction_rrho["gfn_version"],
            sm=instruction_rrho["sm_rrho"],
            solvent=instruction_rrho["solvent"],
        )

        # GmRRHO for SI information:
        if config.evaluate_rrho:
            ensembledata.si["part2"]["G_mRRHO"] = instruction_rrho["method"]
            if "bhess" in ensembledata.si["part2"]["G_mRRHO"]:
                ensembledata.si["part2"]["G_mRRHO"] += " SPH"
        else:
            ensembledata.si["part2"]["G_mRRHO"] = "not included"
        # END SI

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
                            save_errors,
                        )
                    elif not os.path.isdir(tmp_to):
                        print_errors(
                            f"{'ERROR:':{WARNLEN}}Could not create folder {tmp_to}!",
                            save_errors,
                        )
                    print_errors(
                        f"{'ERROR:':{WARNLEN}}Removing conformer {conf.name}!",
                        save_errors,
                    )
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
                    conf.lowlevel_grrho_info["info"] = "failed"
                    store_confs.append(calculate.pop(calculate.index(conf)))
                else:
                    conf.sym = conf.job["symmetry"]
                    conf.linear = conf.job["linear"]
                    conf.symnum = conf.job["symnum"]
                    conf.lowlevel_grrho_info["sym"] = conf.job["symmetry"]
                    conf.lowlevel_grrho_info["linear"] = conf.job["linear"]
                    conf.lowlevel_grrho_info["symnum"] = conf.job["symnum"]
                    conf.lowlevel_grrho_info["rmsd"] = conf.job["rmsd"]
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
                    f"""{conf.get_mrrho(
                        config.temperature,
                        rrho='lowlevel_grrho_info',
                        consider_sym=True):>.8f} """
                    f"S_rot(sym)= {conf.calc_entropy_sym(config.temperature):>.7f}"
                    f""" using= {conf.get_mrrho(
                        config.temperature,
                        rrho='lowlevel_grrho_info',
                        consider_sym=config.consider_sym):>.7f}"""
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
        lambda conf: getattr(conf, "lowlevel_gsolv_info")
        .get("range", {})
        .get(config.temperature, 0.0),
        lambda conf: conf.get_mrrho(
            config.temperature, "lowlevel_grrho_info", config.consider_sym
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
    conf_in_interval(calculate)
    # --------------------------------------------------------------------------
    # ***************************************************************************
    # SD on solvation:
    # DCOSMO-RS_GSOLV
    instruction_gsolv_compare = {
        "func": config.func,
        "basis": getattr(
            config,
            "basis",
            dfa_settings.composite_method_basis.get(config.func, "def2-TZVP"),
        ),
        "charge": config.charge,
        "unpaired": config.unpaired,
        "solvent": config.solvent,
        "sm": "alpb_gsolv",
        "temperature": config.fixed_temperature,
        "energy": 0.0,
        "energy2": 0.0,
        "success": False,
        "gfn_version": config.part2_gfnv,
        "onlyread": config.onlyread,
    }
    instruction_gsolv_compare["trange"] = []
    instruction_gsolv_compare["prepinfo"] = []
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
                    gsolv_min = (
                        getattr(conf, "lowlevel_gsolv_info")
                        .get("range", {})
                        .get(config.fixed_temperature, 0.0)
                    )
                    # gsolv_min = conf.lowlevel_gsolv_info["energy"]
                    gsolv_min_id = conf.id
                    try:
                        dcosmors_gsolv_min = (
                            conf.optimization_info["energy"]
                            - conf.lowlevel_gsolv_info["gas-energy"]
                        )
                    except (TypeError, KeyError):
                        print_errors(
                            f"{'ERROR:':{WARNLEN}}Can't calculate DCOSMO-RS_gsolv. Skipping SD of Gsolv!",
                            save_errors,
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
                    print_errors(
                        f"{'ERROR:':{WARNLEN}}can't copy optimized geometry!",
                        save_errors,
                    )

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

            print(
                f"\nSD of solvation models @ {config.fixed_temperature} K (all units in kcal/mol):"
            )
            print(
                f"CONFX ΔG(COSMO-RS) ΔG(DCOSMO-RS_gsolv) ΔG(ALPB_gsolv)  "
                f"SD(COSMO-RS 40%, DCOSMO-RS_gsolv 40%, ALPB_gsolv 20%)"
            )
            print("".ljust(PLENGTH, "-"))
            pl = max([len(str(conf.id)) for conf in calculate])
            for conf in calculate:
                dgsolv = -gsolv_min + getattr(conf, "lowlevel_gsolv_info").get(
                    "range", {}
                ).get(config.fixed_temperature, 0.0)
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
            conf.calc_free_energy(
                e=e,
                solv=solv,
                rrho=rrho,
                t=temperature,
                consider_sym=config.consider_sym,
            )
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
            avGRRHO += conf.bm_weight * conf.get_mrrho(
                temperature, "lowlevel_grrho_info", config.consider_sym
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
                line = f"{temperature:^15} {avE:>14.7f}  {avG:>14.7f} "
            elif not config.evaluate_rrho:
                line = (
                    f"{temperature:^15} {avE:>14.7f} {avGsolv:>16.7f} " f"{avG:>14.7f} "
                )
            elif config.solvent == "gas":
                line = (
                    f"{temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} " f"{avG:>14.7f} "
                )
        else:
            line = (
                f"{temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                f"{avGsolv:>16.7f} {avG:>14.7f} "
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
        conf.calc_free_energy(
            e=e,
            solv=solv,
            rrho=rrho,
            t=config.temperature,
            consider_sym=config.consider_sym,
        )
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
            config.temperature,
            config.consider_sym,
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
        if sumup >= (config.part2_P_threshold / 100):
            if conf.bm_weight < (1 - (config.part2_P_threshold / 100)):
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
            f"of {config.part2_P_threshold}%:"
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
        config.temperature,
        config.consider_sym,
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
                        conf.free_energy,
                        conf.get_mrrho(
                            config.temperature,
                            rrho="lowlevel_grrho_info",
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

    tmp = int((PLENGTH - len("END of Part2")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part2" + "".rjust(tmp, "<"))
    if config.progress:
        print("#>>># CENSO: Finished part2", file=sys.stderr)
    return config, calculate, store_confs, ensembledata
