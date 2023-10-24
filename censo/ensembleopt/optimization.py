"""
Optimization == part2
performing ensembleopt of the CRE and provide low level free energies.
"""
import os
from typing import List
from functools import reduce

from censo.core import CensoCore
from censo.parallel import ProcessHandler
from censo.part import CensoPart
from censo.datastructure import GeometryData, MoleculeData
from censo.utilities import (
    print,
    timeit,
    DfaHelper,
    format_data,
)
from censo.params import (
    SOLV_MODS,
    GSOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
    AU2KCAL,
)


class Optimization(CensoPart):
    alt_name = "part2"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    __gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    _options = {
        "optimization": {
            "optcycles": {
                "default": 8,
                "range": [
                    1,
                    10
                ]
            },
            # "radsize": { # ???
            #     "default": 10,
            #     "range": [
            #         1,
            #         100
            #     ]
            # },
            "maxcyc": {
                "default": 200,
                "range": [
                    10,
                    1000
                ]
            },
            "threshold": {
                "default": 1.5,
                "range": [
                    0.5,
                    5
                ]
            },
            "hlow": {
                "default": 0.01,
                "range": [
                    0.001,
                    0.1
                ]
            },
            "gradthr": {
                "default": 0.01,
                "range": [
                    0.01,
                    1.0
                ]
            },
            "boltzmannthr": {  # boltzmann sum threshold
                "default": 85.0,
                "range": [
                    1.0,
                    99.9
                ]
            },
            # "spearmanthr": {
            #     "default": 0.0,
            #     "range": [
            #         -1.0,
            #         1.0
            #     ]
            # },
            "func": {
                "default": "r2scan-3c",
                "options": DfaHelper.find_func("optimization")
            },
            "basis": {
                "default": "def2-TZVP",
                "options": BASIS_SETS
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "sm": {
                "default": "smd",
                "options": __solv_mods
            },
            "smgsolv": {
                "default": "smd",
                "options": __gsolv_mods
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "optlevel": {
                "default": "normal",
                "options": [
                    "crude",
                    "sloppy",
                    "loose",
                    "lax",
                    "normal",
                    "tight",
                    "vtight",
                    "extreme",
                ]
            },
            "grid": {
                "default": "high",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": True
            },
            "gcp": {
                "default": True
            },
            "opt_spearman": {
                "default": True
            },
            "crestcheck": {
                "default": False
            }
        },
    }

    _settings = {}

    def __init__(self, core: CensoCore, name: str = "optimization"):
        super().__init__(core, name=name)
        self.confs_nc: List[GeometryData]

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        Optimization of the ensemble at DFT level (possibly with implicit solvation)

        Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
        by calculating the thermodynamics of the optimized geometries

        Alternatively just run the complete geometry optimization for every conformer with xtb as driver (decide with 'opt_spearman')

        TODO - implement regular ensembleopt (no xtb driver)
        TODO - what happens if not a single conformer converges?
        """

        """
        don't use this for now, use new keyword 'maxcyc' instead
        if config.nat > 200:
            # stopcycle = don't optimize more than stopcycle cycles
            stopcycle = config.nat * 2
        else:
            stopcycle = 200
        """
        # NOTE: (IMPORTANT) the following only uses xtb as driver (no native geometry optimizations)

        # print instructions
        self.print_info()

        # setup process handler
        handler = ProcessHandler(self.__instructions)

        # decide for doing spearman ensembleopt or standard ensembleopt (TODO)
        if self.__instructions["opt_spearman"] and len(self.core.conformers) > 1:
            """
            ensembleopt using macrocycles with 'optcycles' microcycles
            """
            # if not self.args.spearmanthr:
            #     # set spearmanthr by number of atoms:
            #     self.spearmanthr = 1 / (exp(0.03 * (self.runinfo["nat"] ** (1 / 4))))

            self.__spearman_opt(handler)
        else:
            """
            do normal geometry optimization
            """
            if not len(self.core.conformers) > 1:
                print(f"Only one conformer ({self.core.conformers[0].name}) is available for ensembleopt.")

            # disable spearman ensembleopt
            print("Spearman ensembleopt turned off.")
            self.__instructions["opt_spearman"] = False

            # run optimizations using xtb as driver
            handler.conformers = [conf.geom for conf in self.core.conformers]
            results_opt = handler.execute(["xtb_opt"], self.dir)

            # update results for each conformer
            for conf in self.core.conformers:
                # store mo_path if 'copy_mo' is enabled
                if self.__instructions.get("copy_mo", None):
                    conf.geom.mo_path = results_opt[id(conf)]["xtb_opt"]["mo_path"]

                # update geometry of the conformer
                conf.geom.xyz = results_opt[id(conf)]["xtb_opt"]["geom"]

                # store results
                conf.results.setdefault(self.__name.lower(), {}).update(results_opt[id(conf)])

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?).
        # we don't do that anymore and just use the final energies from the optimizations, which are done using a solvent model,
        # since this basically makes no difference in comp time
        # do rrho on converged geometries (overwrites previous rrho calculations)
        handler.conformers = [conf.geom for conf in self.core.conformers]
        results_rrho = handler.execute(["xtb_rrho"], self.dir)

        for conf in self.core.conformers:
            conf.results[self.__name.lower()].update(results_rrho[id(conf)])

        # calculate boltzmann weights from gtot values calculated here
        self.core.calc_boltzmannweights(
            self.__instructions.get("temperature", 298.15),
            self.__name.lower()
        )

        # write final results
        self.write_results()

        # dump ensemble
        self.core.dump_ensemble(self.__name.lower())

    def grrho(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy (last step of running ensembleopt) and Gmrrho
        If no GmRRHO is calculated only the most recent DFT energy is returned
        """
        # TODO - implement branch for standard ensembleopt (no ancopt)
        try:
            return conf.results[self.__name.lower()]["xtb_opt"]["energy"] + \
                conf.results[self.__name.lower()]["xtb_rrho"]["energy"]
        except KeyError:
            return conf.results[self.__name.lower()]["xtb_opt"]["energy"]

    def __spearman_opt(self, handler: ProcessHandler):
        # make a separate list of conformers that only includes (considered) conformers that are not converged
        # NOTE: this is a special step only necessary for spearman ensembleopt
        # at this point it's just self.core.conformers
        self.confs_nc = [conf.geom for conf in self.core.conformers]

        stopcond_converged = False
        ncyc = 0
        rrho_done = False
        print(f"Optimization using Spearman threshold, {self.__instructions['optcycles']} cycles per step.")
        print(f"NCYC: {ncyc}")
        while (
                not stopcond_converged
                and ncyc < self.__instructions["maxcyc"]
                # TODO - maybe make this more intelligent:
                # make maxcyc lower and if some apparently relevant conformer doesn't converge within it's chunk,
                # move it to a new chunk and calculate later
        ):
            # NOTE: this loop works through confs_nc, so if the ensembleopt for a conformer is converged, all the following steps will not consider it anymore
            # update conformers for ProcessHandler
            handler.conformers = self.confs_nc

            # run optimizations for 'optcycles' steps
            results_opt = handler.execute(["xtb_opt"], self.dir)

            # put ensembleopt results into conformer objects
            for conf in self.confs_nc:
                for coreconf in self.core.conformers:
                    if id(coreconf) == conf.id:
                        # store mo_path if 'copy_mo' is enabled
                        if self.__instructions.get("copy_mo", None):
                            coreconf.geom.mo_path = results_opt[conf.id]["xtb_opt"]["mo_path"]

                        # update geometry of the conformer
                        coreconf.geom.xyz = results_opt[conf.id]["xtb_opt"]["geom"]
                        conf.xyz = results_opt[conf.id]["xtb_opt"]["geom"]

                        # store results
                        coreconf.results.setdefault(self.__name.lower(), {}).update(results_opt[conf.id])
                        break

            # run xtb_rrho for finite temperature contributions
            # for now only after the first 'optcycles' steps or after at least 6 cycles are done
            # TODO - make this better
            if ncyc + self.__instructions["optcycles"] >= 6 and not rrho_done:

                # evaluate rrho using bhess flag (reset after this)
                tmp = self.__instructions["bhess"]
                self.__instructions["bhess"] = True
                results_rrho = handler.execute(["xtb_rrho"], self.dir)
                self.__instructions["bhess"] = tmp

                # put results into conformer objects
                for conf in self.confs_nc:
                    for coreconf in self.core.conformers:
                        if id(coreconf) == conf.id:
                            coreconf.results[self.__name.lower()].update(results_rrho[conf.id])
                            coreconf.results[self.__name.lower()]["gtot"] = self.grrho(coreconf)
                            break

                # flag to make sure that rrho is only calculated once
                rrho_done = True

            # TODO - crestcheck each iteration if ncyc >= 6

            # kick out conformers above threshold
            threshold = self.__instructions["threshold"] / AU2KCAL

            """
            # TODO TODO TODO TODO - update threshold based on spearman coefficients (???) - leave this out for now
            minthreshold = self.__instructions["threshold"] / AU2KCAL
            # 'decyc' contains the difference between the minimal gtot and the gtot of the conformer for every ensembleopt cycle
            # pseudo: decyc[i] = gtot[i] - mingtot[i]
            # 'toevaluate' contains the indices of the cycles to be evaluated for spearman correlation computation (always the most current cycles)
            # this only goes unitl the third-last cycle evaluated (since deprevious and decurrent are three steps apart (why??)), can also be just one cycle
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
                        f"{i + num_eval:>{digits1}}"
                        f" =  {spearman_v:>.4f}"
                    )
                    evalspearman.append(spearman_v)
                else:
                    print(
                        f"{'':>10} Spearman coeff. from {i:>{digits1}} --> "
                        f"{i + num_eval:>{digits1}}"
                        f" =  {spearman_v:>.4f}"
                    )
            print(
                f"Final averaged Spearman correlation coefficient: "
                f"{(sum(evalspearman) / 2):>.4f}"
            )
            cycle_spearman.append(f"{sum(evalspearman) / 2:.3f}")
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
                        f"Updated ensembleopt threshold to: {ewin:.2f} kcal/mol"
                    )
                else:
                    print(
                        f"Current ensembleopt threshold: {ewin:.2f} kcal/mol"
                    )
            """

            # sort conformers and pick the free enthalpy of the lowest conformer as limit
            self.core.conformers.sort(key=lambda conf: self.grrho(conf))
            limit = min([self.grrho(conf) for conf in self.core.conformers])

            # filter out all conformers above threshold and with a gradient norm smaller than 'gradthr'
            # so that 'filtered' contains all conformers that should not be considered any further
            f = lambda x: self.grrho(x) - limit > threshold and x.results[self.__name.lower()]["xtb_opt"][
                "grad_norm"] < self.__instructions["gradthr"]
            filtered: List[MoleculeData] = [
                conf for conf in filter(
                    f,
                    self.core.conformers
                )
            ]

            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(filtered)

            # also remove conformers from confs_nc
            for conf in filtered:
                print(f"{conf.name} is no longer considered (δG = {(self.grrho(conf) - limit) * AU2KCAL:.2f}).")
                self.confs_nc.remove(conf.geom)

            # update list of converged conformers
            for conf in self.confs_nc:
                if results_opt[conf.id]["xtb_opt"]["converged"]:
                    print(f"{conf.name} converged after {ncyc + results_opt[conf.id]['xtb_opt']['cycles']} steps.")
                    self.confs_nc.remove(conf)

            # check if all (considered) conformers converged - End of While-loop
            stopcond_converged = len(self.confs_nc) == 0

            # TODO - print out information about current state of the ensemble
            self.print_update()

            # update number of cycles
            ncyc += self.__instructions["optcycles"]
            print(f"NCYC: {ncyc}")

    def write_results(self) -> None:
        """
        formatted write of part results (optional)
        """
        # column headers
        headers = [
            "CONF#",
            "E (DFT) (+ ΔGsolv)",
            "ΔE (DFT) (+ δΔGsolv)",
            "GmRRHO",
            "Gtot",
            "ΔGtot",
            "Boltzmann weight",
        ]

        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
            f"\% at {self.__instructions.get('temperature', 298.15)} K",
        ]

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            self.grrho(conf)
            for conf in self.core.conformers
        )

        dftmin = min(
            conf.results[self.__name.lower()]["xtb_opt"]["energy"]
            for conf in self.core.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda
                conf: f"{conf.results[self.__name.lower()]['xtb_opt']['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda
                conf: f"{(conf.results[self.__name.lower()]['xtb_opt']['energy'] - dftmin) * AU2KCAL:.2f}",
            "GmRRHO": lambda conf:
            f"{conf.results[self.__name.lower()]['xtb_rrho']['gibbs'][self.__instructions['temperature']]:.6f}"
            if self.__instructions["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self.__name.lower()]['bmw'] * 100:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # write lines to file
        with open(os.path.join(self.core.workdir, f"{self.__name.lower()}.out"), "w",
                  newline=None) as outfile:
            outfile.writelines(lines)

    def print_update(self) -> None:
        """
        writes information about the current state of the ensemble
        """
        # TODO
        for conf in self.confs_nc:
            for coreconf in self.core.conformers:
                if id(coreconf) == conf.id:
                    print(f"{coreconf.name}: {self.grrho(coreconf)}")
