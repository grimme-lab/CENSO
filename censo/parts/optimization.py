"""
Optimization == part2
performing optimization of the CRE and provide low level free energies.
"""
import os
import sys
from copy import deepcopy
from censo.cfg import PLENGTH, CODING, AU2KCAL, DIGILEN, WARNLEN
from censo.utilities import (
    last_folders,
    print,
    timeit,
    spearman,
)
from typing import List

from censo.core import CensoCore
from censo.settings import CensoSettings
from censo.parallel import ProcessHandler
from censo.parts.part import CensoPart
from censo.datastructure import GeometryData, MoleculeData


class Optimization(CensoPart):
    
    alt_name = "part2"

    def __init__(self, core: CensoCore, settings: CensoSettings):
        super().__init__(core, settings, "optimization")
        self.confs_nc: List[GeometryData]


    @timeit
    def run(self) -> None:
        """
        Optimization of the ensemble at DFT level (possibly with implicit solvation)

        Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
        by calculating the thermodynamics of the optimized geometries

        Alternatively just run the complete optimization for every conformer with xtb as driver (decide with 'opt_spearman')

        TODO - implement regular optimization (no xtb driver)
        TODO - what happens if not a single confomer converges?
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
        handler = ProcessHandler(self._instructions)

        # decide for doing spearman optimization or standard optimization (TODO)
        if self._instructions["opt_spearman"] and len(self.core.conformers) > 1:
            """
            optimization using spearman threshold, updated every 'optcycles' steps
            """
            self.__spearman_opt(handler)
        else:
            """
            do normal optimization
            """
            if not len(self.core.conformers) > 1:
                print(f"Only one conformer ({self.core.conformers[0].name}) is available for optimization.")

            # disable spearman optimization
            print("Spearman optimization turned off.")
            self._instructions["opt_spearman"] = False

            # run optimizations using xtb as driver
            handler.conformers = [conf.geom for conf in self.core.conformers]
            results_opt = handler.execute(["xtb_opt"], folder)

            # update results for each conformer
            for conf in self.core.conformers:
                # store mo_path if 'copy_mo' is enabled
                if self._instructions.get("copy_mo", None):
                    conf.geom.mo_path = results_opt[id(conf)]["xtb_opt"]["mo_path"]

                # update geometry of the conformer
                conf.geom.xyz = results_opt[id(conf)]["xtb_opt"]["geom"]

                # store results
                conf.results.setdefault(self.__class__.__name__.lower(), {}).update(results_opt[id(conf)])

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?).
        # we don't do that anymore and just use the final energies from the optimizations, which are done using a solvent model,
        # since this basically makes no difference in comp time
        # do rrho on converged geometries (overwrites previous rrho calculations)
        handler.conformers = [conf.geom for conf in self.core.conformers]
        results_rrho = handler.execute(["xtb_rrho"], folder)

        for conf in self.core.conformers:
            conf.results[self.__class__.__name__.lower()].update(results_rrho[id(conf)])

        # calculate boltzmann weights from gtot values calculated here
        self.core.calc_boltzmannweights(
            self._instructions.get("temperature", 298.15),
            self.__class__.__name__.lower()
        )

        # TODO - write final results
        self.write_results()


    def grrho(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy (last step of running optimization) and Gmrrho
        """
        # TODO - implement branch for standard optimization (no ancopt)
        # TODO - implement this for usage without rrho
        return conf.results[self.__class__.__name__.lower()]["xtb_opt"]["energy"] + conf.results[self.__class__.__name__.lower()]["xtb_rrho"]["energy"]


    def __spearman_opt(self, handler: ProcessHandler):
        # make a separate list of conformers that only includes (considered) conformers that are not converged
        # NOTE: this is a special step only necessary for spearman optimization
        # at this point it's just self.core.conformers
        self.confs_nc = [conf.geom for conf in self.core.conformers]

        stopcond_converged = False
        ncyc = 0
        rrho_done = False
        print(f"Optimization using Spearman threshold, {self._instructions['optcycles']} cycles per step.")
        print(f"NCYC: {ncyc}")
        while (
                not stopcond_converged
                and ncyc < self._instructions["maxcyc"]
                # TODO - maybe make this more intelligent:
                # make maxcyc lower and if some apparently relevant conformer doesn't converge within it's chunk,
                # move it to a new chunk and calculate later
        ):
            # NOTE: this loop works through confs_nc, so if the optimization for a conformer is converged, all the following steps will not consider it anymore
            # update conformers for ProcessHandler
            handler.conformers = self.confs_nc

            # update number of cycles
            ncyc += self._instructions["optcycles"]

            # run optimizations for 'optcycles' steps
            results_opt = handler.execute(["xtb_opt"], folder)

            # put optimization results into conformer objects
            for conf in self.confs_nc:
                for coreconf in self.core.conformers:
                    if id(coreconf) == conf.id:
                        # store mo_path if 'copy_mo' is enabled
                        if self._instructions.get("copy_mo", None):
                            coreconf.geom.mo_path = results_opt[conf.id]["xtb_opt"]["mo_path"]

                        # update geometry of the conformer
                        coreconf.geom.xyz = results_opt[conf.id]["xtb_opt"]["geom"]
                        conf.xyz = results_opt[conf.id]["xtb_opt"]["geom"]

                        # store results
                        coreconf.results.setdefault(self.__class__.__name__.lower(), {}).update(results_opt[conf.id])
                        break

            # run xtb_rrho for finite temperature contributions
            # for now only after the first 'optcycles' steps or after at least 6 cycles are done
            # TODO - make this better
            if ncyc >= 6 and not rrho_done:
                # NOTE: conforms to general settings (if 'bhess' is false this will run 'ohess' at this point)
                results_rrho = handler.execute(["xtb_rrho"], folder)

                # put results into conformer objects
                for conf in self.confs_nc:
                    for coreconf in self.core.conformers:
                        if id(coreconf) == conf.id:
                            coreconf.results[self.__class__.__name__.lower()].update(results_rrho[conf.id])
                            coreconf.results[self.__class__.__name__.lower()]["gtot"] = self.grrho(coreconf)
                            break

                rrho_done = True

            # TODO - crestcheck each iteration if ncyc >= 6

            # sort conformers
            self.core.conformers.sort(key=lambda conf: conf.results[self.__class__.__name__.lower()]["gtot"])

            # kick out conformers above threshold
            threshold = self._instructions["threshold"] / AU2KCAL

            """
            # TODO TODO TODO TODO - update threshold based on spearman coefficients (???) - leave this out for now
            minthreshold = self._instructions["threshold"] / AU2KCAL
            # 'decyc' contains the difference between the minimal gtot and the gtot of the conformer for every optimization cycle
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
                        f"Updated optimization threshold to: {ewin:.2f} kcal/mol"
                    )
                else:
                    print(
                        f"Current optimization threshold: {ewin:.2f} kcal/mol"
                    )
            """

            # pick the free enthalpy of the lowest conformer
            limit = min([conf.results[self.__class__.__name__.lower()]["gtot"] for conf in self.core.conformers])

            # filter out all conformers above threshold and with a gradient norm smaller than 'gradthr'
            # so that 'filtered' contains all conformers that should not be considered any further
            f = lambda x: self.grrho(x) - limit > threshold and x.results[self.__class__.__name__.lower()]["xtb_opt"][
                "grad_norm"] < self._instructions["gradthr"]
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
            f"\% at {self._instructions.get('temperature', 298.15)} K",
        ]

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            self.grrho(conf)
            for conf in self.core.conformers
        )

        dftmin = min(
            conf.results[self.__class__.__name__.lower()]["xtb_opt"]["energy"]
            for conf in self.core.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['xtb_opt']['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['xtb_opt']['energy'] - dftmin) * AU2KCAL:.2f}",
            "GmRRHO": lambda conf:
            f"{conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}"
            if self._instructions["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['bmw'] * 100:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # append lines to already existing file
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.out"), "a", newline=None) as outfile:
            outfile.writelines(lines)

    def print_update(self) -> None:
        """
        writes information about the current state of the ensemble
        """
        # TODO
        for conf in self.confs_nc:
            for coreconf in self.core.conformers:
                if id(coreconf) == conf.id:
                    print(f"{coreconf.name}: {coreconf.results[self.__class__.__name__.lower()]['xtb_opt']['energy']}")
