"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""
import os
from functools import reduce
from math import exp
from statistics import stdev, mean
from typing import List

from src.core import CensoCore
from src.datastructure import MoleculeData
from src.parallel import execute
from src.params import (
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
    AU2KCAL,
)
from src.part import CensoPart
from src.utilities import (
    print,
    timeit,
    DfaHelper,
    format_data, setup_logger,
)

logger = setup_logger(__name__)


class Optimization(CensoPart):
    alt_name = "part2"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "optcycles": {
            "default": 8,
            "range": [
                1,
                10
            ]
        },
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
        "macrocycles": {
            "default": True
        },
        "crestcheck": {
            "default": False
        },
        "template": {
            "default": False
        },
    }

    _settings = {}

    def __init__(self, core: CensoCore):
        super().__init__(core)
        self.confs_nc: List[MoleculeData]

    @timeit
    @CensoPart._create_dir
    def run(self) -> None:
        """
        Optimization of the ensemble at DFT level (possibly with implicit solvation)

        Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
        by calculating the thermodynamics of the optimized geometries

        Alternatively just run the complete geometry optimization for every conformer with xtb as driver (decide with 'macrocycles')

        TODO - implement regular geometry optimization (no xtb driver)
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

        self._instructions["jobtype"] = ["xtb_opt"]

        # decide for doing spearman ensembleopt or standard geometry optimization (TODO)
        if self._instructions["macrocycles"] and len(self.core.conformers) > 1:
            """
            ensembleopt using macrocycles with 'optcycles' microcycles
            """
            # if not self.args.spearmanthr:
            #     # set spearmanthr by number of atoms:
            #     self.spearmanthr = 1 / (exp(0.03 * (self.runinfo["nat"] ** (1 / 4))))

            self.__macrocycle_opt()
        else:
            """
            do normal geometry optimization
            """
            if not len(self.core.conformers) > 1:
                print(f"Only one conformer ({self.core.conformers[0].name}) is available for optimization.")

            # disable spearman optimization
            print("Macrocycle optimization turned off.")
            self._instructions["macrocycles"] = False

            # run optimizations using xtb as driver
            results_opt = execute(self.core.conformers, self._instructions, self.dir)

            # update results for each conformer
            for conf in self.core.conformers:
                # update geometry of the conformer
                conf.geom.xyz = results_opt[id(conf)]["xtb_opt"]["geom"]

                # store results
                conf.results.setdefault(self._name, {}).update(results_opt[id(conf)])

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?).
        # we don't do that anymore and just use the final energies from the optimizations, which are done using a
        # solvent model, since this basically makes no difference in comp time
        # do rrho on converged geometries (overwrites previous rrho calculations)
        results_rrho = execute(self.core.conformers, self._instructions, self.dir)

        for conf in self.core.conformers:
            conf.results[self._name].update(results_rrho[id(conf)])

        # calculate boltzmann weights from gtot values calculated here
        self.core.calc_boltzmannweights(
            self._instructions.get("temperature", 298.15),
            self._name
        )

        # write final results
        self.write_results()

        # dump ensemble
        self.core.dump_ensemble(self._name)

    def grrho(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy (last step of running optimization) and Gmrrho
        If no GmRRHO is calculated only the most recent DFT energy is returned
        """
        # TODO - implement branch for standard geometry optimization (no ancopt)
        try:
            return conf.results[self._name]["xtb_opt"]["energy"] + \
                conf.results[self._name]["xtb_rrho"]["energy"]
        except KeyError:
            return conf.results[self._name]["xtb_opt"]["energy"]

    def __macrocycle_opt(self):
        # make a separate list of conformers that only includes (considered) conformers that are not converged
        # NOTE: this is a special step only necessary for macrocycle optimization
        # at this point it's just self.core.conformers
        self.confs_nc = self.core.conformers.copy()

        ncyc = 0
        rrho_done = False
        print(f"Optimization using macrocycles, {self._instructions['optcycles']} microcycles per step.")
        print(f"NCYC: {ncyc}")
        while (
                len(self.confs_nc) > 0
                and ncyc < self._instructions["maxcyc"]
        ):
            # NOTE: this loop works through confs_nc, so if the geometry optimization for a conformer is converged,
            # all the following steps will not consider it anymore
            # run optimizations for 'optcycles' steps
            results_opt = execute(self.confs_nc, self._instructions, self.dir)

            # put geometry optimization results into conformer objects
            for conf in self.confs_nc:
                # update geometry of the conformer
                conf.geom.xyz = results_opt[conf.geom.id]["xtb_opt"]["geom"]

                # store results
                conf.results.setdefault(self._name, {}).update(results_opt[conf.geom.id])

            # run xtb_rrho for finite temperature contributions
            # for now only after the first 'optcycles' steps or after at least 6 cycles are done
            # TODO - make this better
            if ncyc + self._instructions["optcycles"] >= 6 and not rrho_done:

                # evaluate rrho using bhess flag (reset after this)
                tmp = self._instructions["bhess"]
                self._instructions["bhess"] = True
                self._instructions["jobtype"] = ["xtb_rrho"]
                results_rrho = execute(self.confs_nc, self._instructions, self.dir)
                self._instructions["jobtype"] = ["xtb_opt"]
                self._instructions["bhess"] = tmp

                # put results into conformer objects
                for conf in self.confs_nc:
                    conf.results[self._name].update(results_rrho[conf.geom.id])
                    conf.results[self._name]["gtot"] = self.grrho(conf)

                # flag to make sure that rrho is only calculated once
                rrho_done = True

            # TODO - crestcheck each iteration if ncyc >= 6

            # kick out conformers above threshold
            threshold = self._instructions["threshold"] / AU2KCAL

            # make threshold fuzzy based on stdev of the gradient norms (adds at most 2 kcal/mol to the threshold)
            if len(self.confs_nc) > 1:
                threshold = threshold + (2 / AU2KCAL) * (1 - exp(
                                 - AU2KCAL * stdev(
                                     [conf.results[self._name]["xtb_opt"]["grad_norm"] for conf in self.confs_nc]
                                 )
                            ))

                print(f"Current fuzzy threshold: {threshold * AU2KCAL:.2f} kcal/mol "
                      f"(min: {self._instructions['threshold']:.2f} kcal/mol).")
            else:
                print(f"Using minimal threshold ({self._instructions['threshold']:.2f} kcal/mol).")

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
                        f"Updated ensembleopt threshold to: {ewin:.2f} kcal/mol"
                    )
                else:
                    print(
                        f"Current ensembleopt threshold: {ewin:.2f} kcal/mol"
                    )
            """

            # sort conformers
            self.core.conformers.sort(key=lambda conf: self.grrho(conf))

            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(
                self.grrho, threshold,
                additional_filter=lambda x:
                x.results[self._name]["xtb_opt"]["grad_norm"] < self._instructions["gradthr"]
            )

            # also remove conformers from confs_nc
            limit = min(self.grrho(conf) for conf in self.core.conformers)
            for conf in self.core.rem:
                if conf in self.confs_nc:
                    print(f"{conf.name} is no longer considered (small gradient and"
                          f" δG = {(self.grrho(conf) - limit) * AU2KCAL:.2f}).")
                    self.confs_nc.remove(conf)

            # update list of converged conformers
            for conf in self.confs_nc:
                if results_opt[conf.geom.id]["xtb_opt"]["converged"]:
                    print(f"{conf.name} converged after {ncyc + results_opt[conf.geom.id]['xtb_opt']['cycles']} steps.")
                    self.confs_nc.remove(conf)

            # TODO - print out information about current state of the ensemble
            self.print_update()

            # update number of cycles
            ncyc += self._instructions["optcycles"]
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
            conf.results[self._name]["xtb_opt"]["energy"]
            for conf in self.core.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda conf: f"{conf.results[self._name]['xtb_opt']['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda conf: f"{(conf.results[self._name]['xtb_opt']['energy'] - dftmin) * AU2KCAL:.2f}",
            "GmRRHO": lambda conf:
            f"{conf.results[self._name]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}"
            if self._instructions["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # write lines to file
        global logger
        logger.debug(f"Writing to {os.path.join(self.core.workdir, f'{self._name}.out')}.")
        with open(os.path.join(self.core.workdir, f"{self._name}.out"), "w", newline=None) as outfile:
            outfile.writelines(lines)

    def print_update(self) -> None:
        """
        writes information about the current state of the ensemble
        """
        # TODO
        for conf in self.confs_nc:
            print(f"{conf.name}: {self.grrho(conf)}")
