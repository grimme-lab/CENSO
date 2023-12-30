"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""
import os
from functools import reduce
from math import exp

from ..core import EnsembleData
from ..datastructure import MoleculeData
from ..parallel import execute
from ..params import (
    SOLV_MODS,
    PROGS,
    BASIS_SETS,
    GRIDOPTIONS,
    GFNOPTIONS,
    AU2KCAL,
    SOLVENTS_DB,
)
from ..part import CensoPart
from ..utilities import (
    print,
    timeit,
    DfaHelper,
    format_data,
    setup_logger,
    mean_similarity,
)

logger = setup_logger(__name__)


class Optimization(CensoPart):
    alt_name = "part2"

    __solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())

    _options = {
        "optcycles": {"default": 8, "range": [1, 10]},
        "maxcyc": {"default": 200, "range": [10, 1000]},
        "threshold": {"default": 1.5, "range": [0.5, 5]},
        "hlow": {"default": 0.01, "range": [0.001, 0.1]},
        "gradthr": {"default": 0.01, "range": [0.01, 1.0]},
        "func": {
            "default": "r2scan-3c",
            "options": DfaHelper.find_func("optimization"),
        },
        "basis": {"default": "def2-TZVP", "options": BASIS_SETS},
        "prog": {"default": "orca", "options": PROGS},
        "sm": {"default": "smd", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
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
            ],
        },
        "grid": {"default": "high", "options": GRIDOPTIONS},
        "run": {"default": True},
        "gcp": {"default": True},
        "macrocycles": {"default": True},
        "crestcheck": {"default": False},
        "template": {"default": False},
    }

    _settings = {}

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

        # Special 'todo-list' for optimization part, contains all unconverged conformers,
        # used in macrocycle optimization
        self.confs_nc: list[MoleculeData]

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

        # decide for doing spearman ensembleopt or standard geometry optimization (TODO)
        if self.get_settings()["macrocycles"] and len(self.ensemble.conformers) > 1:
            """
            ensembleopt using macrocycles with 'optcycles' microcycles
            """
            self.__macrocycle_opt()
        else:
            """
            do normal geometry optimization
            """

            jobtype = ["xtb_opt"]
            prepinfo = self.setup_prepinfo(jobtype)

            if not len(self.ensemble.conformers) > 1:
                print(
                    f"Only one conformer ({self.ensemble.conformers[0].name}) is available for optimization."
                )

            # disable spearman optimization
            print("Macrocycle optimization turned off.")
            self.set_setting("macrocycles", False)

            # run optimizations using xtb as driver
            results_opt, failed = execute(
                self.ensemble.conformers,
                prepinfo,
                self.dir, self.get_settings()["prog"],
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                omp=self.get_general_settings()["omp"],
                maxcores=self.get_general_settings()["maxcores"],
                retry_failed=self.get_general_settings()["retry_failed"],
            )

            # Remove failed conformers
            for confid in failed:
                self.ensemble.remove_conformers(failed)

            # update results for each conformer
            for conf in self.ensemble.conformers:
                # update geometry of the conformer
                conf.geom.xyz = results_opt[id(conf)]["xtb_opt"]["geom"]

                # store results
                conf.results.setdefault(self._name, {}).update(
                    results_opt[id(conf)])

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?).
        # we don't do that anymore and just use the final energies from the optimizations, which are done using a
        # solvent model, since this basically makes no difference in comp time
        # do rrho on converged geometries (overwrites previous rrho calculations)
        jobtype = ["xtb_rrho"]
        prepinfo = self.setup_prepinfo(jobtype)
        results_rrho = execute(
            self.ensemble.conformers,
            prepinfo,
            self.dir, self.get_settings()["prog"],
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            omp=self.get_general_settings()["omp"],
            maxcores=self.get_general_settings()["maxcores"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        for confid in failed:
            # Remove failed conformers
            self.ensemble.remove_conformers(failed)

        for conf in self.ensemble.conformers:
            conf.results[self._name].update(results_rrho[id(conf)])

        # calculate boltzmann weights from gtot values calculated here
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        # write final results
        self.write_results()

        # dump ensemble
        self.ensemble.dump_ensemble(self._name)

    def grrho(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy (last step of running optimization) and Gmrrho
        If no GmRRHO is calculated only the most recent DFT energy is returned

        Args:
            conf: conformer object

        Returns:
            float: Gtot (in the unit of the results)
        """
        # TODO - implement branch for standard geometry optimization (no ancopt)
        try:
            return (
                conf.results[self._name]["xtb_opt"]["energy"]
                + conf.results[self._name]["xtb_rrho"]["energy"]
            )
        except KeyError:
            return conf.results[self._name]["xtb_opt"]["energy"]

    def setup_prepinfo(self, jobtype: list[str]) -> dict[str, any]:
        prepinfo = {jt: {} for jt in jobtype}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        if "xtb_opt" in jobtype:
            prepinfo["xtb_opt"] = {
                "func_name": DfaHelper.get_name(
                    self.get_settings()["func"], self.get_settings()["prog"]
                ),
                "func_type": DfaHelper.get_type(
                    self.get_settings()["func"]),
                "disp": DfaHelper.get_disp(
                    self.get_settings()["func"]),
                "basis": self.get_settings()["basis"],
                "grid": self.get_settings()["grid"],
                "template": self.get_settings()["template"],
                "gcp": self.get_settings()["gcp"],
                "sm": self.get_settings()["sm"],
                "solvent_key_prog": SOLVENTS_DB.get(self.get_general_settings()["solvent"])[self.get_settings()["sm"]][1],
                "optcycles": self.get_settings()["optcycles"],
                "hlow": self.get_settings()["hlow"],
                "optlevel": self.get_settings()["optlevel"],
                "macrocycles": self.get_settings()["macrocycles"],
            }

        if "xtb_rrho" in jobtype:
            prepinfo["xtb_rrho"] = {
                "gfnv": self.get_settings()["gfnv"],
                "solvent_key_xtb": SOLVENTS_DB.get(self.get_general_settings()["solvent"])["xtb"][1],
            }

        return prepinfo

    def __macrocycle_opt(self):
        """
        Macrocycle optimization using xtb as driver. Also calculates GmRRHO for finite temperature contributions
        and uses adaptive threshold based on mean trajectory similarity.
        """
        # make a separate list of conformers that only includes (considered) conformers that are not converged
        # NOTE: this is a special step only necessary for macrocycle optimization
        # at this point it's just self.ensemble.conformers, it is basically a todo-list
        self.confs_nc = self.ensemble.conformers.copy()

        jobtype = ["xtb_opt"]
        prepinfo = self.setup_prepinfo(jobtype)

        ncyc = 0
        rrho_done = False
        print(
            f"Optimization using macrocycles, {self.get_settings()['optcycles']} microcycles per step."
        )
        print(f"NCYC: {ncyc}")
        while len(self.confs_nc) > 0 and ncyc < self.get_settings()["maxcyc"]:
            # NOTE: this loop works through confs_nc, so if the geometry optimization for a conformer is converged,
            # and therefore removed from our 'todo-list', all the following steps will not consider it anymore
            # run optimizations for 'optcycles' steps
            results_opt, failed = execute(
                self.ensemble.conformers,
                prepinfo,
                self.dir, self.get_settings()["prog"],
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                omp=self.get_general_settings()["omp"],
                maxcores=self.get_general_settings()["maxcores"],
                retry_failed=self.get_general_settings()["retry_failed"],
            )

            # put geometry optimization results into conformer objects
            for conf in self.confs_nc:
                # update geometry of the conformer
                conf.geom.xyz = results_opt[conf.geom.id]["xtb_opt"]["geom"]

                # store results
                conf.results.setdefault(self._name, {}).update(
                    results_opt[conf.geom.id]
                )

            # run xtb_rrho for finite temperature contributions
            # for now only after the first 'optcycles' steps or after at least 6 cycles are done
            # TODO - make this better
            if ncyc + self.get_settings()["optcycles"] >= 6 and not rrho_done:
                # evaluate rrho using bhess flag (reset after this)
                # TODO - this smells a bit
                tmp = self.get_general_settings()["bhess"]
                self.set_general_setting("bhess", True)
                jobtype = ["xtb_rrho"]
                results_rrho, failed = execute(
                    self.ensemble.conformers,
                    self.setup_prepinfo(jobtype),
                    self.dir, self.get_settings()["prog"],
                    copy_mo=self.get_general_settings()["copy_mo"],
                    balance=self.get_general_settings()["balance"],
                    omp=self.get_general_settings()["omp"],
                    maxcores=self.get_general_settings()["maxcores"],
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Reset
                jobtype = ["xtb_opt"]
                self.set_general_setting("bhess", tmp)

                # put results into conformer objects
                for conf in self.confs_nc:
                    conf.results[self._name].update(results_rrho[conf.geom.id])
                    conf.results[self._name]["gtot"] = self.grrho(conf)

                # flag to make sure that rrho is only calculated once
                rrho_done = True

            # TODO - crestcheck each iteration if ncyc >= 6

            # sort conformers
            self.ensemble.conformers.sort(key=lambda conf: self.grrho(conf))

            threshold = self.get_settings()["threshold"] / AU2KCAL

            # threshold increase based on mean trajectory similarity
            if len(self.ensemble.conformers) > 1:
                trajectories = []
                for conf in self.ensemble.conformers:
                    # If conformers already converged before, use the constant final energy as trajectory
                    if conf not in self.confs_nc:
                        trajectories.append(
                            [
                                conf.results[self._name]["xtb_opt"]["energy"]
                                for _ in range(self.get_settings()["optcycles"])
                            ]
                        )
                    # If conformers converged in this macrocycle, use the trajectory from the results and pad it with
                    # the final energy
                    elif results_opt[conf.geom.id]["xtb_opt"]["converged"]:
                        trajectories.append(
                            results_opt[conf.geom.id]["xtb_opt"]["ecyc"] +
                            [
                                conf.results[self._name]["xtb_opt"]["energy"]
                                for _ in range(self.get_settings()["optcycles"] -
                                               len(results_opt[conf.geom.id]["xtb_opt"]["ecyc"]))
                            ]
                        )
                    # If conformers did not converge, use the trajectory from the results
                    else:
                        trajectories.append(
                            conf.results[self._name]["xtb_opt"]["ecyc"])

                mu_sim = mean_similarity(trajectories)
                threshold += (1.5 / AU2KCAL) * \
                    max(1 - exp(-AU2KCAL * mu_sim), 0.0)
                logger.debug(
                    f"Mean trajectory similarity: {AU2KCAL * mu_sim:.2f} kcal/mol"
                )

            logger.info(f"Threshold: {threshold * AU2KCAL:.2f} kcal/mol")

            # remove converged conformers from 'todo-list'
            for conf in list(
                filter(
                    lambda x: x.results[self._name]["xtb_opt"]["converged"],
                    self.confs_nc,
                )
            ):
                print(
                    f"{conf.name} converged after {ncyc + results_opt[conf.geom.id]['xtb_opt']['cycles']} steps."
                )
                self.confs_nc.remove(conf)

            # update the conformer list (remove conf if below threshold and gradient too small for all microcycles in
            # this macrocycle)
            self.ensemble.update_conformers(
                self.grrho,
                threshold,
                additional_filter=lambda x: all(
                    gn < self.get_settings()["gradthr"]
                    for gn in x.results[self._name]["xtb_opt"]["gncyc"]
                )
                # x.results[self._name]["xtb_opt"]["grad_norm"] < self.get_settings()["gradthr"]
            )

            # make sure that all the conformers, that are not converged but filtered out, are also removed
            # from self.confs_nc
            limit = min(self.grrho(conf)
                        for conf in self.ensemble.conformers)
            for conf in self.ensemble.rem:
                if conf in self.confs_nc:
                    print(
                        f"{conf.name} is no longer considered (gradient too small and"
                        f" ΔG = {(self.grrho(conf) - limit) * AU2KCAL:.2f})."
                    )
                    self.confs_nc.remove(conf)

            # TODO - print out information about current state of the ensemble
            self.print_update()

            # update number of cycles
            ncyc += self.get_settings()["optcycles"]
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
            f"% at {self.get_general_settings().get('temperature', 298.15)} K",
        ]

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(self.grrho(conf) for conf in self.ensemble.conformers)

        dftmin = min(
            conf.results[self._name]["xtb_opt"]["energy"]
            for conf in self.ensemble.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda conf: f"{conf.results[self._name]['xtb_opt']['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda conf: f"{(conf.results[self._name]['xtb_opt']['energy'] - dftmin) * AU2KCAL:.2f}",
            "GmRRHO": lambda conf: f"{conf.results[self._name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
            if self.get_general_settings()["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
        }

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # write lines to file
        logger.debug(
            f"Writing to {os.path.join(self.ensemble.workdir, f'{self._name}.out')}."
        )
        with open(
            os.path.join(self.ensemble.workdir, f"{self._name}.out"), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results of this part to a json file
        self.write_json()

    def print_update(self) -> None:
        """
        writes information about the current state of the ensemble
        """
        headers = [
            "CONF#",
            "E (DFT) (+ ΔGsolv)",
            "ΔE (DFT) (+ δΔGsolv)",
            "GmRRHO",
            "Gtot",
            "ΔGtot",
            "grad_norm",
        ]
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
            "[Eh/a0]",
        ]

        limit = min(conf.results[self._name]['xtb_opt']['energy']
                    for conf in self.ensemble.conformers)
        limit2 = min(self.grrho(conf) for conf in self.ensemble.conformers)

        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda conf: f"{conf.results[self._name]['xtb_opt']['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda conf: f"{(conf.results[self._name]['xtb_opt']['energy'] - limit) * AU2KCAL:.2f}",
            "GmRRHO": lambda conf: f"{conf.results[self._name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
            if self.get_general_settings()["evaluate_rrho"]
            else "---",
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - limit2) * AU2KCAL:.2f}",
            "grad_norm": lambda conf: f"{conf.results[self._name]['xtb_opt']['grad_norm']:.6f}",
        }
        rows = [
            printmap[header](conf) for header in headers
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        for line in lines:
            print(line)
