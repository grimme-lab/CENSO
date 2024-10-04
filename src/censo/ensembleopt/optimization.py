"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""

import os
from functools import reduce

from .optimizer import EnsembleOptimizer
from ..ensembledata import EnsembleData
from ..datastructure import MoleculeData
from ..parallel import execute
from ..params import SOLV_MODS, PROGS, GFNOPTIONS, AU2KCAL, PLENGTH
from ..utilities import print, format_data, h1
from ..logging import setup_logger

logger = setup_logger(__name__)


class Optimization(EnsembleOptimizer):
    _part_no = "2"

    __solv_mods = tuple(
        t
        for t in reduce(lambda x, y: x + y, SOLV_MODS.values())
        if t not in ("cosmors", "cosmors-fine")
    )

    _grid = "high"

    _options = {
        "optcycles": {"default": 8},
        "maxcyc": {"default": 200},
        "threshold": {"default": 1.5},
        "hlow": {"default": 0.01},
        "gradthr": {"default": 0.01},
        "func": {"default": "r2scan-3c"},
        "basis": {"default": "def2-TZVP"},
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
        "run": {"default": True},
        "macrocycles": {"default": True},
        "crestcheck": {"default": False},
        "template": {"default": False},
        "constrain": {"default": False},
        "xtb_opt": {"default": False},
    }

    _settings = {}

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

        # Special 'todo-list' for optimization part, contains all unconverged conformers,
        # used in macrocycle optimization
        self.confs_nc: list[MoleculeData] = None

        # Attribute to store path to constraints file if used
        self.constraints = None

    def optimize(self, ncores: int, cut: bool = True) -> None:
        """
        Optimization of the ensemble at DFT level (possibly with implicit solvation)

        Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
        by calculating the thermodynamics of the optimized geometries

        Alternatively just run the complete geometry optimization for every conformer with xtb as driver (decide with 'macrocycles')
        """
        """
        don't use this for now, use new keyword 'maxcyc' instead
        if config.nat > 200:
            # stopcycle = don't optimize more than stopcycle cycles
            stopcycle = config.nat * 2
        else:
            stopcycle = 200
        """
        jobtype = ["xtb_opt"] if self.get_settings()["xtb_opt"] else ["opt"]
        # NOTE: (IMPORTANT) the following only uses xtb as driver (no native geometry optimizations)
        # Check for constraint file
        if self.get_settings()["constrain"]:
            assert os.path.isfile(
                os.path.join(self.ensemble.workdir, "constraints.xtb")
            )
            print("Found constraints-file constraints.xtb ...")
            self.constraints = os.path.join(self.ensemble.workdir, "constraints.xtb")

        # Use macrocycle optimization only if there is more than one conformer
        if self.get_settings()["macrocycles"] and len(self.ensemble.conformers) > 1:
            # ensembleopt using macrocycles with 'optcycles' microcycles
            self.__macrocycle_opt(ncores, cut)
        else:
            # do complete geometry optimization
            if not len(self.ensemble.conformers) > 1:
                print(
                    f"Only one conformer ({self.ensemble.conformers[0].name}) is available for optimization."
                )

            # disable spearman optimization
            print("Macrocycle optimization turned off.")
            self.set_setting("macrocycles", False)

            prepinfo = self.setup_prepinfo(jobtype)

            # run optimizations using xtb as driver
            success, results_opt, failed = execute(
                self.ensemble.conformers,
                self.dir,
                self.get_settings()["prog"],
                prepinfo,
                jobtype,
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                omp=self.get_general_settings()["omp"],
                maxcores=ncores,
                retry_failed=self.get_general_settings()["retry_failed"],
            )

            # Remove failed conformers
            self.ensemble.remove_conformers(failed)

            # update results for each conformer
            for conf in self.ensemble.conformers:
                # update geometry of the conformer
                conf.geom.xyz = results_opt[conf.name][jobtype[0]]["geom"]

        # Handle unconverged conformers (TODO)
        unconverged = self.confs_nc or [
            conf
            for conf in self.ensemble.conformers
            if not conf.results[self._name][jobtype[0]]["converged"]
        ]

        if len(unconverged) == len(self.ensemble.conformers):
            # TODO - important - is this somehow recoverable?
            raise RuntimeError(
                "Cannot continue because there is no conformer with converged geometry optimization."
            )
        elif len(unconverged) > 0:
            logger.warning(
                "The geometry optimization of at least one conformer did not converge."
            )
            print("Unconverged conformers:")
            for conf in unconverged:
                print(
                    f"{conf.name}, grad_norm: {conf.results[self._name][jobtype[0]]['grad_norm']}"
                )
            print("The unconverged conformers will now be removed from consideration.")
            self.ensemble.remove_conformers([conf.name for conf in unconverged])

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?).
        # we don't do that anymore and just use the final energies from the optimizations, which are done using a
        # solvent model, since this basically makes no difference in comp time
        # do rrho on converged geometries (overwrites previous rrho calculations)
        jobtype = ["xtb_rrho"]
        prepinfo = self.setup_prepinfo(jobtype)
        success, _, failed = execute(
            self.ensemble.conformers,
            self.dir,
            self.get_settings()["prog"],
            prepinfo,
            jobtype,
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            omp=self.get_general_settings()["omp"],
            maxcores=ncores,
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

        for conf in self.ensemble.conformers:
            conf.results[self._name]["gtot"] = self.grrho(conf)

        # sort conformers list with optimization key (gtot)
        self.ensemble.conformers.sort(key=lambda conf: conf.results[self._name]["gtot"])

        # calculate boltzmann weights from gtot values calculated here
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings().get("temperature", 298.15), self._name
        )

        # write final results
        self.write_results()

    def grrho(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy (last step of running optimization) and Gmrrho
        If no GmRRHO is calculated only the most recent DFT energy is returned

        Args:
            conf: conformer object

        Returns:
            float: Gtot (in the unit of the results)
        """
        jobtype = "xtb_opt" if self.get_settings()["xtb_opt"] else "opt"

        try:
            return (
                conf.results[self._name][jobtype]["energy"]
                + conf.results[self._name]["xtb_rrho"]["energy"]
            )
        except KeyError:
            return conf.results[self._name][jobtype]["energy"]

    def __macrocycle_opt(self, ncores: int, cut: bool):
        """
        Macrocycle optimization using xtb as driver. Also calculates GmRRHO for finite temperature contributions
        and uses adaptive threshold based on mean trajectory similarity.
        """
        # make a separate list of conformers that only includes (considered) conformers that are not converged
        # NOTE: this is a special step only necessary for macrocycle optimization
        # at this point it's just self.ensemble.conformers, it is basically a todo-list
        self.confs_nc = self.ensemble.conformers.copy()

        jobtype = ["xtb_opt"] if self.get_settings()["xtb_opt"] else ["opt"]

        prepinfo = self.setup_prepinfo(jobtype)

        ncyc = 0
        rrho_done = False
        print(
            f"Optimization using macrocycles, {self.get_settings()['optcycles']} microcycles per step."
        )
        nconv = 0
        ninit = len(self.confs_nc)
        while len(self.confs_nc) > 0 and ncyc < self.get_settings()["maxcyc"]:
            # NOTE: this loop works through confs_nc, so if the geometry optimization for a conformer is converged,
            # and therefore removed from our 'todo-list', all the following steps will not consider it anymore
            # run optimizations for 'optcycles' steps
            success, results_opt, failed = execute(
                self.confs_nc,
                self.dir,
                self.get_settings()["prog"],
                prepinfo,
                jobtype,
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                omp=self.get_general_settings()["omp"],
                maxcores=ncores,
                retry_failed=self.get_general_settings()["retry_failed"],
                update=False,
            )

            # Remove failed conformers
            self.ensemble.remove_conformers(failed)
            for conf in filter(lambda x: x.name in failed, self.confs_nc):
                self.confs_nc.remove(conf)

            # put geometry optimization results into conformer objects
            for conf in self.confs_nc:
                # update geometry of the conformer
                conf.geom.xyz = results_opt[conf.name][jobtype[0]]["geom"]

                conf.results.setdefault(self._name, {}).setdefault(jobtype[0], {})

                # Update the values for "energy", "grad_norm", "converged", "geom"
                for key in ["energy", "grad_norm", "converged", "geom"]:
                    conf.results[self._name][jobtype[0]][key] = results_opt[conf.name][
                        jobtype[0]
                    ][key]

                # Add the number of cycles
                conf.results[self._name][jobtype[0]].setdefault("cycles", 0)
                conf.results[self._name][jobtype[0]]["cycles"] += results_opt[
                    conf.name
                ][jobtype[0]]["cycles"]

                # Extend the energy and grad_norm lists
                for key in ["ecyc", "gncyc"]:
                    conf.results[self._name][jobtype[0]].setdefault(key, []).extend(
                        results_opt[conf.name][jobtype[0]][key]
                    )

            # run xtb_rrho for finite temperature contributions
            # for now only after the first 'optcycles' steps or after at least 6 cycles are done
            # TODO - make this better
            if ncyc + self.get_settings()["optcycles"] >= 6 and not rrho_done:
                # evaluate rrho using bhess flag (reset after this)
                # TODO - this smells a bit
                tmp = self.get_general_settings()["bhess"]
                tmp2 = jobtype
                self.set_general_setting("bhess", True)
                jobtype = ["xtb_rrho"]
                success, _, failed = execute(
                    self.confs_nc,
                    self.dir,
                    self.get_settings()["prog"],
                    self.setup_prepinfo(jobtype),
                    jobtype,
                    copy_mo=self.get_general_settings()["copy_mo"],
                    balance=self.get_general_settings()["balance"],
                    omp=self.get_general_settings()["omp"],
                    maxcores=ncores,
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # TODO - what to do about failed conformers?

                # Reset
                jobtype = tmp2
                self.set_general_setting("bhess", tmp)

                # put results into conformer objects
                for conf in self.confs_nc:
                    conf.results[self._name]["gtot"] = self.grrho(conf)

                # flag to make sure that rrho is only calculated once
                rrho_done = True

            # TODO - crestcheck each iteration if ncyc >= 6

            # sort conformers
            self.ensemble.conformers.sort(key=lambda conf: self.grrho(conf))

            # remove converged conformers from 'todo-list'
            for conf in list(
                filter(
                    lambda x: x.results[self._name][jobtype[0]]["converged"],
                    self.confs_nc,
                )
            ):
                print(
                    f"{conf.name} converged after {ncyc + results_opt[conf.name][jobtype[0]]['cycles']} steps."
                )
                self.confs_nc.remove(conf)
                nconv += 1

            if cut:
                threshold = self.get_settings()["threshold"] / AU2KCAL

                # threshold increase based on number of converged conformers
                # NOTE: it is important to work with the results of the current macrocycle,
                # since in the results dict of the MoleculeData objects all the results
                # from previous cycles are stored
                if len(self.ensemble.conformers) > 1:
                    n = 1
                    threshold += (
                        n
                        * (
                            self.get_settings()["threshold"]
                            - self.get_settings()["threshold"] * nconv / ninit
                        )
                        / AU2KCAL
                    )

                logger.info(f"Threshold: {threshold * AU2KCAL:.2f} kcal/mol")

                # TODO - count removed conformers as converged?
                # update the conformer list (remove conf if below threshold and gradient too small for all microcycles in
                # this macrocycle)
                for conf in self.ensemble.update_conformers(
                    self.grrho,
                    threshold,
                    additional_filter=lambda x: all(
                        gn < self.get_settings()["gradthr"]
                        for gn in x.results[self._name][jobtype[0]]["gncyc"]
                    ),
                    # x.results[self._name]["xtb_opt"]["grad_norm"] < self.get_settings()["gradthr"]
                ):
                    print(f"No longer considering {conf}.")

                # make sure that all the conformers, that are not converged but filtered out, are also removed
                # from self.confs_nc
                limit = min(self.grrho(conf) for conf in self.ensemble.conformers)
                for conf in self.ensemble.rem:
                    if conf in self.confs_nc:
                        print(
                            f"{conf.name} is no longer considered (gradient too small and"
                            f" ΔG = {(self.grrho(conf) - limit)
                                     * AU2KCAL:.2f})."
                        )
                        self.confs_nc.remove(conf)

            # update number of cycles
            ncyc += self.get_settings()["optcycles"]

            # Print out information about current state of the ensemble
            self.print_opt_update(ncyc)

    def write_results(self) -> None:
        """
        formatted write of part results (optional)
        """
        print(h1(f"{self._name.upper()} RESULTS"))
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
        jobtype = "xtb_opt" if self.get_settings()["xtb_opt"] else "opt"

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(self.grrho(conf) for conf in self.ensemble.conformers)

        # Minimal pure DFT energy
        dftmin = min(
            conf.results[self._name][jobtype]["energy"]
            for conf in self.ensemble.conformers
        )

        # Define what gets printed for which header
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda conf: f"{conf.results[self._name][jobtype]['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda conf: f"{(conf.results[self._name][jobtype]
                ['energy'] - dftmin) * AU2KCAL:.2f}",
            "GmRRHO": lambda conf: (
                f"{conf.results[self._name]['xtb_rrho']['gibbs']
                [self.get_general_settings()['temperature']]:.6f}"
                if self.get_general_settings()["evaluate_rrho"]
                else "---"
            ),
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self._name]['bmw'] * 100:.2f}",
        }

        # Create rows via the printmap
        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        # Format everything into a table
        lines = format_data(headers, rows, units=units)

        # list the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble on optimized geometries:\n"
        )
        lines.append(
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n"
        )

        # calculate averaged free enthalpy
        avG = sum(
            [
                conf.results[self._name]["bmw"] * conf.results[self._name]["gtot"]
                for conf in self.ensemble.conformers
            ]
        )

        # calculate averaged free energy
        avE = sum(
            [
                conf.results[self._name]["bmw"]
                * conf.results[self._name][jobtype]["energy"]
                for conf in self.ensemble.conformers
            ]
        )

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part2==\n"
        )
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # write lines to file
        filename = f"{self._part_no}_{self._name.upper()}.out"
        logger.debug(f"Writing to {os.path.join(self.ensemble.workdir, filename)}.")
        with open(
            os.path.join(self.ensemble.workdir, filename), "w", newline=None
        ) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results of this part to a json file
        self.write_json()

    def print_opt_update(self, ncyc) -> None:
        """
        writes information about the current state of the ensemble in form of a table
        """
        print(h1(f"{self._name.upper()} CYCLE {ncyc} UPDATE"))
        # Define headers for the table
        headers = [
            "CONF#",
            "Gtot",
            "ΔGtot",
            "grad_norm",
            "converged",
        ]

        # Define the units for the table
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh/a0]",
            "",
        ]
        jobtype = "xtb_opt" if self.get_settings()["xtb_opt"] else "opt"

        # Lower limit for the free enthalpy
        limit = min(self.grrho(conf) for conf in self.ensemble.conformers)

        # Define what gets printed for which header
        printmap = {
            "CONF#": lambda conf: conf.name,
            "Gtot": lambda conf: f"{self.grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.grrho(conf) - limit) * AU2KCAL:.2f}",
            "grad_norm": lambda conf: f"{conf.results[self._name][jobtype]['grad_norm']:.6f}",
            "converged": lambda conf: f"{conf.results[self._name][jobtype]['converged']}",
        }

        # Create rows via the printmap
        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        # Format everything into a table
        lines = format_data(headers, rows, units=units)

        # Print the lines
        for line in lines:
            print(line, end="")

        print("".ljust(PLENGTH, "-"))
