"""
Optimization == part2
performing geometry optimization of the CRE and provide low level free energies.
"""

import os

from .optimizer import EnsembleOptimizer
from ..ensembledata import EnsembleData
from ..datastructure import MoleculeData
from ..parallel import execute
from ..params import Params
from ..utilities import print, format_data, h1, DfaHelper
from ..logging import setup_logger

logger = setup_logger(__name__)


class Optimization(EnsembleOptimizer):
    __solv_mods = {
        prog: tuple(
            t for t in Params.SOLV_MODS[prog] if t not in ("cosmors", "cosmors-fine")
        )
        for prog in Params.PROGS
    }

    _grid = "high"

    _options = {
        "optcycles": {"default": 8},
        "maxcyc": {"default": 200},
        "threshold": {"default": 1.5},
        "hlow": {"default": 0.01},
        "gradthr": {"default": 0.01},
        "func": {
            "default": "r2scan-3c",
            "options": {prog: DfaHelper.get_funcs(prog) for prog in Params.PROGS},
        },
        "basis": {"default": "def2-TZVP"},
        "prog": {"default": "tm", "options": Params.PROGS},
        "sm": {"default": "dcosmors", "options": __solv_mods},
        "gfnv": {"default": "gfn2", "options": Params.GFNOPTIONS},
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
        "xtb_opt": {"default": True},
    }

    _settings = {}

    def _optimize(self, cut: bool = True) -> None:
        """
        Optimization of the ensemble at DFT level (possibly with implicit solvation)

        Uses xtb as driver for orca/tm, calculates hessians from orca/tm single-points and reevaluates ensemble every 'optcycles' steps
        by calculating the thermodynamics of the optimized geometries

        Alternatively just run the complete geometry optimization for every conformer with xtb as driver (decide with 'macrocycles')
        """
        # Special 'todo-list' for optimization part, contains all unconverged conformers,
        # used in macrocycle optimization
        self.confs_nc: list = None

        # Attribute to store path to constraints file if used
        self.constraints = None

        """
        don't use this for now, use new keyword 'maxcyc' instead
        if config.nat > 200:
            # stopcycle = don't optimize more than stopcycle cycles
            stopcycle = config.nat * 2
        else:
            stopcycle = 200
        """

        # Set jobtype depending on program (tm has not tm-pure geom opt yet)
        jobtype = None
        if self.get_settings()["prog"] == "orca":
            jobtype = ["xtb_opt"] if self.get_settings()["xtb_opt"] else ["opt"]
        else:
            if not self.get_settings()["xtb_opt"]:
                logger.warning(
                    "TURBOMOLE-driven geometry optimization not available yet. Switching back to ANCOPT as driver."
                )
            jobtype = ["xtb_opt"]
        assert jobtype is not None

        # NOTE: (IMPORTANT) the following only uses xtb as driver (no native geometry optimizations)
        # Check for constraint file
        if self.get_settings()["constrain"]:
            assert os.path.isfile(os.path.join(os.getcwd(), "constraints.xtb"))
            print("Found constraints-file constraints.xtb ...")
            self.constraints = os.path.join(os.getcwd(), "constraints.xtb")

        # Use macrocycle optimization only if there is more than one conformer
        if self.get_settings()["macrocycles"] and len(self.ensemble.conformers) > 1:
            # ensembleopt using macrocycles with 'optcycles' microcycles
            self.__macrocycle_opt(cut)
        else:
            # do complete geometry optimization
            if not len(self.ensemble.conformers) > 1:
                print(
                    f"Only one conformer ({self.ensemble.conformers[0].name}) is available for optimization."
                )

            # disable spearman optimization
            print("Macrocycle optimization turned off.")
            self.set_setting("macrocycles", False)

            prepinfo = self._setup_prepinfo(jobtype)

            # Run complete geometry optimization
            results_opt, failed = execute(
                self.ensemble.conformers,
                self.dir,
                self.get_settings()["prog"],
                prepinfo,
                jobtype,
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                omp=self.get_general_settings()["omp"],
                retry_failed=self.get_general_settings()["retry_failed"],
            )

            # Remove failed conformers
            self.ensemble.remove_conformers(failed)

            # update results for each conformer
            for conf in self.ensemble.conformers:
                self.results[conf.name].update(results_opt[conf.name])
                # update geometry of the conformer
                conf.geom.xyz = results_opt[conf.name][jobtype[0]]["geom"]

        # Handle unconverged conformers (WIP)
        unconverged = self.confs_nc or [
            conf.name
            for conf in self.ensemble.conformers
            if not self.results[conf.name][jobtype[0]]["converged"]
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
            for confname in unconverged:
                print(
                    f"{confname}, grad_norm: {self.results[confname][jobtype[0]]['grad_norm']}"
                )
            print("The unconverged conformers will now be removed from consideration.")
            self.ensemble.remove_conformers(unconverged)

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?).
        # we don't do that anymore and just use the final energies from the optimizations, which are done using a
        # solvent model, since this basically makes no difference in comp time
        # do rrho on converged geometries (overwrites previous rrho calculations)
        jobtype = ["xtb_rrho"]
        prepinfo = self._setup_prepinfo(jobtype)
        results, failed = execute(
            self.ensemble.conformers,
            self.dir,
            self.get_settings()["prog"],
            prepinfo,
            jobtype,
            copy_mo=self.get_general_settings()["copy_mo"],
            balance=self.get_general_settings()["balance"],
            omp=self.get_general_settings()["omp"],
            retry_failed=self.get_general_settings()["retry_failed"],
        )

        # Remove failed conformers
        self.ensemble.remove_conformers(failed)

        # TODO - Add the possibility to explicitely calculate solvation contributions

        # Update results
        for conf in self.ensemble.conformers:
            self.results[conf.name].update(results[conf.name])
            self.results[conf.name]["gtot"] = self._grrho(conf)

        # sort conformers list with optimization key (gtot)
        self.ensemble.conformers.sort(key=lambda conf: self.results[conf.name]["gtot"])

        # calculate boltzmann weights from gtot values calculated here
        self._calc_boltzmannweights()

        # write final results
        self._write_results()

    def _grrho(self, conf: MoleculeData) -> float:
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
                self.results[conf.name][jobtype]["energy"]
                + self.results[conf.name]["xtb_rrho"]["energy"]
            )
        except KeyError:
            return self.results[conf.name][jobtype]["energy"]

    def __macrocycle_opt(self, cut: bool):
        """
        Macrocycle optimization using xtb as driver. Also calculates GmRRHO for finite temperature contributions
        and uses adaptive threshold based on mean trajectory similarity.
        """
        # make a separate list of conformers that only includes (considered) conformers that are not converged
        # NOTE: this is a special step only necessary for macrocycle optimization
        # at this point it's just self.ensemble.conformers, it is basically a todo-list
        self.confs_nc = [conf.name for conf in self.ensemble.conformers]

        jobtype = ["xtb_opt"] if self.get_settings()["xtb_opt"] else ["opt"]

        prepinfo = self._setup_prepinfo(jobtype)

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
            results_opt, failed = execute(
                [
                    conf
                    for conf in self.ensemble.conformers
                    if conf.name in self.confs_nc
                ],
                self.dir,
                self.get_settings()["prog"],
                prepinfo,
                jobtype,
                copy_mo=self.get_general_settings()["copy_mo"],
                balance=self.get_general_settings()["balance"],
                omp=self.get_general_settings()["omp"],
                retry_failed=self.get_general_settings()["retry_failed"],
            )

            # Remove failed conformers
            self.ensemble.remove_conformers(failed)
            for conf in failed:
                self.confs_nc.remove(conf.name)

            # put geometry optimization results into conformer objects
            for conf in filter(
                lambda x: x.name in self.confs_nc, self.ensemble.conformers
            ):
                # update geometry of the conformer
                conf.geom.xyz = results_opt[conf.name][jobtype[0]]["geom"]

                self.results.setdefault(conf.name, {}).setdefault(jobtype[0], {})

                # Update the values for "energy", "grad_norm", "converged", "geom"
                # NOTE: this replaces the default self.results.update, why see below
                for key in ["energy", "grad_norm", "converged", "geom"]:
                    self.results[conf.name][jobtype[0]][key] = results_opt[conf.name][
                        jobtype[0]
                    ][key]

                # NOTE: Because of the following two steps it is necessary not to blindly call self.results.update,
                # since values would be overwritten and information lost
                # Add the number of cycles
                self.results[conf.name][jobtype[0]].setdefault("cycles", 0)
                self.results[conf.name][jobtype[0]]["cycles"] += results_opt[conf.name][
                    jobtype[0]
                ]["cycles"]

                # Extend the energy and grad_norm lists
                for key in ["ecyc", "gncyc"]:
                    self.results[conf.name][jobtype[0]].setdefault(key, []).extend(
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
                results, failed = execute(
                    self.confs_nc,
                    self.dir,
                    self.get_settings()["prog"],
                    self._setup_prepinfo(jobtype),
                    jobtype,
                    copy_mo=self.get_general_settings()["copy_mo"],
                    balance=self.get_general_settings()["balance"],
                    omp=self.get_general_settings()["omp"],
                    retry_failed=self.get_general_settings()["retry_failed"],
                )

                # Remove failed conformers
                self.ensemble.remove_conformers(failed)

                # Reset
                jobtype = tmp2
                self.set_general_setting("bhess", tmp)

                # Update results
                for conf in filter(
                    lambda x: x.name in self.confs_nc, self.ensemble.conformers
                ):
                    self.results[conf.name].update(results[conf.name])
                    self.results[conf.name]["gtot"] = self._grrho(conf)

                # flag to make sure that rrho is only calculated once
                rrho_done = True

            # TODO - crestcheck each iteration if ncyc >= 6
            # => replace this with molbar later

            # sort conformers
            self.ensemble.conformers.sort(self._grrho)

            # remove converged conformers from 'todo-list'
            for confname in list(
                filter(
                    lambda x: self.results[x][jobtype[0]]["converged"],
                    self.confs_nc,
                )
            ):
                print(
                    f"{confname} converged after {ncyc + results_opt[confname][jobtype[0]]['cycles']} steps."
                )
                self.confs_nc.remove(confname)
                nconv += 1

            if cut:
                threshold = self.get_settings()["threshold"] / Params.AU2KCAL

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
                        / Params.AU2KCAL
                    )

                logger.info(f"Threshold: {threshold * Params.AU2KCAL:.2f} kcal/mol")

                # TODO - count removed conformers as converged?
                # update the conformer list (remove conf if below threshold and gradient too small for all microcycles in
                # this macrocycle)
                # NOTE: we need to look at results_opt here to look at only the current macrocycle
                limit = min(self._grrho(conf) for conf in self.ensemble.conformers)
                filtered = list(
                    filter(
                        lambda conf: all(
                            gn < self.get_settings()["gradthr"]
                            for gn in results_opt[conf.name][jobtype[0]]["gncyc"]
                        )
                        and self._grrho(conf) - limit < threshold,
                        self.ensemble.conformers,
                    )
                )
                self.ensemble.remove_conformers([conf.name for conf in filtered])
                for conf in filtered:
                    # print(f"No longer considering {conf}.")
                    print(
                        f"No longer considering {conf.name} (gradient too small and"
                        f" ΔG = {(self._grrho(conf) - limit) * Params.AU2KCAL:.2f})."
                    )
                    # make sure that all the conformers, that are not converged but filtered out, are also removed
                    # from self.confs_nc
                    self.confs_nc.remove(conf.name)

            # update number of cycles
            ncyc += self.get_settings()["optcycles"]

            # Print out information about current state of the ensemble
            self.__print_opt_update(ncyc)

    def _write_results(self) -> None:
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
        gtotmin = min(self._grrho(conf) for conf in self.ensemble.conformers)

        # Minimal pure DFT energy
        dftmin = min(
            self.results[conf.name][jobtype]["energy"]
            for conf in self.ensemble.conformers
        )

        # Define what gets printed for which header
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (DFT) (+ ΔGsolv)": lambda conf: f"{self.results[conf.name][jobtype]['energy']:.6f}",
            "ΔE (DFT) (+ δΔGsolv)": lambda conf: f"{(self.results[conf.name][jobtype]['energy'] - dftmin) * Params.AU2KCAL:.2f}",
            "GmRRHO": lambda conf: (
                f"{self.results[conf.name]['xtb_rrho']['gibbs'][self.get_general_settings()['temperature']]:.6f}"
                if self.get_general_settings()["evaluate_rrho"]
                else "---"
            ),
            "Gtot": lambda conf: f"{self._grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self._grrho(conf) - gtotmin) * Params.AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{self.results[conf.name]['bmw'] * 100:.2f}",
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
                self.results[conf.name]["bmw"] * self.results[conf.name]["gtot"]
                for conf in self.ensemble.conformers
            ]
        )

        # calculate averaged free energy
        avE = sum(
            [
                self.results[conf.name]["bmw"]
                * self.results[conf.name][jobtype]["energy"]
                for conf in self.ensemble.conformers
            ]
        )

        # append the lines for the free energy/enthalpy
        lines.append(
            f"{self.get_general_settings().get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part2==\n"
        )
        lines.append("".ljust(int(Params.PLENGTH), "-") + "\n\n")

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        # write lines to file
        filename = f"{self._part_nos[self._name]}_{self._name.upper()}.out"
        logger.debug(f"Writing to {os.path.join(os.getcwd(), filename)}.")
        with open(os.path.join(os.getcwd(), filename), "w", newline=None) as outfile:
            outfile.writelines(lines)

        # Additionally, write the results of this part to a json file
        self._write_json()

    def __print_opt_update(self, ncyc) -> None:
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
        limit = min(self._grrho(conf) for conf in self.ensemble.conformers)

        # Define what gets printed for which header
        printmap = {
            "CONF#": lambda conf: conf.name,
            "Gtot": lambda conf: f"{self._grrho(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self._grrho(conf) - limit) * Params.AU2KCAL:.2f}",
            "grad_norm": lambda conf: f"{self.results[conf.name][jobtype]['grad_norm']:.6f}",
            "converged": lambda conf: f"{self.results[conf.name][jobtype]['converged']}",
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

        print("".ljust(Params.PLENGTH, "-"))
