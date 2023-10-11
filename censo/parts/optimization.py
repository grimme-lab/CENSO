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
        # NOTE: (IMPORTANT) this does not work for standard optimization (opt) yet!!!

        handler = ProcessHandler(self._instructions)

        self.print_info()

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
            # NOTE: this is a special step only necessary for spearman threshold optimization
            # at this point it's just self.core.conformers
            self.confs_nc = [conf.geom for conf in self.core.conformers]

            stopcond_converged = False
            ncyc = 0
            rrho_done = False
            print(f"Optimization using Spearman threshold, {self._instructions['optcycles']} cycles per step.")
            while (
                not stopcond_converged 
                and ncyc < self._instructions["maxcyc"]
                # TODO - maybe make this more intelligent: 
                # make maxcyc lower and if some apparently relevant conformer doesn't converge within it's chunk, 
                # move it to a new chunk and calculate later
            ):
                print(f"NCYC: {ncyc}")
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
                                coreconf.resutts[self.__class__.__name__.lower()]["gtot"] = self.key(conf)
                                break

                    rrho_done = True

                # TODO - crestcheck each iteration if ncyc >= 6

                # TODO - recalculate fuzzy threshold

                # kick out conformers above threshold
                threshold = self._instructions.get("threshold", None)
                if not threshold is None:
                    # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
                    limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
                    
                    # filter out all conformers above threshold and with a gradient norm smaller than 'gradthr'
                    # so that 'filtered' contains all conformers that should not be considered any further
                    f = lambda x: self.key(x) > limit + threshold and x.results[self.__class__.__name__.lower()]["xtb_opt"]["grad_norm"] < self._instructions["gradthr"], 
                    filtered = [
                        conf for conf in filter(
                            f,
                            self.core.conformers
                        )
                    ]
                    
                    # update the conformer list in core (remove conf if below threshold)
                    self.core.update_conformers(filtered)  

                    # also remove conformers from confs_nc
                    for conf in filtered:
                        self.confs_nc.remove(conf.geom)
                else:
                    # TODO - error handling
                    raise RuntimeError

                # update list of converged conformers
                for conf in self.confs_nc:
                    if results_opt[conf.id]["xtb_opt"]["converged"]:
                        self.confs_nc.remove(conf)

                # check if all (considered) conformers converged - End of While-loop
                stopcond_converged = len(self.confs_nc) == 0

                # TODO - print out information about current state of the ensemble
                self.write_update()
        else:
            """do normal optimization - not implemented yet"""
            # TODO
            pass

        # NOTE: old censo did a single-point after all optimizations were done (to include gsolv?). 
        # we don't do that anymore and just use the final energies from the optimizations, which are done using solventmodel,
        # since this basically makes no difference in comp time
        # TODO - do rrho on converged geometries
        handler.conformers = self.core.conformers
        results_rrho = handler.execute(["xtb_rrho"], os.path.join(folder, "rrho"))

        for conf in self.core.conformers:
            conf.results[self.__class__.__name__.lower()].update(results_rrho[id(conf)])

        # TODO - boltzmann calculations
        self.core.calc_boltzmannweights(
            self._instructions.get("temperature", 298.15),
            self.__class__.__name__.lower()
        )

        # TODO - write final results
        self.write_results()


    def key(self, conf: MoleculeData) -> float:
        """
        Calculate Gtot from DFT energy (last step of running optimization) and Gmrrho
        """
        # TODO - implement branch for standard optimization (no ancopt)
        return conf.results[self.__class__.__name__.lower()]["xtb_opt"]["energy"] + conf.results[self.__class__.__name__.lower()]["xtb_rrho"]["energy"]


    def write_results(self) -> None:
        """
        formatted write of part results (optional)
        """
        # TODO
        print(f"Conformers {[conf.name for conf in self.core.conformers]} remain.")


    def write_update(self) -> None:
        """
        writes information about the current state of the ensemble
        """
        # TODO
        for conf in self.confs_nc:
            print(f"{conf.name}: {conf.results[self.__class__.__name__.lower()]['xtb_opt']['energy']}")