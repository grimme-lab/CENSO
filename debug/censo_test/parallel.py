"""
Performs the parallel execution of the QM calls.
"""
import os
from typing import Any, Dict, List
from concurrent.futures import ProcessPoolExecutor

from censo_test.qm_processor import ProcessorFactory
from censo_test.utilities import print
from censo_test.datastructure import GeometryData, MoleculeData
from censo_test.settings import CensoSettings


class ProcessHandler:

    def __init__(self, settings: CensoSettings, conformers: List[GeometryData]):
        """
        creates instance with 
        also sets the number of processes (nprocs) and cores per processes (omp) to optimized vaules if possible
        """
        # stores the conformers
        self.conformers: List[GeometryData] = conformers
        
        # stores the processor
        self._processor = None
        
        # get number of cores
        try:
            self._ncores = os.cpu_count()
        except AttributeError:
            raise AttributeError("Could not determine number of cores.")
        
        self._nprocs = None
        self._omp = None
        if (
            getattr(settings.settings_current.byname("balance"), "value", False) 
            or not self._balance_load()
        ):
            print("\nCould not readjust maxprocs and omp.\n")
            # get number of processes
            self._nprocs = getattr(
                settings.settings_current.byname("maxprocs"), 
                "value",
                None
            )
            
            # get number of cores per process
            self._omp = getattr(
                settings.settings_current.byname("omp"),
                "value",
                None
            )
    
    
    def _balance_load(self) -> bool:
        """
        distribute computational load between cores
        keeping the number of processes below number of conformers
        trying to utilize as many cores as possible
        """
        # check if load can be rebalanced
        if (
            not self._ncores is None 
            and not (
                self._nprocs is None or self._omp is None
            )
            and (
                self._ncores > len(self.conformers) 
                or self._ncores > self._nprocs * self._omp
            )
        ):
            # you want to create a number of parallel processes that in the best case is a divisor of nconf
            # the number of processes should allow for utilization of all cores
            try:
                # finds the largest divisor of ncores that is less or equal than nconf
                self._nprocs = max([i for i in range(self._ncores) if self._ncores % i and i <= len(self.conformers)])
            except ValueError:
                print("There was an error while determining the number of processes in load balancing.") # TODO
                return False
            
            # should always be divisible, int casting only as safeguard
            self._omp = int(self._ncores / self._nprocs)
            
            return True
        else:
            return False
                

    def execute(self, jobtype: List[str], instructions: Dict[str, Any]):
        """
        creates and executes the processes
        returns the results sorted by conformer, divided into jobtypes
        """
        # TODO - 'smart balancing'
        # try to get program from instructions
        prog = instructions.get("prog", None)
        
        if prog is None:
            raise Exception # TODO
        
        # initialize the processor for the respective program (depends on part)
        # and set the jobtype as well as instructions
        self._processor = ProcessorFactory.create_processor(prog)
        self._processor.jobtype = jobtype
        self._processor.instructions = instructions
        
        results: Dict[int, Dict] = {}
        
        # execute processes for conformers
        # TODO - set PARNODES
        with ProcessPoolExecutor(max_workers=self._nprocs) as executor:
            resiter = executor.map(self._processor.run, self.conformers)

            # sort results by conformer and then by jobtype
            for res in resiter:
                if res["confid"] not in results.keys():
                    results[res["confid"]] = {}
                    
                results[res["confid"]]["success"] = True
                for job in jobtype:
                    results[res["confid"]][job] = res["confid"][job]
                    if not res["confid"]["success"]:
                        """job failed""" # TODO
                        results[res["confid"]]["success"] = False
            
        return results