"""
Performs the parallel execution of the QM calls.
"""
from functools import reduce
import os
from typing import Any, Dict, List
from concurrent.futures import ProcessPoolExecutor
from pprint import pprint

from censo.procfact import ProcessorFactory
from censo.utilities import print
from censo.datastructure import GeometryData, MoleculeData
from censo.settings import CensoSettings


class ProcessHandler:

    def __init__(self, settings: CensoSettings, conformers: List[GeometryData]):
        """
        creates instance with 
        also sets the number of processes (nprocs) and cores per processes (omp) to optimized vaules if possible
        """
        
        self.settings = settings

        # stores the conformers
        self.conformers: List[GeometryData] = conformers
        
        # stores the processor
        self.__processor = None
        
        # get number of cores
        try:
            self.__ncores = os.cpu_count()
        except AttributeError:
            raise AttributeError("Could not determine number of cores.")
        
        self.__nprocs = None
        self.__omp = None
        if (
            getattr(self.settings.settings_current.byname("balance"), "value", False) 
            or not self.__balance_load()
        ):
            print("\nNot readjusting maxprocs and omp.\n")
            # get number of processes
            self.__nprocs = getattr(
                self.settings.settings_current.byname("maxprocs"), 
                "value",
                None
            )
            
            # get number of cores per process
            self.__omp = getattr(
                self.settings.settings_current.byname("omp"),
                "value",
                None
            )
    
    
    def __balance_load(self) -> bool:
        """
        distribute computational load between cores
        keeping the number of processes below number of conformers
        trying to utilize as many cores as possible
        """
        # check if load can be rebalanced
        if (
            not self.__ncores is None 
            and not (
                self.__nprocs is None or self.__omp is None
            )
            and (
                self.__ncores > len(self.conformers) 
                or self.__ncores > self.__nprocs * self.__omp
            )
        ):
            # you want to create a number of parallel processes that in the best case is a divisor of nconf
            # the number of processes should allow for utilization of all cores
            try:
                # finds the largest divisor of ncores that is less or equal than nconf
                self.__nprocs = max([i for i in range(self.__ncores) if self.__ncores % i and i <= len(self.conformers)])
            except ValueError:
                print("There was an error while determining the number of processes in load balancing.") # TODO
                return False
            
            # should always be divisible, floor division only as safeguard
            self.__omp = self.__ncores // self.__nprocs
            
            return True
        else:
            return False
                

    def execute(self, jobtype: List[str], instructions: Dict[str, Any], workdir: str):
        """
        creates and executes the processes
        returns the results sorted by conformer, divided into jobtypes

        jobtype: list of ordered jobtypes (e.g. [xtb_sp, xtb_gsolv])
        instructions: dict with settings from CensoSettings (specific for part, see CensoPart)
        workdir: absolute path to folder where calculations should be executed in
        """
        # TODO - 'smart balancing'
        # set cores per process in instructions
        instructions["omp"] = self.__omp

        # try to get program from instructions
        prog = instructions.get("prog", None)
        
        if prog is None:
            raise Exception # TODO
        
        # initialize the processor for the respective program (depends on part)
        # and set the jobtype as well as instructions, also pass workdir to compute in
        # also pass a lookup dict to the processor so it can set the solvent for the program call correctly
        self.__processor = ProcessorFactory.create_processor(
            prog, 
            self.settings.external_paths,
            solvents_dict=self.settings.solvents_db[instructions["solvent"]], 
            dfa_settings=self.settings.dfa_settings
        )
        self.__processor.jobtype = jobtype
        self.__processor.instructions = instructions
        self.__processor.workdir = workdir

        # execute processes for conformers
        # TODO - set PARNODES
        with ProcessPoolExecutor(max_workers=self.__nprocs) as executor:
            resiter = executor.map(self.__processor.run, self.conformers)
         
        # returns merged result dicts
        # structure of results: 
        #   e.g. {id(conf): {"xtb_sp": {"success": ..., "energy": ...}, ...}, ...}
        results = reduce(lambda x, y: {**x, **y}, resiter)
        pprint(results)
        
        # assert that there is a result for every conformer
        try:
            assert all([id(conf) in results.keys() for conf in self.conformers])  
        except AssertionError:
            raise KeyError("There is a mismatch between conformer ids and returned results. Cannot find at least one conformer id in results.")
        
        return results