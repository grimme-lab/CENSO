"""
Performs the parallel execution of the QM calls.
"""
from functools import reduce
import os
from typing import Any, Dict, List, Tuple
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
            # TODO - error handling
            raise AttributeError("Could not determine number of cores.")
        
        # get number of processes
        self.__nprocs = getattr(
            self.settings.settings_current.byname("maxprocs"), 
            "value"
        )
        
        # get number of cores per process
        self.__omp = getattr(
            self.settings.settings_current.byname("omp"),
            "value"
        )
                

    def execute(self, jobtype: List[str], instructions: Dict[str, Any], workdir: str):
        """
        creates and executes the processes
        returns the results sorted by conformer, divided into jobtypes

        jobtype: list of ordered jobtypes (e.g. [xtb_sp, xtb_gsolv])
        instructions: dict with settings from CensoSettings (specific for part, see CensoPart)
        workdir: absolute path to folder where calculations should be executed in
        """
        # try to get program from instructions
        prog = instructions.get("prog", None)
        
        if prog is None:
            raise Exception # TODO - error handling

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
        
        # TODO - set PARNODES
        if self.settings.settings_current.get_setting(str, "general", "balance"):
            # execute processes for conformers with load balancing enabled
            # divide the conformers into chunks, work through them one after the other
            chunks, procs = self.__chunkate()
                
            # execute each chunk with dynamic queue
            results = {}
            for i, chunk in enumerate(chunks):
                # TODO - this is not very nice
                self.__processor.instructions["omp"] = self.__ncores // procs[i]
                self.__omp = self.__ncores // procs[i]
                self.__nprocs = procs[i]

                # collect results by iteratively merging
                results = {**results, **self.__dqp(chunk)}
        else:
            # execute processes without load balancing, taking the user-provided settings
            results = self.__dqp(self.conformers)

        # assert that there is a result for every conformer
        # TODO - error handling
        try:
            assert all([conf.id in results.keys() for conf in self.conformers])  
        except AssertionError:
            raise KeyError("There is a mismatch between conformer ids and returned results. Cannot find at least one conformer id in results.")
        
        return results


    def __dqp(self, confs: List[GeometryData]) -> Dict[str, Any]:
        """
        D ynamic Q ueue P rocessing
        parallel execution of processes with settings defined in self.__processor.instructions
        """
        
        # execute calculations for given list of conformers
        with ProcessPoolExecutor(max_workers=self.__nprocs) as executor:
            resiter = executor.map(self.__processor.run, confs) 
        
        # returns merged result dicts
        # structure of results: 
        #   e.g. {id(conf): {"xtb_sp": {"success": ..., "energy": ...}, ...}, ...}
        return reduce(lambda x, y: {**x, **y}, resiter)
        
        
    def __chunkate(self) -> Tuple[List[Any]]:
        """
        distribute conformers until none are left
        group chunks by number of processes
        """ 
        chunks, procs = []
        
        i = 0
        pold = -1
        nconf, lconf = len(self.conformers)
        while nconf > 0:
            if nconf >= self.__ncores:
                p = self.__ncores
            elif nconf < self.__ncores:
                if self.__ncores % nconf:
                    p = nconf
                else:
                    # largest integer smaller than self.__ncores that divides nconf
                    p = max([j for j in range(1, self.__ncores) if self.__ncores % j == 0 and j <= nconf])
            
            if p != pold:
                # if number of processes is different than before add new chunk
                chunks.append(self.conformers[lconf-nconf:lconf-nconf+p])
                procs.append(p)
                i += 1
            else:
                # if number of processes didn't change, merge with the previous chunk
                chunks[i-1].extend(self.conformers[lconf-nconf:lconf-nconf+p])
            
            pold = p
            nconf -= p

        return chunks, procs