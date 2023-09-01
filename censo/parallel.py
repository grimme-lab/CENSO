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

    def __init__(self, instructions: Dict[str, Any], conformers: List[GeometryData]):
        """
        Initializes the process handler

        instructions: Single-level dictionary with settings for the calculation (specific for each part)
        conformers: List of GeometryData objects, the calculations are run using these
        """
        
        self.__instructions = instructions
        self.__paths = instructions.get("paths")

        # stores the conformers
        self.__conformers: List[GeometryData] = conformers
        
        # stores the processor
        self.__processor = None
        
        # get number of cores
        try:
            self.__ncores = os.cpu_count()
        except AttributeError:
            # TODO - error handling
            raise AttributeError("Could not determine number of cores.")
        
        # get number of processes
        self.__nprocs = self.__instructions.get("maxprocs"), 
        
        # get number of cores per process
        self.__omp = self.__instructions.get("omp"),
                

    def execute(self, jobtype: List[str], workdir: str):
        """
        creates and executes the processes
        returns the results sorted by conformer, divided into jobtypes

        jobtype: list of ordered jobtypes (e.g. [xtb_sp, xtb_gsolv])
        workdir: absolute path to folder where calculations should be executed in
        """
        # try to get program from instructions
        prog = self.__instructions.get("prog", None)
        
        if prog is None:
            raise Exception # TODO - error handling

        # initialize the processor for the respective program (depends on part)
        # and set the jobtype as well as instructions, also pass workdir to compute in
        # also pass a lookup dict to the processor so it can set the solvent for the program call correctly
        self.__processor = ProcessorFactory.create_processor(
            prog, 
            self.__paths,
            self.__instructions,
            jobtype,
            workdir
        )
        
        # TODO - set PARNODES
        if self.__instructions["balance"]:
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
            results = self.__dqp(self.__conformers)

        # assert that there is a result for every conformer
        # TODO - error handling
        try:
            assert all([conf.id in results.keys() for conf in self.__conformers])  
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
        Distributes conformers until none are left.
        Groups chunks by the number of processes.

        Returns:
            Tuple of lists: chunks - the distributed conformers,
                            procs - the number of processes for each chunk
        """ 
        chunks, procs = [], []  # Initialize empty lists to store chunks and process numbers
        
        i = 0  # Initialize a counter variable to keep track of the current chunk
        pold = -1  # Initialize a variable to store the previous number of processes used
        nconf, lconf = len(self.__conformers)  # Get the total number of conformers
        
        while nconf > 0:  # Loop until all conformers are distributed
            if nconf >= self.__ncores:  # If there are more conformers than the number of available cores
                p = self.__ncores  # Set the number of processes to the number of available cores
            elif nconf < self.__ncores:  # If there are fewer conformers than the number of available cores
                if self.__ncores % nconf:  # Check if the number of cores is not divisible by the number of conformers
                    p = nconf  # Set the number of processes to the number of conformers
                else:
                    # Find the largest integer smaller than the number of cores that divides the number of conformers
                    p = max([j for j in range(1, self.__ncores) if self.__ncores % j == 0 and j <= nconf])
            
            if p != pold:  # If the number of processes is different than before
                # Add a new chunk to the list of chunks, containing a subset of conformers
                chunks.append(self.__conformers[lconf-nconf:lconf-nconf+p])
                # Add the number of processes used for the chunk to the list of process numbers
                procs.append(p)
                i += 1  # Increment the counter variable to keep track of the current chunk
            else:  # If the number of processes didn't change
                # Merge the current chunk with the previous chunk
                chunks[i-1].extend(self.__conformers[lconf-nconf:lconf-nconf+p])
            
            pold = p  # Update the previous number of processes used
            nconf -= p  # Subtract the number of processes used from the remaining conformers

        return chunks, procs  # Return the list of chunks and the list of process numbers