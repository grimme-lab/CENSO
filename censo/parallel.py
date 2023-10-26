"""
Performs the parallel execution of the QM calls.
"""
from functools import reduce
import os
from typing import Any, Dict, List
from concurrent.futures import ProcessPoolExecutor
import atexit

from censo.procfact import ProcessorFactory
from censo.utilities import print
from censo.datastructure import GeometryData
from censo.params import OMPMIN, OMPMAX


# TODO - handle failed jobs
class ProcessHandler:
    # there might be a good argument for making this a mostly static class
    # how to handle jobs that take too long?
    # you don't want to terminate jobs in our case, since this might waste a lot of time
    # instead maybe the now free workers could be used to continue with the workload after some time has passed
    # failed jobs should be restarted (in a later chunk?) with special flags
    # TODO - what to do if a failed job is not easily remedied?
    def __init__(self, instructions: Dict[str, Any], conformers: List[GeometryData] = None):
        """
        Initializes the process handler

        instructions: Single-level dictionary with settings for the calculation (specific for each part)
        conformers (optional): List of GeometryData objects, the calculations are run using these
        """

        # all-in-one package for settings, paths, and miscellaneous information 
        self.__instructions = instructions

        # stores the conformers
        # NOTE: why not just inject core here as well? because for e.g. 'optimization' you don't actually use all of the core.conformers for every step,
        # since some might converge sooner. this is more flexible
        self.__conformers: List[GeometryData] = conformers

        # stores the processor
        self.__processor = None

        # get number of cores
        try:
            self.__ncores = os.cpu_count()
        except AttributeError:
            print(f"WARNING: Could not determine number of cores. Automatic load balancing will be disabled.")
            self.__instructions["balance"] = False

        # get number of processes
        self.__nprocs = self.__instructions["procs"]

        # get number of cores per process
        self.__omp = self.__instructions["omp"]

    def execute(self, jobtype: List[str], workdir: str):
        """
        creates and executes the processes
        returns the results sorted by conformer, divided into jobtypes

        jobtype: list of ordered jobtypes (e.g. [xtb_sp, xtb_gsolv])
        workdir: absolute path to folder where calculations should be executed in
        NOTE: automatically creates subdirs for each jobtype
        """
        # try to get program from instructions
        prog = self.__instructions.get("prog", None)

        if prog is None:
            raise RuntimeError("Could not determine program from instructions.")

        # initialize the processor for the respective program (depends on part)
        # and set the jobtype as well as instructions, also pass workdir to compute in
        # also pass a lookup dict to the processor so it can set the solvent for the program call correctly
        # TODO - maybe find a way such that the processor doesn't have to be reinitialized all the time
        self.__processor = ProcessorFactory.create_processor(
            prog,
            self.__instructions,
            jobtype,
            workdir
        )

        if self.__instructions["balance"]:
            # execute processes for conformers with load balancing enabled
            # divide the conformers into chunks, work through them one after the other
            chunks, procs = self.__chunkate()

            # execute each chunk with dynamic queue
            # TODO - add capability to use free workers after some time has elapsed (average time for each conformers per number of cores)
            # to avoid single conformers to clog up the queue for a chunk (basically move conformers between chunks?)
            results = {}
            for i, chunk in enumerate(chunks):
                # TODO - this is not very nice
                self.__processor.instructions["omp"] = self.__ncores // procs[i]
                self.__omp = self.__ncores // procs[i]
                self.__nprocs = procs[i]

                # collect results by iterative merging
                results = {**results, **self.__dqp(chunk)}
        else:
            # execute processes without automatic load balancing, taking the user-provided settings
            results = self.__dqp(self.__conformers)

        # assert that there is a result for every conformer
        try:
            assert all([conf.id in results.keys() for conf in self.__conformers])
        except AssertionError:
            raise RuntimeError(
                "There is a mismatch between conformer ids and returned results. Cannot find at least one conformer id in results.")

        return results

    @property
    def conformers(self):
        return self.__conformers

    @conformers.setter
    def conformers(self, conformers: List[GeometryData]):
        # TODO - include check?
        self.__conformers = conformers

    def __dqp(self, confs: List[GeometryData]) -> Dict[str, Any]:
        """
        D ynamic Q ueue P rocessing
        parallel execution of processes with settings defined in self.__processor.instructions
        """

        # execute calculations for given list of conformers
        with ProcessPoolExecutor(max_workers=self.__nprocs) as executor:
            # make sure that the executor exits gracefully on termination
            # TODO - is using wait=False a good option here?
            # should be fine since workers will kill programs with SIGTERM
            # wait=True leads to the workers waiting for their current task to be finished before terminating
            atexit.register(executor.shutdown, wait=False)

            # execute processes
            resiter = executor.map(self.__processor.run, confs)

            # returns merged result dicts
        # structure of results: 
        #   e.g. {id(conf): {"xtb_sp": {"success": ..., "energy": ...}, ...}, ...}
        return reduce(lambda x, y: {**x, **y}, resiter)

    def __chunkate(self) -> tuple[list[list[GeometryData]], list[int | Any]]:
        """
        Distributes conformers until none are left.
        Groups chunks by the number of processes.

        Each process shouldn't use less than OMPMIN cores.

        Returns:
            Tuple of lists: chunks - the distributed conformers,
                            procs - the number of processes for each chunk
        """

        # Initialize empty lists to store chunks and process numbers
        chunks, procs = [], []

        # Initialize a counter variable to keep track of the current chunk
        i = 0

        # Initialize a variable to store the previous number of processes used
        pold = -1

        # Get the total number of conformers
        nconf, lconf = len(self.__conformers), len(self.__conformers)

        # Calculate the maximum and minimum number of processes
        maxprocs = self.__ncores // OMPMIN
        minprocs = max(1, self.__ncores // OMPMAX)

        # Loop until all conformers are distributed
        while nconf > 0:

            # If the number of conformers is greater than or equal to the maximum number of processes
            if nconf >= maxprocs:
                p = maxprocs
            # If the number of conformers is equal to the maximum number of processes
            elif nconf == maxprocs:
                p = nconf
            # If the number of conformers is less than the minimum number of processes
            elif nconf < minprocs:
                p = minprocs
            else:
                # Find the largest number of processes that is less than or equal to the number of conformers
                # and is a divisor of the total number of cores (you basically never want to waste capacity)
                p = max([j for j in range(minprocs, maxprocs) if self.__ncores % j == 0 and j <= nconf])

            # If the number of processes is different than before
            if p != pold:
                # Add a new chunk to the list of chunks, containing a subset of conformers
                chunks.append(self.__conformers[lconf - nconf:lconf - nconf + p])
                # Add the number of processes used for the chunk to the list of process numbers
                procs.append(p)
                # Increment the counter variable to keep track of the current chunk
                i += 1
            else:
                # If the number of processes didn't change, merge the current chunk with the previous chunk
                chunks[i - 1].extend(self.__conformers[lconf - nconf:lconf - nconf + p])

            # Update the previous number of processes used
            pold = p
            # Subtract the number of processes used from the remaining conformers
            nconf -= p

        # Return the list of chunks and the list of process numbers
        return chunks, procs
