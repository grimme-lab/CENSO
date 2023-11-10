"""
Performs the parallel execution of the QM calls.
"""
from functools import reduce
import os
from typing import Any, Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor
import atexit

from censo.procfact import ProcessorFactory
from censo.qm_processor import QmProc
from censo.utilities import print
from censo.datastructure import GeometryData
from censo.params import OMPMIN, OMPMAX


# extend ProcessPoolExecutor?
# fixed number of workers
# keep track of resources
# dynamic load balancing
# TODO - handle failed jobs
class ProcessHandler:
    # there might be a good argument for making this a mostly static class
    # how to handle jobs that take too long?
    # you don't want to terminate jobs in our case, since this might waste a lot of time
    # instead maybe the now free workers could be used to continue with the workload after some time has passed
    # failed jobs should be restarted (in a later chunk?) with special flags
    # TODO - what to do if a failed job is not easily remedied?
    # get number of cores
    try:
        __ncores: int = os.cpu_count()
    except AttributeError:
        __ncores = None

    __processor: QmProc
    __nprocs: int
    __omp: int

    def __init__(self, instructions: Dict[str, Any]):
        """
        Initializes the process handler

        instructions: Single-level dictionary with settings for the calculation (specific for each part)
        conformers (optional): List of GeometryData objects, the calculations are run using these
        """

        # all-in-one package for settings, paths, and miscellaneous information 
        self.__instructions = instructions

        # disable automatic load balancing if number of cores could not be determined
        if self.__ncores is None:
            print(f"WARNING: Could not determine number of cores. Automatic load balancing will be disabled.")
            self.__instructions["balance"] = False

        # get number of processes
        self.__nprocs = self.__instructions["procs"]

        # get number of cores per process
        self.__omp = self.__instructions["omp"]

    def execute(self, jobtype: List[str], workdir: str, conformers: List[GeometryData]) -> Dict[int, Any]:
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

        # TODO - add capability to rerun failed jobs
        if self.__instructions["balance"]:
            # execute processes for conformers with load balancing enabled
            # divide the conformers into chunks, work through them one after the other
            chunks, procs = self.__chunkate(conformers)

            # execute each chunk with dynamic queue
            # TODO - add capability to use free workers after some time has elapsed (average time for each conformers per number of cores)
            # to avoid single conformers to clog up the queue for a chunk (basically move conformers between chunks?)
            results = {}
            metadata = {}
            for i, chunk in enumerate(chunks):
                # TODO - this is not very nice
                self.__processor.instructions["omp"] = self.__ncores // procs[i]
                self.__omp = self.__ncores // procs[i]
                self.__nprocs = procs[i]

                # collect results by iterative merging
                tmpres, tmpmeta = self.__dqp(chunk)
                results = {**results, **tmpres}
                metadata = {**metadata, **tmpmeta}
        else:
            # execute processes without automatic load balancing, taking the user-provided settings
            results, metadata = self.__dqp(conformers)

        # assert that there is a result for every conformer
        try:
            assert all([conf.id in results.keys() for conf in conformers])
        except AssertionError:
            raise RuntimeError(
                "There is a mismatch between conformer ids and returned results. Cannot find at least one conformer id in results.")

        # if 'copy_mo' is enabled, try to get the mo_path from the results and store it in the respective GeometryData object
        if self.__instructions["copy_mo"]:
            for conf in conformers:
                if metadata[conf.id]["sp"]["mo_path"] is not None:
                    conf.mo_path = metadata[conf.id]["sp"]["mo_path"]

        # TODO - create a new list of failed jobs that should be restarted with special flags
        if self.__instructions["retry_failed"]:
            failed_jobs = {}
            errors = {}
            for conf in conformers:
                failed_jobs[conf.id] = [job for job in metadata[conf.id].keys() if metadata[conf.id][job]["success"] is False]
                errors[conf.id] = [metadata[conf.id][job]["error"] for job in metadata[conf.id].keys() if metadata[conf.id][job]["success"] is False]

            retry = []
            for conf, jobs in failed_jobs.items():
                for job in jobs:
                    if (job == "sp" or job == "gsolv") and errors[conf][job] == "SCF not converged":
                        retry.append(conf)

        return results

    @staticmethod
    def __dqp(confs: List[GeometryData]) -> Tuple[Dict[int, Any], Dict[int, Any]]:
        """
        D ynamic Q ueue P rocessing
        parallel execution of processes with settings defined in __processor.instructions
        """

        # execute calculations for given list of conformers
        with ProcessPoolExecutor(max_workers=ProcessHandler.__nprocs) as executor:
            # make sure that the executor exits gracefully on termination
            # TODO - is using wait=False a good option here?
            # should be fine since workers will kill programs with SIGTERM
            # wait=True leads to the workers waiting for their current task to be finished before terminating
            atexit.register(executor.shutdown, wait=False)

            # execute tasks
            resiter = executor.map(ProcessHandler.__processor.run, confs)

            results, metadata = zip(*[res for res in resiter])

            # returns merged result dicts
        # structure of results: 
        #   e.g. {id(conf): {"xtb_sp": {"success": ..., "energy": ...}, ...}, ...}
        return reduce(lambda x, y: {**x, **y}, results), reduce(lambda x, y: {**x, **y}, metadata)

    @staticmethod
    def __chunkate(conformers: list[GeometryData]) -> tuple[list[list[GeometryData]], list[int | Any]]:
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
        nconf, lconf = len(conformers), len(conformers)

        # Calculate the maximum and minimum number of processes
        maxprocs = ProcessHandler.__ncores // OMPMIN
        minprocs = max(1, ProcessHandler.__ncores // OMPMAX)

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
                p = max([j for j in range(minprocs, maxprocs) if ProcessHandler.__ncores % j == 0 and j <= nconf])

            # If the number of processes is different than before
            if p != pold:
                # Add a new chunk to the list of chunks, containing a subset of conformers
                chunks.append(conformers[lconf - nconf:lconf - nconf + p])
                # Add the number of processes used for the chunk to the list of process numbers
                procs.append(p)
                # Increment the counter variable to keep track of the current chunk
                i += 1
            else:
                # If the number of processes didn't change, merge the current chunk with the previous chunk
                chunks[i - 1].extend(conformers[lconf - nconf:lconf - nconf + p])

            # Update the previous number of processes used
            pold = p
            # Subtract the number of processes used from the remaining conformers
            nconf -= p

        # Return the list of chunks and the list of process numbers
        return chunks, procs
