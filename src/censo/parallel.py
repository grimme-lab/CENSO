"""
Performs the parallel execution of the QM calls.
"""
import multiprocessing
import os
import signal
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List

from .datastructure import MoleculeData, ParallelJob
from .params import OMPMIN, OMPMAX
from .procfact import ProcessorFactory
from .qm_processor import QmProc
from .utilities import setup_logger

logger = setup_logger(__name__)

# get number of cores
ncores = os.cpu_count()


def execute(conformers: List[MoleculeData], instructions: Dict[str, Any], workdir: str) -> Dict[int, Any]:
    # try to get program from instructions
    prog = instructions.get("prog", None)

    if prog is None:
        raise RuntimeError("Could not determine program from instructions.")

    # initialize the processor for the respective program
    processor = ProcessorFactory.create_processor(
        prog,
        instructions,
        workdir,
    )

    balance = instructions["balance"]

    # adjust value for the number of available cores
    global ncores
    ncores = min(ncores, instructions["maxcores"])

    # create jobs, by default with the omp setting from instructions
    jobs = [ParallelJob(conf.geom, instructions["jobtype"], instructions["omp"]) for conf in conformers]

    if instructions["copy_mo"]:
        # check for the most recent mo files for each conformer
        # TODO - how would this work when multiple different programs are supported?
        for job in jobs:
            try:
                job.mo_guess = next(c for c in conformers if c.geom.id == job.conf.id).mo_paths[-1]
            except IndexError:
                pass

    # set cores per process for each job 
    if balance:
        set_omp_chunking(jobs)

    # execute the jobs
    # TODO - let this exit gracefully when process is killed
    # RuntimeError: cannot schedule new futures after shutdown
    jobs = dqp(jobs, processor)

    # if 'copy_mo' is enabled, try to get the mo_path from metadata and store it in the respective conformer object
    if instructions["copy_mo"]:
        mo_paths = {job.conf.id: job.meta["mo_path"] for job in jobs}
        for conf in conformers:
            if mo_paths[conf.geom.id] is not None:
                conf.mo_paths.append(mo_paths[conf.geom.id])

    # create a new list of failed jobs that should be restarted with special flags
    if instructions["retry_failed"]:
        # determine failed jobs
        logger.debug("Retrying failed jobs...")
        failed_jobs = [i for i, job in enumerate(jobs) if any(not job.meta[jt]["success"] for jt in job.jobtype)]

        if len(failed_jobs) != 0:

            # contains jobs that should be retried (depends if the error can be handled or not)
            retry = []

            # determine flags for jobs based on error messages
            for failed_job in failed_jobs:
                for jt in jobs[failed_job].jobtype:
                    # for now only sp and gsolv calculations are caught
                    if not jobs[failed_job].meta[jt]["success"] and jt in ["sp", "gsolv"]:
                        if jobs[failed_job].meta[jt]["error"] == "SCF not converged":
                            retry.append(failed_job)
                            jobs[failed_job].flags[jt] = "scf_not_converged"
                    # remove all successful jobs from jobtype to avoid re-execution
                    elif jobs[failed_job].meta[jt]["success"]:
                        jobs[failed_job].jobtype.remove(jt)

            # execute jobs that should be retried
            logger.info(f"Failed jobs: {len(failed_jobs)}. Restarting {len(retry)} jobs.")

            if len(retry) > 0:
                set_omp_chunking([jobs[i] for i in retry])
                for i, job in zip([i for i in retry], dqp([jobs[i] for i in retry], processor)):
                    jobs[i] = job

            # any jobs that still failed will lead to the conformer to be removed from the list (TODO)
            for job in jobs:
                if not all(job.meta[jt]["success"] for jt in job.jobtype):
                    logger.warning(f"Removed {job.conf.name} from conformer list due to failed jobs.")
                    conformers.remove(next(c for c in conformers if c.geom.id == job.conf.id))

            # again, try to get the mo_path from metadata and store it in the respective conformer object
            if instructions["copy_mo"]:
                mo_paths = {job.conf.id: job.meta["mo_path"] for job in [jobs[i] for i in retry]}
                for conf in conformers:
                    if mo_paths.get(conf.geom.id, None) is not None:
                        conf.mo_paths.append(mo_paths[conf.geom.id])
        else:
            logger.info("All jobs executed successfully.")

    # collect all results from the job objects
    return {job.conf.id: job.results for job in jobs}


def reduce_cores(free_cores: multiprocessing.Value, omp: int, enough_cores: multiprocessing.Condition):
    # acquire lock on the condition and wait until enough cores are available
    with enough_cores:
        enough_cores.wait_for(lambda: free_cores.value >= omp)
        free_cores.value -= omp
        logger.debug(f"Free cores decreased {free_cores.value + omp} -> {free_cores.value}.")


def increase_cores(free_cores: multiprocessing.Value, omp: int, enough_cores: multiprocessing.Condition):
    # acquire lock on the condition and increase the number of cores, notifying one waiting process
    with enough_cores:
        free_cores.value += omp
        logger.debug(f"Free cores increased {free_cores.value - omp} -> {free_cores.value}.")
        enough_cores.notify()


def handle_sigterm(signum, frame, executor):
    logger.critical("Received SIGTERM. Terminating.")
    executor.shutdown(wait=False)


def dqp(jobs: List[ParallelJob], processor: QmProc) -> list[ParallelJob]:
    """
    D ynamic Q ueue P rocessing
    """

    global ncores

    with multiprocessing.Manager() as manager:
        # execute calculations for given list of conformers
        with ProcessPoolExecutor(max_workers=ncores // min(job.omp for job in jobs)) as executor:
            # make sure that the executor exits gracefully on termination
            # TODO - is using wait=False a good option here?
            # should be fine since workers will kill programs with SIGTERM
            # wait=True leads to the workers waiting for their current task to be finished before terminating
            # Register the signal handler
            signal.signal(signal.SIGTERM, lambda signum, frame: handle_sigterm(signum, frame, executor))

            # define shared variables that can be safely asynchronously accessed
            free_cores = manager.Value(int, ncores)
            enough_cores = manager.Condition()

            # sort the jobs by the number of cores used
            # (the first item will be the one with the lowest number of cores)
            jobs.sort(key=lambda x: x.omp)

            tasks = []
            for i in range(len(jobs)):
                # TODO - something to readjust omp based on expected time to finish and the timings of other jobs
                # try to reduce the number of cores by job.omp, if there are not enough cores available we wait
                reduce_cores(free_cores, jobs[i].omp, enough_cores)

                # submit the job
                tasks.append(executor.submit(processor.run, jobs[i]))
                # NOTE: explanation of the lambda: the first argument passed to the done_callback is always the future
                # itself, it is not assigned (_), the second parameter is the number of openmp threads of the job (i.e.
                # job.omp) if this is not specified like this (omp=jobs[i].omp) the done_callback will instead use the
                # omp of the current item in the for-iterator (e.g. the submitted job has omp=4, but the current jobs[i]
                # has omp=7, so the callback would use 7 instead of 4)
                tasks[-1].add_done_callback(lambda _, omp=jobs[i].omp: increase_cores(free_cores, omp, enough_cores))

            # wait for all jobs to finish and collect results
            results = [task.result() for task in as_completed(tasks)]

    return results


def set_omp_chunking(jobs: list[ParallelJob]) -> None:
    """
    Determines and sets the number of cores that are supposed to be used for every job.
    This method is efficient if it can be assumed that the jobs take roughly the same amount of time each.
    Each job shouldn't use less than OMPMIN cores.
    """
    global ncores  # Access the global variable ncores

    # Get the total number of jobs
    jobs_left, tot_jobs = len(jobs), len(jobs)

    # Calculate the maximum and minimum number of processes (number of jobs that can be executed simultaneously)
    maxprocs = ncores // OMPMIN  # Calculate the maximum number of processes
    minprocs = max(1, ncores // OMPMAX)  # Calculate the minimum number of processes

    # Loop until all jobs are distributed
    while jobs_left > 0:
        if jobs_left >= maxprocs:
            p = maxprocs  # Set the number of processes to the maximum if there are enough jobs left
        elif jobs_left < minprocs:
            p = minprocs  # Set the number of processes to the minimum if there are less jobs left than minprocs
        else:
            # Find the largest number of processes that evenly divides the remaining jobs
            p = max([j for j in range(minprocs, maxprocs) if ncores % j == 0 and j <= jobs_left])

        # Set the number of cores for each job for as many jobs as possible before moving onto the next omp value
        while jobs_left - p >= 0:
            for job in jobs[tot_jobs - jobs_left:tot_jobs - jobs_left + p]:
                job.omp = ncores // p  # Set the number of cores for each job
            jobs_left -= p  # Decrement the number of remaining jobs
