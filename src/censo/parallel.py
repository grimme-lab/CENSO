"""
Performs the parallel execution of the QM calls.
"""

import multiprocessing
import os
import signal
from concurrent.futures import ProcessPoolExecutor, as_completed

from .datastructure import MoleculeData, ParallelJob
from .logging import setup_logger
from .params import OMPMAX, OMPMIN, ENVIRON
from .procfact import ProcessorFactory
from .qm_processor import QmProc
from .tm_processor import TmProc

logger = setup_logger(__name__)

# get number of cores
ncores = os.cpu_count()


def execute(
    conformers: list[MoleculeData],
    workdir: str,
    prog: str,
    prepinfo: dict[str, dict],
    jobtype: list[str],
    copy_mo: bool = False,
    retry_failed: bool = False,
    balance: bool = True,
    maxcores: int = OMPMIN,
    omp: int = OMPMIN,
    update: bool = True,
) -> tuple[bool, dict, list]:
    """
    Manages parallel execution of external program calls. Sets cores used per job, checks requirements,
    can copy MO-files, and retry failed jobs.

    Args:
        conformers (list[MoleculeData]): List of conformers for which jobs will be created and executed.
        workdir (str): Working directory.
        prog (str): Name of the program to be used.
        copy_mo (bool, optional): Whether to copy the MO-files from the previous calculation.
        retry_failed (bool, optional): Whether to retry failed jobs.
        balance (bool, optional): Whether to balance the number of cores used per job.
        maxcores (int, optional): Maximum number of cores to be used.
        omp (int, optional): Number of cores to be used per job.
        update (bool, optional): Wether to update the results dict for each conformer.

    Returns:
        tuple[bool, dict, list]: Whether the execution was successful, a dictionary containing the results, and a list
            of conformers that could not be recovered.
    """

    # Check first if there are any conformers at all
    try:
        assert len(conformers) > 0
    except AssertionError as e:
        raise AssertionError("No jobs to compute!") from e

    # Create jobs from conformers data
    jobs = prepare_jobs(conformers, prepinfo, jobtype)

    # initialize the processor for the respective program
    processor = ProcessorFactory.create_processor(
        prog,
        workdir,
    )

    # Adjust value for the number of available cores
    global ncores
    ncores = min(ncores, maxcores)

    # processor.check_requirements(jobs)

    # Set processor to copy the MO-files
    processor.copy_mo = copy_mo

    # check for the most recent mo files for each conformer
    # TODO - how would this work when multiple different programs are supported?
    for job in jobs:
        try:
            job.mo_guess = next(
                c for c in conformers if c.name == job.conf.name
            ).mo_paths[-1]
        except IndexError:
            pass

    # set cores per process for each job
    # NOTE: since parallelization in tm is controlled using environment variables we cannot use automatic load balancing
    if balance and not isinstance(processor, TmProc):
        set_omp_chunking(jobs)
    else:
        if omp < OMPMIN:
            logger.warning(
                f"User OMP setting is below the minimum value of {OMPMIN}. Using {OMPMIN} instead."
            )
            omp = OMPMIN
        elif omp > ncores:
            logger.warning(
                f"Value of {omp} for OMP is larger than the number of available cores {ncores}. Using OMP = {ncores}."
            )
            omp = ncores

        if isinstance(processor, TmProc):
            if balance:
                logger.warning(
                    f"Automatic load balancing is not supported for TURBOMOLE. Using OMP = {omp} for all jobs instead."
                )

            # Configure environment variables
            ENVIRON["PARA_ARCH"] = "SMP"
            ENVIRON["PARNODES"] = str(omp)

        for job in jobs:
            job.omp = omp

    # execute the jobs
    jobs = dqp(jobs, processor)

    # Try to get the mo_path from metadata and store it in the respective conformer object
    mo_paths = {job.conf.name: job.meta["mo_path"] for job in jobs}
    for conf in conformers:
        if mo_paths[conf.name] is not None:
            conf.mo_paths.append(mo_paths[conf.name])

    if retry_failed:
        retried, failed_confs = retry_failed_jobs(jobs, processor, balance)

        # Again, try to get the mo_path from metadata and store it in the respective conformer object
        mo_paths = {
            job.conf.name: job.meta["mo_path"] for job in [jobs[i] for i in retried]
        }
        for conf in conformers:
            if mo_paths.get(conf.name, None) is not None:
                conf.mo_paths.append(mo_paths[conf.name])

    # special warning if all jobs failed
    if len(jobs) == len(failed_confs):
        logger.warning("All jobs failed and could not be recovered!")

    # update results for all conformers (if enabled)
    if update:
        conformers = sorted(conformers, key=lambda x: x.name)
        jobs = sorted(jobs, key=lambda x: x.conf.name)
        for conf, job in zip(conformers, jobs):
            conf.results.setdefault(job.prepinfo["partname"], {}).update(job.results)

    return (
        len(conformers) != len(failed_confs),
        {job.conf.name: job.results for job in jobs},
        failed_confs,
    )


def prepare_jobs(
    conformers: list[MoleculeData], prepinfo: dict[str, dict], jobtype: list[str]
) -> list[ParallelJob]:
    # create jobs from conformers
    jobs = [ParallelJob(conf.geom, jobtype) for conf in conformers]

    # put settings into jobs
    for job in jobs:
        job.prepinfo.update(prepinfo)

    return jobs


def reduce_cores(
    free_cores: multiprocessing.Value, omp: int, enough_cores: multiprocessing.Condition
):
    # acquire lock on the condition and wait until enough cores are available
    with enough_cores:
        enough_cores.wait_for(lambda: free_cores.value >= omp)
        free_cores.value -= omp
        logger.debug(
            f"Free cores decreased {free_cores.value + omp} -> {free_cores.value}."
        )


def increase_cores(
    free_cores: multiprocessing.Value, omp: int, enough_cores: multiprocessing.Condition
):
    # acquire lock on the condition and increase the number of cores, notifying one waiting process
    with enough_cores:
        free_cores.value += omp
        logger.debug(
            f"Free cores increased {free_cores.value - omp} -> {free_cores.value}."
        )
        enough_cores.notify()


def handle_sigterm(signum, frame, executor):
    logger.critical("Received SIGTERM. Terminating.")
    executor.shutdown(wait=False)


def dqp(jobs: list[ParallelJob], processor: QmProc) -> list[ParallelJob]:
    """
    D ynamic Q ueue P rocessing
    """

    global ncores

    with multiprocessing.Manager() as manager:
        # execute calculations for given list of conformers
        with ProcessPoolExecutor(
            max_workers=ncores // min(job.omp for job in jobs)
        ) as executor:
            # make sure that the executor exits gracefully on termination
            # TODO - is using wait=False a good option here?
            # should be fine since workers will kill programs with SIGTERM
            # wait=True leads to the workers waiting for their current task to be finished before terminating
            # Register the signal handler
            signal.signal(
                signal.SIGTERM,
                lambda signum, frame: handle_sigterm(signum, frame, executor),
            )

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

                try:
                    # submit the job
                    tasks.append(executor.submit(processor.run, jobs[i]))
                    # NOTE: explanation of the lambda: the first argument passed to the done_callback is always the future
                    # itself, it is not assigned (_), the second parameter is the number of openmp threads of the job (i.e.
                    # job.omp) if this is not specified like this (omp=jobs[i].omp) the done_callback will instead use the
                    # omp of the current item in the for-iterator (e.g. the submitted job has omp=4, but the current jobs[i]
                    # has omp=7, so the callback would use 7 instead of 4)
                    tasks[-1].add_done_callback(
                        lambda _, omp=jobs[i].omp: increase_cores(
                            free_cores, omp, enough_cores
                        )
                    )
                except RuntimeError:
                    # Makes this exit gracefully in case that the main process is killed
                    return None

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
    # Calculate the minimum number of processes
    minprocs = max(1, ncores // OMPMAX)

    # Loop until all jobs are distributed
    while jobs_left > 0:
        if jobs_left >= maxprocs:
            p = maxprocs  # Set the number of processes to the maximum if there are enough jobs left
        elif minprocs <= jobs_left < maxprocs:
            # Find the largest number of processes that evenly divides the remaining jobs
            p = max(
                [
                    j
                    for j in range(minprocs, maxprocs)
                    if ncores % j == 0 and j <= jobs_left
                ]
            )
        else:
            # There are not enough jobs left for at least minprocs processes
            for job in jobs[tot_jobs - jobs_left : tot_jobs]:
                job.omp = (
                    ncores // minprocs
                )  # Set the number of cores for each job to the maximum value
            jobs_left -= jobs_left
            continue

        # Set the number of cores for each job for as many jobs as possible before moving onto the next omp value
        while jobs_left - p >= 0:
            for job in jobs[tot_jobs - jobs_left : tot_jobs - jobs_left + p]:
                job.omp = ncores // p  # Set the number of cores for each job
            jobs_left -= p  # Decrement the number of remaining jobs


def retry_failed_jobs(
    jobs: list[ParallelJob], processor: QmProc, balance: bool
) -> tuple[list[int], list[str]]:
    """
    Tries to recover failed jobs.

    Args:
        jobs (list[ParallelJob]): List of jobs.
        processor (QmProc): Processor object.

    Returns:
        tuple[list[int], list[str]]: List of indices of jobs that should be retried, list of names of conformers
            that could not be recovered.
    """
    # determine failed jobs
    logger.debug("Checking for failed jobs...")
    failed_jobs = [
        i
        for i, job in enumerate(jobs)
        if any(not job.meta[jt]["success"] for jt in job.jobtype)
    ]

    if len(failed_jobs) != 0:
        # create a new list of failed jobs that should be restarted with special flags
        # contains jobs that should be retried (depends on wether the error can be handled or not)
        retry = []

        # determine flags for jobs based on error messages
        for failed_job in failed_jobs:
            handled_errors = ["scf_not_converged", "Previous calculation failed"]

            # list of jobtypes that should be removed from the jobtype list
            jtremove = []
            for jt in jobs[failed_job].jobtype:
                if not jobs[failed_job].meta[jt]["success"]:
                    if jobs[failed_job].meta[jt]["error"] in handled_errors:
                        retry.append(failed_job)
                        jobs[failed_job].flags[jt] = jobs[failed_job].meta[jt]["error"]
                # store all successful jobtypes to be removed later
                elif jobs[failed_job].meta[jt]["success"]:
                    jtremove.append(jt)

            # remove all successful jobs from jobtype to avoid re-execution
            for jt in jtremove:
                jobs[failed_job].jobtype.remove(jt)

        # execute jobs that should be retried
        logger.info(
            f"Number of failed jobs: {len(failed_jobs)}. Restarting {len(retry)} jobs."
        )

        if len(retry) > 0:
            # Rebalancing necessary
            if balance:
                set_omp_chunking([jobs[i] for i in retry])

            for i, job in zip(
                [i for i in retry], dqp([jobs[i] for i in retry], processor)
            ):
                jobs[i] = job

        # any jobs that still failed will lead to the conformer being marked as unrecoverable
        failed_confs = []
        for job in jobs:
            if not all(job.meta[jt]["success"] for jt in job.jobtype):
                logger.warning(
                    f"{job.conf.name} job recovery failed. Error: {job.meta[jt]['error']}. Check output files."
                )
                failed_confs.append(job.conf.name)
            else:
                logger.info(f"Successfully retried job for {job.conf.name}.")
    else:
        retry = []
        failed_confs = []
        logger.info("All jobs executed successfully.")

    return retry, failed_confs
