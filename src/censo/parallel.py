from inspect import ismethod
import traceback
import os
from typing import Literal, Callable
from concurrent.futures import CancelledError
from dask.distributed import (
    Variable,
    Client,
    LocalCluster,
    as_completed as dask_as_completed,
)


from .molecules import MoleculeData
from .logging import setup_logger
from .params import OMPMAX_DEFAULT, QmProg, OMPMIN_DEFAULT
from .config.job_config import XTBJobConfig, SPJobConfig
from .processing.results import Result, MetaData
from .processing.job import JobContext
from .config import GenericConfig


logger = setup_logger(__name__)


def get_cluster(maxcores: int | None = None, ompmin: int = OMPMIN_DEFAULT):
    """
    Set up Dask parallel execution environment.

    :param maxcores: Total number of cores made available for the LocalCluster.
    :type maxcores: int | None
    :param ompmin: Minimum number of cores per job.
    :type ompmin: int
    :return: Yields LocalCluster.
    :rtype: LocalCluster
    """
    nnodes = int(os.environ.get("SLURM_NNODES", "1"))
    slurm_ntasks_str = os.environ.get("SLURM_NTASKS", None)
    slurm_ntasks = None

    if slurm_ntasks_str is not None:
        slurm_ntasks = int(slurm_ntasks_str)
        if nnodes > 1:
            slurm_ntasks = slurm_ntasks // nnodes
        logger.debug(f"Found SLURM_NTASKS={slurm_ntasks} (for each node).")

    if maxcores is None:
        ncores = slurm_ntasks if slurm_ntasks is not None else os.cpu_count()
        if ncores is None:
            raise RuntimeError("Could not determine number of available cores.")
    elif slurm_ntasks is not None and maxcores <= slurm_ntasks:
        ncores = maxcores
    elif slurm_ntasks is not None and maxcores > slurm_ntasks:
        raise RuntimeError(
            "Requested more cores than slurm tasks available for this node."
        )
    else:
        ncores = maxcores

    # TODO: SSHCluster for multiple nodes
    logger.debug(f"Setting up LocalCluster with {ncores} cores.")
    cluster = LocalCluster(
        n_workers=1,
        threads_per_worker=ncores // ompmin,
        resources={"CPU": ncores},  # One worker receives all CPUs
    )
    return cluster


def set_omp(
    jobs: list[JobContext],
    balance: bool,
    omp: int,
    ncores: int,
):
    """
    Sets the number of cores to be used for every job.

    :param jobs: List of parallel jobs.
    :type jobs: list[JobContext]
    :param balance: Whether to balance.
    :type balance: bool
    :param omp: Default OMP value.
    :type omp: int
    :param ncores: Total cores.
    :type ncores: int
    :return: None
    :rtype: None
    """
    if balance:
        jobs_left, tot_jobs = len(jobs), len(jobs)

        if tot_jobs >= ncores // omp:
            # Assign minimum OMP for as many jobs as possible
            for i in range(tot_jobs - tot_jobs % omp):
                jobs[i].omp = omp

            jobs_left = tot_jobs % omp

        if jobs_left > 0:
            if jobs_left == 1:
                jobs[-1].omp = ncores
            else:
                while jobs_left > 0:
                    omps = list(range(omp, min(OMPMAX_DEFAULT, ncores) + 1))
                    jobs_in_parallel = {o: ncores // o for o in omps}

                    # Remove OMP values that would lead to idle groups
                    omps = [o for o in omps if jobs_left - jobs_in_parallel[o] >= 0]

                    divisors = [o for o in omps if ncores % o == 0 and o != ncores]
                    if len(divisors) > 0:
                        tmp_omp = min(divisors)
                    else:
                        # Not possible to evenly divide cores
                        # Pick OMP value with least idle cores that is not ncores
                        # Chances are that it is one of the smallest OMP values
                        omps = [o for o in omps if o != ncores]
                        if len(omps) > 0:
                            tmp_omp = min(omps, key=lambda x: ncores % x)
                        else:
                            tmp_omp = ncores

                    jobs[tot_jobs - jobs_left].omp = tmp_omp
                    jobs_left -= jobs_in_parallel[tmp_omp]
    else:
        for job in jobs:
            job.omp = omp


def prepare_jobs(
    conformers: list[MoleculeData],
    prog: str,
    ncores: int,
    omp: int,
    from_part: str,
    balance: bool = True,
    copy_mo: bool = True,
) -> list[JobContext]:
    """
    Prepares the jobs from the conformers data.

    :param conformers: List of conformers.
    :type conformers: list[MoleculeData]
    :param prog: Program name.
    :type prog: str
    :param ncores: Total cores per node.
    :type ncores: int
    :param omp: Default OMP value.
    :type omp: int
    :param from_part: From part.
    :type from_part: str
    :param balance: Whether to balance.
    :type balance: bool
    :param copy_mo: Whether to copy MO.
    :type copy_mo: bool
    :return: List of parallel jobs.
    :rtype: list[JobContext]
    """
    # Create JobContext instances from the conformers, sharing the prepinfo dict
    jobs = [
        JobContext(conf.geom, conf.charge, conf.unpaired, omp) for conf in conformers
    ]
    if copy_mo:
        for job in jobs:
            # Insert mo guess file path
            if prog == "orca":
                job.mo_guess = next(
                    (
                        mo_path
                        for conf in conformers
                        for mo_path in conf.mo_paths["orca"]
                        if conf.name == job.conf.name and ".gbw" in mo_path
                    ),
                    None,
                )
            elif prog == "tm":
                job.mo_guess = next(
                    (
                        mo_path
                        for conf in conformers
                        for mo_path in conf.mo_paths["tm"]
                        if conf.name == job.conf.name
                        and any(kw in mo_path for kw in ["alpha", "beta", "mos"])
                    ),
                    None,
                )
            job.from_part = from_part

    if balance:
        set_omp(jobs, balance, omp, ncores)

    return sorted(jobs, key=lambda job: job.omp)


def execute[
    T: GenericConfig, U: Result
](
    conformers: list[MoleculeData],
    task: Callable[[JobContext, T], tuple[U, MetaData]],
    job_config: XTBJobConfig | SPJobConfig,
    prog: QmProg | Literal["xtb"],
    from_part: str,
    client: Client,
    ignore_failed: bool = True,
    balance: bool = True,
    copy_mo: bool = True,
) -> dict[str, U]:
    """
    Executes the parallel tasks using Dask distributed execution.

    :param conformers: List of conformers for which jobs will be created and executed.
    :type conformers: list[MoleculeData]
    :param task: Callable to be mapped onto the list of jobs created from conformers. This should always be a processor function.
    :type task: Callable[[JobContext, T], tuple[U, MetaData]]
    :param job_config: Instructions for the execution of the task.
    :type job_config: XTBJobConfig | SPJobConfig
    :param prog: Name of the QM program.
    :type prog: QmProg | Literal["xtb"]
    :param from_part: From part.
    :type from_part: str
    :param client: dask.distributed.Client for parallel execution.
    :type client: Client
    :param ignore_failed: Whether to ignore failed jobs.
    :type ignore_failed: bool
    :param balance: Whether to balance the number of cores used per job.
    :type balance: bool
    :param copy_mo: Whether to store the paths to the MO files for reuse.
    :type copy_mo: bool
    :return: Job results.
    :rtype: dict[str, U]
    """
    # Setting the collective cancellation variable
    cancel_var = None
    if ismethod(task):
        # Only meaningful if the task is a method of a GenericProc
        # but this is difficult to annotate with the present architecture
        proc = task.__self__
        if hasattr(proc, "cancel_var"):
            proc_id = id(proc)

            # Cooperative cancellation variable
            cancel_var_name = f"censo_cancel_{proc_id}"
            cancel_var = Variable(cancel_var_name, client=client)
            cancel_var.set(False)
            setattr(proc, "cancel_var", cancel_var)

    logger.debug("Using Dask parallel environment...")

    # Get settings from client
    scheduler_info = client.scheduler_info()
    nnodes = int(scheduler_info["n_workers"])
    total_threads = scheduler_info["total_threads"]
    threads_per_worker = total_threads // nnodes
    total_cores = sum(
        [worker["resources"]["CPU"] for _, worker in scheduler_info["workers"].items()]
    )
    ncores = total_cores // nnodes
    omp = total_cores // threads_per_worker

    # Prepare jobs for parallel execution by assigning number of cores per job etc.
    jobs: list[JobContext] = prepare_jobs(
        conformers,
        prog,
        ncores,
        omp,
        from_part,
        balance=balance,
        copy_mo=copy_mo,
    )

    # Submit jobs directly
    futures = []
    for job in jobs:
        resources = {"CPU": int(job.omp)}
        fut = client.submit(
            task,
            job,
            job_config,
            resources=resources,
        )
        futures.append(fut)

    failed_confs: list[MoleculeData] = []
    results: dict[str, U] = {}
    broken = False

    logger.debug("Waiting for jobs to complete...")
    try:
        completed_iter = dask_as_completed(futures)

        # Process results as they come in
        for completed_future in completed_iter:
            conf_name = None
            try:
                result, meta = completed_future.result()
                conf_name = meta.conf_name
                conf = next(c for c in conformers if c.name == conf_name)
                if not meta.success:
                    failed_confs.append(conf)
                    logger.warning(
                        f"{meta.conf_name} job failed. Error: {meta.error}. Check output files."
                    )
                else:
                    results[meta.conf_name] = result
                    if prog in QmProg:
                        conf.mo_paths[prog].append(result.mo_path)
            except CancelledError:
                logger.debug("Future cancelled.")
                broken = True
            except Exception:
                # Fatal error: set cancel flag and best-effort cancel all remaining futures
                logger.critical(
                    "Encountered exception in task. Cancelling remaining tasks."
                )
                if cancel_var is not None:
                    cancel_var.set(True)
                try:
                    client.cancel(futures)
                except Exception:
                    logger.debug(
                        "client.cancel raised while trying to cancel futures.",
                        exc_info=True,
                    )
                tb = traceback.format_exc()
                logger.debug(tb)
                # Try to identify failure
                if conf_name is None:
                    try:
                        idx = futures.index(completed_future)
                        if 0 <= idx < len(jobs):
                            conf_name = jobs[idx].conf.name
                    except Exception:
                        pass
                logger.debug(f"Job failed for {conf_name}:\n{tb}")
                broken = True
                break
    except KeyboardInterrupt:
        logger.warning("Received keyboard interrupt. Cancelling remaining tasks...")
        if cancel_var is not None:
            cancel_var.set(True)
        try:
            client.cancel(futures)
        except Exception:
            logger.debug(
                "client.cancel raised during KeyboardInterrupt handling.", exc_info=True
            )
        raise

    if broken:
        raise RuntimeError("Exception encountered in parallel execution.")

    # Summary
    if len(failed_confs) == len(conformers):
        raise RuntimeError(
            "All jobs failed to execute. Please check your setup and output files.\n"
        )

    if len(failed_confs) > 0:
        logger.info(f"Number of failed jobs: {len(failed_confs)}.")
        logger.info("Failed conformers:")
        logger.info(", ".join(c.name for c in failed_confs))
        if not ignore_failed:
            raise RuntimeError(f"Stopping CENSO due to failed jobs in {from_part}.")
    else:
        logger.info("All jobs executed successfully.")

    return results
