import traceback
import os
from typing import Literal
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path
from concurrent.futures import CancelledError
from uuid import uuid4
from dask.distributed import Variable, Client, as_completed as dask_as_completed

from .molecules import GeometryData, MoleculeData
from .logging import setup_logger
from .params import QmProg
from .config.job_config import XTBJobConfig, SPJobConfig
from .config.parallel_config import ParallelConfig
from .config.job_config import QmResult, MetaData


logger = setup_logger(__name__)


@contextmanager
def setup_parallel(ncores: int, threads_per_worker: int):
    """
    Set up Dask parallel execution environment.

    :param ncores: Total number of cores available.
    :param threads_per_worker: Number of threads per worker.
    :return: Yields LocalCluster and Client.
    """
    from dask.distributed import LocalCluster, Client

    cluster = LocalCluster(
        n_workers=1,
        threads_per_worker=threads_per_worker,
        resources={"CPU": ncores},  # One worker receives all CPUs
    )
    client = Client(cluster)
    try:
        yield cluster, client
    finally:
        client.close()
        cluster.close()


class ParallelJob:
    """
    Represents a job for parallel execution.
    """

    def __init__(self, conf: GeometryData, charge: int, unpaired: int, omp: int):
        """
        Initialize parallel job.

        :param conf: Geometry data.
        :param charge: Charge.
        :param unpaired: Unpaired electrons.
        :param omp: Number of cores.
        """
        # conformer for the job
        self.conf: GeometryData = conf

        # number of cores to use
        self.omp: int = omp

        # stores path to an mo file which is supposed to be used as a guess
        # In case of open shell tm calculation this can be a tuple of files
        self.mo_guess: Path | str | tuple[str | Path, str | Path] | None = None

        self.from_part: str = ""

        self.charge: int = charge
        self.unpaired: int = unpaired

        # stores all flags for the jobtypes
        self.flags: dict[str, str] = {}


def set_omp(
    jobs: list[ParallelJob],
    balance: bool,
    omp: int,
    ncores: int,
    ompmin: int,
    ompmax: int,
):
    """
    Sets the number of cores to be used for every job.

    :param jobs: List of parallel jobs.
    :param balance: Whether to balance.
    :param omp: OMP value.
    :param ncores: Total cores.
    :param ompmin: Minimum OMP.
    :param ompmax: Maximum OMP.
    :return: None
    """
    if balance:
        jobs_left, tot_jobs = len(jobs), len(jobs)

        # Maximum/minimum number of parallel processes
        maxprocs = ncores // ompmin
        minprocs = max(1, ncores // ompmax)

        while jobs_left > 0:
            # Find the optimal number of parallel processes to maximize core utilization while avoiding
            # to parallelize on too many cores (efficiency drops off quickly)
            p = (
                maxprocs
                if jobs_left >= maxprocs
                else max(
                    j
                    for j in range(minprocs, maxprocs)
                    if ncores % j == 0 and j <= jobs_left
                )
            )

            # Assign determined OMP value to jobs
            while jobs_left - p >= 0:
                for job in jobs[tot_jobs - jobs_left : tot_jobs - jobs_left + p]:
                    job.omp = ncores // p
                jobs_left -= p
    else:
        for job in jobs:
            job.omp = omp


def prepare_jobs(
    conformers: list[MoleculeData],
    prog: str,
    ncores: int,
    omp: int,
    from_part: str,
    ompmin: int,
    ompmax: int,
    balance: bool = True,
    copy_mo: bool = True,
) -> list[ParallelJob]:
    """
    Prepares the jobs from the conformers data.

    :param conformers: List of conformers.
    :param prog: Program name.
    :param ncores: Total cores.
    :param omp: OMP value.
    :param from_part: From part.
    :param ompmin: Minimum OMP.
    :param ompmax: Maximum OMP.
    :param balance: Whether to balance.
    :param copy_mo: Whether to copy MO.
    :return: List of parallel jobs.
    """
    # Create ParallelJob instances from the conformers, sharing the prepinfo dict
    jobs = [
        ParallelJob(conf.geom, conf.charge, conf.unpaired, ompmin)
        for conf in conformers
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
        set_omp(jobs, balance, omp, ncores, ompmin, ompmax)

    return sorted(jobs, key=lambda job: job.omp)


def execute[T: QmResult](
    conformers: list[MoleculeData],
    task: Callable[..., tuple[T, MetaData]],
    job_config: XTBJobConfig | SPJobConfig,
    prog: QmProg | Literal["xtb"],
    from_part: str,
    parallel_config: ParallelConfig | None,
    client: Client,
    ignore_failed: bool = True,
    balance: bool = True,
    copy_mo: bool = True,
) -> dict[str, T]:
    """
    Executes the parallel tasks using Dask distributed execution.

    :param conformers: List of conformers for which jobs will be created and executed.
    :param task: Callable to be mapped onto the list of jobs created from conformers.
    :param job_config: Instructions for the execution of the task.
    :param prog: Name of the QM program.
    :param from_part: From part.
    :param parallel_config: Parallel configuration.
    :param client: dask.distributed.Client for parallel execution.
    :param ignore_failed: Whether to ignore failed jobs.
    :param balance: Whether to balance the number of cores used per job.
    :param copy_mo: Whether to store the paths to the MO files for reuse.
    :return: Job results.
    """
    # Set up parallel config if not configured
    if parallel_config is None:
        ncores = os.cpu_count()
        if ncores is None:
            raise RuntimeError(
                "ParallelConfig not provided and could not determine number of available cores."
            )
        parallel_config = ParallelConfig(ncores=ncores, omp=1)

    logger.debug("Using Dask parallel environment...")
    # Prepare jobs for parallel execution by assigning number of cores per job etc.
    jobs: list[ParallelJob] = prepare_jobs(
        conformers,
        prog,
        parallel_config.ncores,
        parallel_config.omp,
        from_part,
        parallel_config.ompmin,
        parallel_config.ompmax,
        balance=balance,
        copy_mo=copy_mo,
    )

    # Cooperative cancellation variable (name deterministic/configurable)
    cancel_var_name = getattr(parallel_config, "cancel_variable_name", None)
    if not cancel_var_name:
        cancel_var_name = f"censo_cancel_{uuid4().hex}"
    cancel_var = Variable(cancel_var_name, client=client)
    cancel_var.set(False)

    # Set cancellation variable in processor
    if hasattr(task, "__self__"):
        proc = task.__self__
        proc.cancel_var = cancel_var

    # Submit jobs directly; tasks must accept (job, job_config) and may poll Variable(cancel_var_name)
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
    results: dict[str, T] = {}
    broken = False

    logger.debug("Waiting for jobs to complete...")
    try:
        completed_iter = dask_as_completed(futures)

        # Process results as they come in
        logger.debug("Iterating over results...")
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

    # Summary behavior (same as current)
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
