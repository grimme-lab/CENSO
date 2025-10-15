import traceback
import os
import subprocess
from typing import Literal
from collections.abc import Callable
from pathlib import Path
from concurrent.futures import CancelledError, ThreadPoolExecutor, as_completed
from uuid import uuid4
from dask.distributed import (
    Variable,
    Client,
    as_completed as dask_as_completed,
    SSHCluster,
    LocalCluster,
)

from .molecules import GeometryData, MoleculeData
from .logging import setup_logger
from .params import QmProg
from .config.job_config import XTBJobConfig, SPJobConfig
from .config.parallel_config import ParallelConfig
from .config.job_config import QmResult, MetaData


logger = setup_logger(__name__)


def _parse_slurm_nodelist(nodelist: str) -> list[str]:
    """
    Parse SLURM_NODELIST environment variable to extract node hostnames.

    Handles formats like:
    - "node001"
    - "node[001-003]"
    - "node001,node002,node003"
    - "node[001-002,004]"

    :param nodelist: The SLURM_NODELIST string
    :return: List of node hostnames
    """
    nodes = []

    # Split by commas first, but be careful with commas inside brackets
    parts = []
    current_part = ""
    bracket_depth = 0

    for char in nodelist:
        if char == "[":
            bracket_depth += 1
            current_part += char
        elif char == "]":
            bracket_depth -= 1
            current_part += char
        elif char == "," and bracket_depth == 0:
            parts.append(current_part)
            current_part = ""
        else:
            current_part += char

    if current_part:
        parts.append(current_part)

    for part in parts:
        part = part.strip()
        if "[" in part and "]" in part:
            # Handle bracketed ranges like "node[001-003]"
            prefix = part.split("[")[0]
            range_part = part.split("[")[1].split("]")[0]

            # Split ranges by comma
            ranges = range_part.split(",")
            for r in ranges:
                if "-" in r:
                    # Handle range like "001-003"
                    start, end = r.split("-")
                    start_num = int(start)
                    end_num = int(end)
                    width = len(start)
                    for i in range(start_num, end_num + 1):
                        nodes.append(f"{prefix}{i:0{width}d}")
                else:
                    # Single number
                    width = len(r)
                    nodes.append(f"{prefix}{int(r):0{width}d}")
        else:
            # Simple hostname
            nodes.append(part)

    return nodes


def _ssh_check_one_host(host: str, timeout: float, connect_timeout: int) -> bool:
    """
    Return True if an SSH call to `host` succeeds within `timeout`.
    """
    cmd = [
        "ssh",
        "-o", "BatchMode=yes",
        "-o", f"ConnectTimeout={connect_timeout}",
        "-o", "StrictHostKeyChecking=no",
        host,
        "true",
    ]
    try:
        completed = subprocess.run(cmd, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return completed.returncode == 0
    except subprocess.TimeoutExpired:
        logger.debug("SSH timeout for host %s", host)
        return False
    except Exception:
        logger.debug("SSH error for host %s", host, exc_info=True)
        return False


def _test_ssh_connectivity(nodes: list[str], timeout: float = 5.0, max_workers: int = 16) -> bool:
    """
    Concurrently verify SSH connectivity to all nodes using ThreadPoolExecutor.

    Returns True only if every node check succeeds. On the first failing check this
    function returns False (get_client will fall back to LocalCluster).
    """
    if not nodes:
        return False

    connect_timeout = max(1, int(timeout))
    workers = min(max_workers, len(nodes))

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_ssh_check_one_host, h, timeout, connect_timeout): h for h in nodes}
        for fut in as_completed(futures):
            try:
                ok = fut.result()
            except Exception:
                ok = False
            if not ok:
                return False
    return True


def get_client(ncores: int, threads_per_worker: int):
    """
    Set up Dask parallel execution environment.

    :param ncores: Total number of cores available.
    :param threads_per_worker: Number of threads per worker.
    :return: Yields LocalCluster or SSHCluster and Client.
    """
    # Check if we're in a Slurm job
    slurm_job_id = os.environ.get("SLURM_JOBID")
    slurm_nodelist = os.environ.get("SLURM_NODELIST")

    if slurm_job_id and slurm_nodelist:
        # We're in a Slurm job - use SSHCluster for multi-node
        logger.info(f"Detected Slurm job {slurm_job_id}")

        try:
            nodes = _parse_slurm_nodelist(slurm_nodelist)
            logger.debug(f"Parsed nodes from SLURM_NODELIST: {nodes}")

            if len(nodes) > 1:
                # Test SSH connectivity before creating SSHCluster
                logger.debug("Testing SSH connectivity to nodes...")
                if not _test_ssh_connectivity(nodes):
                    logger.warning(
                        "SSH connectivity test failed. Falling back to LocalCluster."
                    )
                    cluster = LocalCluster(
                        n_workers=1,
                        threads_per_worker=threads_per_worker,
                        resources={"CPU": ncores},
                    )
                else:
                    # Multi-node job - use SSHCluster
                    logger.info(f"Creating SSHCluster with {len(nodes)} nodes")
                    cluster = SSHCluster(
                        [nodes[0]] + nodes,  # double count the first node as scheduler
                        connect_options={"known_hosts": None},  # Skip host key checking
                        worker_options={
                            "nthreads": threads_per_worker,
                            "resources": {"CPU": ncores},
                        },
                        scheduler_options={"port": 0},  # Use any available port
                    )
                    logger.info(
                        f"Created SSHCluster with {len(nodes)} nodes. Scheduler: {nodes[0]}."
                    )
            else:
                # Single node in Slurm - fall back to LocalCluster
                logger.debug("Single node in Slurm job, using LocalCluster")
                cluster = LocalCluster(
                    n_workers=1,
                    threads_per_worker=threads_per_worker,
                    resources={"CPU": ncores},
                )
        except Exception as e:
            logger.warning(
                f"Failed to create SSHCluster: {e}. Falling back to LocalCluster."
            )
            cluster = LocalCluster(
                n_workers=1,
                threads_per_worker=threads_per_worker,
                resources={"CPU": ncores},
            )
    else:
        # Not in Slurm - use LocalCluster
        cluster = LocalCluster(
            n_workers=1,
            threads_per_worker=threads_per_worker,
            resources={"CPU": ncores},
        )

    client = cluster.get_client()
    return client, cluster


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
