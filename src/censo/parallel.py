import traceback
from dataclasses import dataclass, field
from typing import Callable, Literal
from contextlib import contextmanager
from pathlib import Path
from multiprocessing.managers import SyncManager
from multiprocessing import Manager
from concurrent.futures import Future, ProcessPoolExecutor, as_completed

from .molecules import Atom, GeometryData, MoleculeData
from .logging import setup_logger
from .params import QmProg, OMPMIN, OMPMAX, ENVIRON
from .config.job_config import XTBJobConfig, SPJobConfig

logger = setup_logger(__name__)


@contextmanager
def setup_managers(max_workers: int, ncores: int):
    executor: ProcessPoolExecutor = ProcessPoolExecutor(max_workers=max_workers)
    manager: SyncManager = Manager()
    resource_manager: ResourceMonitor = ResourceMonitor(manager, ncores)
    try:
        yield executor, manager, resource_manager
    finally:
        executor.shutdown(False, cancel_futures=True)
        manager.shutdown()


class ResourceMonitor:
    def __init__(self, manager: SyncManager, ncores: int):
        self.__free_cores = manager.Value(int, ncores)
        self.__enough_cores = manager.Condition()

    @contextmanager
    def occupy_cores(self, ncores: int):
        try:
            with self.__enough_cores:
                self.__enough_cores.wait_for(lambda: self.__free_cores.value >= ncores)
                self.__free_cores.value -= ncores
            yield
        finally:
            with self.__enough_cores:
                self.__free_cores.value += ncores
                self.__enough_cores.notify()  # TODO: try this with notify_all instead


@dataclass
class QmResult:
    """Base class for results returned by QM calculations."""

    mo_path: str | tuple[str, str] = ""


@dataclass
class SPResult(QmResult):
    """Results class for single-point calculations."""

    energy: float = 0.0


@dataclass
class GsolvResult(QmResult):
    """Results class for Gsolv calculations."""

    gsolv: float = 0.0
    energy_gas: float = 0.0
    energy_solv: float = 0.0


@dataclass
class RRHOResult(QmResult):
    """Results class for mRRHO calculations."""

    energy: float = 0.0
    rmsd: float = 0.0
    gibbs: dict[float, float] = field(default_factory=dict)
    enthalpy: dict[float, float] = field(default_factory=dict)
    entropy: dict[float, float] = field(default_factory=dict)
    symmetry: str = "c1"
    symnum: int = 1
    linear: bool = False


@dataclass
class OptResult(QmResult):
    """Results class for geometry optimizations."""

    # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
    ecyc: list[float] = field(default_factory=list)
    # 'gncyc' contains the gradient norms for all cycles
    gncyc: list[float] = field(default_factory=list)
    # 'energy' contains the final energy of the optimization (converged or unconverged)
    energy: float = 0.0
    # 'geom' stores the optimized geometry in GeometryData.xyz format
    geom: list[Atom] = field(default_factory=list)
    grad_norm: float = 0.0
    cycles: int = 0
    converged: bool = False


@dataclass
class NMRResult(QmResult):
    """Results class for NMR calculations."""

    shieldings: list[tuple[int, float]] = field(default_factory=list)
    couplings: list[tuple[tuple[int, int], float]] = field(default_factory=list)


@dataclass
class UVVisResult(QmResult):
    """Results class for UVVis calculations."""

    excitations: list[dict[str, float]] = field(default_factory=list)


@dataclass
class MetaData:
    conf_name: str
    success: bool = False
    error: str = ""


class ParallelJob:

    def __init__(self, conf: GeometryData, charge: int, unpaired: int):
        # conformer for the job
        self.conf: GeometryData = conf

        # number of cores to use
        self.omp: int = OMPMIN

        # stores path to an mo file which is supposed to be used as a guess
        # In case of open shell tm calculation this can be a tuple of files
        self.mo_guess: Path | str | tuple[str | Path, str | Path] | None = None

        self.from_part: str = ""

        self.charge: int = charge
        self.unpaired: int = unpaired

        # stores all flags for the jobtypes
        self.flags: dict[str, str] = {}


def set_omp_tmproc(jobs: list[ParallelJob], balance: bool, omp: int, ncores: int):
    """
    Configures the number of cores for TURBOMOLE processor or when automatic load balancing ('balancing 2.0')
    is disabled.
    """
    if balance:
        logger.warning(
            "Load balancing 2.0 is not supported for TURBOMOLE. Falling back to old behaviour."
        )
        # Set the number of cores per process to the (effective) minimum because this is expected to be most efficient
        omp = ncores // len(jobs) if len(jobs) < ncores // OMPMIN else OMPMIN
        ENVIRON["PARA_ARCH"] = "SMP"
        ENVIRON["PARNODES"] = str(omp)
        for job in jobs:
            job.omp = omp
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
        for job in jobs:
            job.omp = omp


def set_omp(jobs: list[ParallelJob], balance: bool, omp: int, ncores: int):
    """
    Sets the number of cores to be used for every job.
    """
    if balance:
        jobs_left, tot_jobs = len(jobs), len(jobs)

        # Maximum/minimum number of parallel processes
        maxprocs = ncores // OMPMIN
        minprocs = max(1, ncores // OMPMAX)

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
    omp: int,
    ncores: int,
    from_part: str,
    balance: bool = True,
    copy_mo: bool = True,
) -> list[ParallelJob]:
    """
    Prepares the jobs from the conformers data.
    """
    # Create ParallelJob instances from the conformers, sharing the prepinfo dict
    jobs = [ParallelJob(conf.geom, conf.charge, conf.unpaired) for conf in conformers]
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

    # Check if the the execution uses TM, because here OMP cannot be assigned on a job-variable basis
    if balance and prog != "tm":
        set_omp(jobs, balance, omp, ncores)
    # Similar steps will be taken if balancing 2.0 is disabled
    else:
        set_omp_tmproc(jobs, balance, omp, ncores)

    return jobs


def execute[T: QmResult](
    conformers: list[MoleculeData],
    task: Callable[..., tuple[T, MetaData]],
    job_config: XTBJobConfig | SPJobConfig,
    prog: QmProg | Literal["xtb"],
    ncores: int,
    omp: int,
    from_part: str,
    balance: bool = True,
    copy_mo: bool = True,
) -> tuple[dict[str, T], list[MoleculeData]]:
    """
    Executes the parallel tasks using a managed environment.

    Args:
        conformers (list[MoleculeData]): List of conformers for which jobs will be created and executed.
        task (Callable): Callable to be mapped onto the list of jobs created from conformers.
        job_config (XTBJobConfig | SPJobConfig): instructions for the execution of the task.
        prog (QmProg | "xtb"): Name of the QM program.
        balance (bool, optional): Whether to balance the number of cores used per job. Defaults to True.
        copy_mo (bool, optional): Whether to store the paths to the MO files for reuse. Defaults to True.

    Returns:
        dict[str, QmResult]: Job results.
        list[MoleculeData]: List of failed conformers.
    """
    # Initialize lists to store failed conformers and results
    failed_confs: list[MoleculeData] = []
    results: dict[str, T] = {}

    # Create a managed context to be able to share the number of free cores between processes
    logger.debug("Setting up parallel environment...")
    try:
        with setup_managers(ncores // OMPMIN, ncores) as (executor, _, resources):
            # Prepare jobs for parallel execution by assigning number of cores per job etc.
            jobs: list[ParallelJob] = prepare_jobs(
                conformers,
                prog,
                ncores,
                omp,
                from_part,
                balance=balance,
                copy_mo=copy_mo,
            )

            # Execute the jobs
            tasks: list[Future[tuple[T, MetaData]]] = []
            for job in jobs:
                tasks.append(executor.submit(task, job, job_config, resources))

            logger.debug("Waiting for jobs to complete...")
            try:
                completed = as_completed(tasks)

                # Process results as they come in
                logger.debug("Iterating over results...")
                for completed_task in completed:
                    conf_name = None
                    try:
                        # Get result or raise exception if task failed
                        result, meta = completed_task.result()
                        conf_name = meta.conf_name

                        # Find the corresponding conformer
                        conf = next(c for c in conformers if c.name == conf_name)

                        # Check for success
                        if not meta.success:
                            failed_confs.append(conf)
                            logger.warning(
                                f"{meta.conf_name} job failed. Error: {meta.error}. Check output files."
                            )
                        else:
                            results[meta.conf_name] = result
                            if prog in QmProg:
                                conf.mo_paths[prog].append(result.mo_path)

                    except Exception as e:
                        # Cancel remaining jobs
                        for t in tasks:
                            if not t.done():
                                t.cancel()

                        tb = traceback.format_exc()

                        # If we don't know which conformer failed, try to figure it out
                        # by finding the index of the task in the tasks list
                        if conf_name is None:
                            idx = tasks.index(completed_task)
                            if 0 <= idx < len(jobs):
                                conf_name = jobs[idx].conf.name

                        # Record the error for this conformer
                        logger.debug(f"Job failed for {conf_name}:\n{tb}")
                        raise RuntimeError(
                            "One or more tasks raised an exception. Set loglevel to DEBUG for details."
                        )

            except KeyboardInterrupt:
                logger.warning(
                    "Received keyboard interrupt. Cancelling remaining tasks..."
                )
                for t in tasks:
                    if not t.done():
                        t.cancel()
                # Re-raise to allow proper program termination
                raise

    except Exception as e:
        # Handle setup errors
        logger.error(f"Error setting up parallel execution environment: {str(e)}")
        # Re-raise with clearer context
        raise RuntimeError(
            f"Failed to set up parallel execution environment: {str(e)}"
        ) from e

    # Summarize results
    if len(failed_confs) == len(conformers):
        # Provide more detailed error information when all jobs fail
        raise RuntimeError(
            f"All jobs failed to execute. Please check your setup and output files.\n"
        )

    if len(failed_confs) > 0:
        logger.info(f"Number of failed jobs: {len(failed_confs)}.")
    else:
        logger.info("All jobs executed successfully.")

    return results, failed_confs
