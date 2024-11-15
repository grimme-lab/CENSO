import multiprocessing
import signal

from .datastructure import ParallelJob, MoleculeData
from .logging import setup_logger
from .params import Config
from .tm_processor import TmProc
from .utilities import Factory
from .qm_processor import QmProc

logger = setup_logger(__name__)


class ParallelExecutor:
    """
    Manages parallel execution of external program calls. Handles the setup, execution, and cleanup of parallel tasks.
    """

    # NOTE: it might be possible to merge QmProc and this into a class to manage everything holistically

    def __init__(
        self,
        ncores: int,
        prog: str,
        workdir: str,
    ):
        """
        Initializes the ParallelExecutor with the given parameters.

        Args:
            ncores (int): Number of available cores.
            workdir (str): Working directory.
            prog (str): Name of the program to be used.
        """
        self.__MAXCORES = ncores
        self.__workdir = workdir
        self.__prog = prog

        self.__jobs: list[ParallelJob] = []

        self.__ncores = ncores
        self.__processor: QmProc = Factory.create(
            self.__prog, self.__workdir, self.__ncores
        )

    def __call__(self) -> None:
        """
        Executes the parallel tasks.
        """
        with multiprocessing.Pool(
            processes=self.__ncores // min([job.omp for job in self.__jobs] + [1])
        ) as pool:
            # Make sure SIGTERM is propagated
            # FIXME: will this raise an exception an break the context manager?
            signal.signal(
                signal.SIGTERM,
                lambda sig, frame: self.__handle_sigterm(sig, frame, pool),
            )

            # Register tasks in the pool
            # NOTE: this processes tasks synchronously so that the caller is blocked until execution is completed for all tasks (async would continue)
            pool.map(self.__processor.run, self.__jobs)

        logger.debug(f"Finished execution of {len(self.__jobs)} jobs.")

    @property
    def jobs(self) -> list[ParallelJob]:
        return self.__jobs

    @jobs.setter
    def jobs(self, jobs: list[ParallelJob]):
        assert all(isinstance(job, ParallelJob) for job in jobs)
        self.__jobs = jobs

    @property
    def failed_jobs(self) -> list[ParallelJob]:
        return [
            job
            for job in self.__jobs
            if any(not job.meta[jt]["success"] for jt in job.jobtype)
        ]

    def prepare_jobs(
        self,
        conformers: list[MoleculeData],
        jobtype: list[str],
        prepinfo: dict,
        balance: bool = True,
        copy_mo: bool = True,
    ):
        """
        Prepares the jobs from the conformers data.
        """
        self.__processor.copy_mo = copy_mo

        # Create ParallelJob instances from the conformers, sharing the prepinfo dict
        jobs = [ParallelJob(conf.geom, jobtype) for conf in conformers]
        for job in jobs:
            job.prepinfo.update(prepinfo)

            # Also insert mo guess file path
            if self.__prog == "orca":
                job.mo_guess = next(
                    (
                        conf.mo_paths
                        for conf in conformers
                        for mo_path in conf.mo_paths
                        if conf.name == job.conf.name and ".gbw" in mo_path
                    ),
                    None,
                )
            elif self.__prog == "tm":
                job.mo_guess = next(
                    (
                        conf.mo_paths
                        for conf in conformers
                        for mo_path in conf.mo_paths
                        if conf.name == job.conf.name
                        and any(kw in mo_path for kw in ["alpha", "beta", "mos"])
                    ),
                    None,
                )

        self.__jobs = jobs

        # Check if the the execution uses TM, because here OMP cannot be assigned on a job-variable basis
        if balance and not isinstance(self.__processor, TmProc):
            self.__set_omp_chunking()
        # Similar steps will be taken if balancing 2.0 is disabled
        else:
            self.__set_omp_tmproc()

    def __set_omp_tmproc(self):
        """
        Configures the number of cores for TURBOMOLE processor or when automatic load balancing ('balancing 2.0')
        is disabled.
        """
        if self.__balance and isinstance(self.__processor, TmProc):
            logger.warning(
                "Load balancing 2.0 is not supported for TURBOMOLE. Falling back to old behaviour."
            )
            # Set the number of cores per process to the (effective) minimum because this is expected to be most efficient
            omp = (
                Config.NCORES // len(self.__jobs)
                if len(self.__jobs) < Config.NCORES // Config.OMPMIN
                else Config.OMPMIN
            )
            Config.ENVIRON["PARA_ARCH"] = "SMP"
            Config.ENVIRON["PARNODES"] = str(omp)
            for job in self.__jobs:
                job.omp = omp
        else:
            omp = Config.OMP
            if omp < Config.OMPMIN:
                logger.warning(
                    f"User OMP setting is below the minimum value of {Config.OMPMIN}. Using {Config.OMPMIN} instead."
                )
                omp = Config.OMPMIN
            elif omp > Config.NCORES:
                logger.warning(
                    f"Value of {omp} for OMP is larger than the number of available cores {Config.NCORES}. Using OMP = {Config.NCORES}."
                )
                omp = Config.NCORES
            for job in self.__jobs:
                job.omp = omp

    def __set_omp_chunking(self):
        """
        Sets the number of cores to be used for every job.
        """
        jobs_left, tot_jobs = len(self.__jobs), len(self.__jobs)

        # Maximum/minimum number of parallel processes
        maxprocs = Config.NCORES // Config.OMPMIN
        minprocs = max(1, Config.NCORES // Config.OMPMAX)

        while jobs_left > 0:
            # Find the optimal number of parallel processes to maximize core utilization while avoiding
            # to parallelize on too many cores (efficiency drops off quickly)
            p = (
                maxprocs
                if jobs_left >= maxprocs
                else max(
                    j
                    for j in range(minprocs, maxprocs)
                    if Config.NCORES % j == 0 and j <= jobs_left
                )
            )

            # Assign determined OMP value to jobs
            while jobs_left - p >= 0:
                for job in self.__jobs[tot_jobs - jobs_left : tot_jobs - jobs_left + p]:
                    job.omp = Config.NCORES // p
                jobs_left -= p

    def __handle_sigterm(self, signum: int, frame, pool: multiprocessing.Pool):
        """
        Handles SIGTERM signal to terminate gracefully.

        Args:
            signum (int): Signal number.
            frame (frame): Current stack frame.
        """
        logger.critical("Received SIGTERM. Terminating.")
        # Immediately terminate all subprocesses
        pool.terminate()


def execute(
    conformers,
    workdir,
    prog,
    prepinfo,
    jobtype,
    copy_mo=False,
    balance=True,
):
    """
    Executes the parallel tasks using the ParallelExecutor.

    Args:
        conformers (list[MoleculeData]): List of conformers for which jobs will be created and executed.
        workdir (str): Working directory.
        prog (str): Name of the program to be used.
        prepinfo (dict): Preparation information for the jobs.
        jobtype (list): List of job types.
        copy_mo (bool, optional): Whether to copy the MO-files from the previous calculation.
        balance (bool, optional): Whether to balance the number of cores used per job.

    Returns:
        list: List of results.
    """
    executor = ParallelExecutor(Config.NCORES, prog, workdir)

    # Prepare jobs for parallel execution by assigning number of cores per job etc.
    executor.prepare_jobs(
        conformers, jobtype, prepinfo, balance=balance, copy_mo=copy_mo
    )

    # Execute the jobs
    executor()

    # determine failed jobs
    logger.debug("Checking for failed jobs...")

    if len(executor.failed_jobs) > 0:
        logger.info(f"Number of failed jobs: {len(executor.failed_jobs)}.")

        """
    if retry_failed and len(failed_jobs) > 0:
        # create a new list of failed jobs that should be restarted with special flags
        # contains jobs that should be retried (depends on wether the error can be handled or not)
        # retry = []

        # determine flags for jobs based on error messages
        for failed_job in failed_jobs:
            handled_errors = ["scf_not_converged", "Previous calculation failed"]

            # list of jobtypes that should be removed from the jobtype list
            jtremove = []
            for jt in failed_job.jobtype:
                if not failed_job.meta[jt]["success"]:
                    if failed_job.meta[jt]["error"] in handled_errors:
                        retry.append(failed_job)
                        failed_job.flags[jt] = failed_job.meta[jt]["error"]
                # store all successful jobtypes to be removed later
                elif failed_job.meta[jt]["success"]:
                    jtremove.append(jt)

            # remove all successful jobs from jobtype to avoid re-execution
            for jt in jtremove:
                failed_job.jobtype.remove(jt)


        if len(retry) > 0:
            # execute jobs that should be retried
            executor.jobs = retry
            _ = executor.execute()
        else:
            logger.info("No failed jobs can be recovered.")
        """

        # any jobs that still failed will lead to the conformer being marked as unrecoverable
        failed_confs = [job.conf.name for job in executor.failed_jobs]
        for job in executor.failed_jobs:
            errors = {jt: job.meta[jt]["error"] for jt in job.meta}
            logger.warning(
                f"{job.conf.name} job failed. Errors: {errors}. Check output files."
            )
    else:
        logger.debug("All jobs executed successfully.")

    return {job.conf.name: job.results for job in executor.jobs}, failed_confs
