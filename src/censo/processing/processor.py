"""
Generic processor class. Only abstract methods: sp, opt (single-point and geometry optimization).
Also implements important methods inherited by all other processors, e.g. making external program calls and creating job dirs.
"""

from abc import abstractmethod
import os
from pathlib import Path
import subprocess
import time
from typing import final
from dask.distributed import Variable, Lock


from ..config.job_config import (
    SPJobConfig,
    OptJobConfig,
)
from .job import (
    JobContext,
)
from .results import (
    OptResult,
    MetaData,
    SPResult,
)
from ..params import (
    WARNLEN,
    Prog,
    ENVIRON,
)
from ..utilities import printf
from ..logging import setup_logger

logger = setup_logger(__name__)


class GenericProc:
    """
    Generic processor class
    """

    progname: Prog = Prog.GENERIC

    # rotational entropy from symmetry
    # https://cccbdb.nist.gov/thermo.asp
    _rot_sym_num = {
        "c1": 1,
        "ci": 1,
        "cs": 1,
        "c2": 2,
        "c3": 3,
        "c4": 4,
        "c5": 5,
        "c6": 6,
        "c7": 7,
        "c8": 8,
        "c9": 9,
        "c10": 10,
        "c11": 11,
        "s4": 2,
        "s6": 3,
        "s8": 4,
        "d2": 4,
        "d3": 6,
        "d4": 8,
        "d5": 10,
        "d6": 12,
        "d7": 14,
        "d8": 16,
        "d9": 18,
        "d10": 20,
        "t": 12,
        "th": 12,
        "td": 12,
        "o": 24,
        "oh": 24,
        "ih": 60,
    }

    # @staticmethod
    # @final
    # def _run(f):
    #     """
    #     Wrapper function to manage resources and create job directory.
    #
    #     :param f: Callable that returns a tuple (results and metadata).
    #     :return: The wrapped function.
    #     """
    #
    #     @functools.wraps(f)
    #     def wrapper(
    #         self,
    #         job: JobContext,
    #         job_config,
    #         **kwargs,
    #     ):
    #         jobtype = f.__name__
    #         jobdir = Path(self._workdir) / job.conf.name / jobtype
    #         jobdir.mkdir(exist_ok=True, parents=True)
    #
    #         logger.info(
    #             f"{f'worker{os.getpid()}:':{WARNLEN}}Running "
    #             + f"{jobtype} calculation using {self.__class__.__name__} in {jobdir} on {job.omp} cores using {self.progname.value}."
    #         )
    #         printf(
    #             f"Running {jobtype} calculation for {job.conf.name} using {self.progname.value}."
    #         )
    #
    #         result, meta = f(self, job, jobdir, job_config, **kwargs)
    #
    #         return result, meta
    #
    #     return wrapper
    #
    @final
    def _setup(self, job: JobContext, jobtype: str):
        """
        Setup function to create job directory and print info before task execution.

        :param job: Job context
        :type job: JobContext
        :param jobtype: Type of job
        :type jobtype: str
        :returns: Path to the job directory
        :rtype: Path
        """
        jobdir = Path(self._workdir) / job.conf.name / jobtype
        jobdir.mkdir(exist_ok=True, parents=True)

        logger.info(
            f"{f'worker{os.getpid()}:':{WARNLEN}}Running "
            + f"{jobtype} calculation using {self.__class__.__name__} in {jobdir} on {job.omp} cores using {self.progname.value}."
        )
        printf(
            f"Running {jobtype} calculation for {job.conf.name} using {self.progname.value}."
        )

        return jobdir

    def __init__(self, workdir: Path):
        """
        QM processor base class containing only xtb-related functions.

        :param workdir: Working directory.
        """
        self._workdir: Path = workdir
        self._lock: Lock = Lock()
        self.cancel_var: Variable | None = None

    @final
    def _make_call(
        self,
        call: list[str],
        outputpath: str,
        jobdir: str | Path,
        env: dict[str, str] | None = None,
    ) -> tuple[int, str]:
        """
        Make a call to an external program and write output into outputfile.

        :param call: List containing the call args to the external program.
        :param outputpath: Path to the outputfile.
        :param jobdir: Path to the jobdir.
        :return: Tuple of (returncode of the external program or -1 in case of exception, stderr output).
        """
        try:
            # call external program and write output into outputfile
            with open(outputpath, "w", newline=None) as outputfile:
                logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running {call}...")

                # create subprocess for external program
                sub = subprocess.Popen(
                    call,
                    shell=False,
                    stderr=subprocess.PIPE,
                    cwd=jobdir,
                    stdout=outputfile,
                    env=env or ENVIRON,
                )

                logger.debug(
                    f"{f'worker{os.getpid()}:':{WARNLEN}}Started (PID: {sub.pid})."
                )

                # Cooperative cancellation using dask Variable
                while sub.poll() is None:
                    with self._lock:
                        if self.cancel_var and self.cancel_var.get():
                            logger.info(
                                f"{f'worker{os.getpid()}:':{WARNLEN}}Cancelling subprocess {sub.pid} due to cancellation signal."
                            )
                            sub.terminate()
                            break
                    time.sleep(0.1)

                # wait for process to finish
                _, stderr = sub.communicate()
                errors = stderr.decode(errors="replace")
                returncode = sub.returncode
                if returncode != 0:
                    logger.info(
                        f"{f'worker{os.getpid()}:':{WARNLEN}} Returncode: {returncode} Errors:\n{errors}"
                    )

            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Done.")

            return returncode, errors
        except Exception as e:
            logger.error(f"{f'worker{os.getpid()}:':{WARNLEN}}Error: {e}")
            return -1, str(e)

    @final
    def _get_sym_num(self, sym: str | None = None, linear: bool = False) -> int:
        """
        Get rotational symmetry number from Schoenfließ symbol.

        :param sym: Symmetry symbol.
        :param linear: Whether the molecule is linear.
        :return: Symmetry number.
        """
        if sym is None:
            sym = "c1"
        symnum = 1
        if linear and "c" in sym.lower()[0]:
            symnum = 1
            return symnum
        elif linear and "d" in sym.lower()[0]:
            symnum = 2
            return symnum
        for key in self._rot_sym_num:
            if key in sym.lower():
                symnum = self._rot_sym_num.get(key, 1)
                break
        return symnum

    @abstractmethod
    def sp(self, job: JobContext, config: SPJobConfig) -> tuple[SPResult, MetaData]:
        """
        Perform single-point calculation.

        :param job: Parallel job.
        :param config: SP configuration.
        :return: Tuple of (SP result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    def opt(self, job: JobContext, config: OptJobConfig) -> tuple[OptResult, MetaData]:
        """
        Perform geometry optimization.

        :param job: Parallel job.
        :param config: Optimization configuration.
        :return: Tuple of (optimization result, metadata).
        """
        raise NotImplementedError
