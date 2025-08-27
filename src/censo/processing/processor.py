"""
Generic processor class. Only abstract methods: sp, opt (single-point and geometry optimization).
Also implements important methods inherited by all other processors, e.g. making external program calls and creating job dirs.
"""

from abc import abstractmethod
import functools
import os
from pathlib import Path
import subprocess
from collections.abc import Callable
from threading import Event
import time
from typing import final


from ..config.job_config import (
    SPJobConfig,
    XTBJobConfig,
    OptJobConfig,
)
from ..parallel import (
    ResourceMonitor,
    ParallelJob,
)
from ..config.job_config import (
    OptResult,
    QmResult,
    MetaData,
    SPResult,
)
from ..params import (
    WARNLEN,
    ENVIRON,
    Prog,
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

    # TODO: find better way to annotate this more clearly
    @staticmethod
    @final
    def _run[
        T: QmResult,
        U: XTBJobConfig | SPJobConfig,
    ](
        f: Callable[
            ...,
            tuple[T, MetaData],
        ],
    ) -> Callable[..., tuple[T, MetaData]]:
        """
        Wrapper function to manager resources and create job directory.
        Takes a callable as input which returns a tuple (results and metadata),
        and returns a the wrapped function.
        """

        @functools.wraps(f)
        def wrapper(
            self,
            job: ParallelJob,
            job_config: U,
            resources: ResourceMonitor,
            **kwargs,
        ) -> tuple[T, MetaData]:
            with resources.occupy_cores(job.omp):
                jobtype = f.__name__
                jobdir = Path(self._workdir) / job.conf.name / jobtype
                jobdir.mkdir(exist_ok=True, parents=True)
                logger.info(
                    f"{f'worker{os.getpid()}:':{WARNLEN}}Running "
                    + f"{jobtype} calculation using {self.__class__.__name__} in {jobdir} on {job.omp} cores using {self.progname.value}."
                )
                printf(
                    f"Running {jobtype} calculation for {job.conf.name} using {self.progname.value}."
                )
                result, meta = f(self, job, jobdir, job_config, **kwargs)
            return result, meta

        return wrapper

    def __init__(self, workdir: Path):
        """QM processor base class containing only xtb-related functions."""
        self._workdir: Path = workdir
        self.stop_event: Event | None = None

    @final
    def _make_call(
        self, call: list[str], outputpath: str, jobdir: str | Path
    ) -> tuple[int, str]:
        """
        Make a call to an external program and write output into outputfile.

        Args:
            call (list): list containing the call args to the external program
            outputpath (str): path to the outputfile
            jobdir (str): path to the jobdir

        Returns:
            returncode (int): returncode of the external program
            errors (str): stderr output
        """
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
                env=ENVIRON,
            )

            logger.debug(
                f"{f'worker{os.getpid()}:':{WARNLEN}}Started (PID: {sub.pid})."
            )

            # TODO: is this really required?
            # make sure to send SIGTERM to subprocess if program is quit
            # signal.signal(
            #     signal.SIGTERM, lambda signum, frame: handle_sigterm(signum, frame, sub)
            # )

            if self.stop_event is not None:
                while sub.poll() is None:
                    if self.stop_event.is_set():
                        logger.info(
                            f"{f'worker{os.getpid()}:':{WARNLEN}}Terminating subprocess {sub.pid}."
                        )
                        sub.terminate()
                        break
                    time.sleep(1)

            # wait for process to finish
            _, errors = sub.communicate()
            errors = errors.decode(errors="replace")
            returncode = sub.returncode
            logger.debug(
                f"{f'worker{os.getpid()}:':{WARNLEN}} Returncode: {returncode} Errors:\n{errors}"
            )

            # unregister SIGTERM handler
            # signal.signal(signal.SIGTERM, signal.SIG_DFL)

        logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Done.")

        return returncode, errors

    @final
    def _get_sym_num(self, sym: str | None = None, linear: bool = False) -> int:
        """Get rotational symmetry number from Schoenfließ symbol"""
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
    @_run
    def sp(
        self, job: ParallelJob, jobdir: Path | str, config: SPJobConfig, **kwargs
    ) -> tuple[SPResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @_run
    def opt(
        self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        raise NotImplementedError
