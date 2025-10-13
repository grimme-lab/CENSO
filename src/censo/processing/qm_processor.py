"""
Base class for all QM-based processors, e.g. ORCA, TM (or QChem, ...).
"""

from abc import abstractmethod
from pathlib import Path


from .processor import GenericProc
from ..config.job_config import (
    NMRJobConfig,
    SPJobConfig,
    UVVisJobConfig,
    OptJobConfig,
    XTBOptJobConfig,
)
from ..parallel import (
    ParallelJob,
)
from ..config.job_config import (
    GsolvResult,
    NMRResult,
    OptResult,
    MetaData,
    SPResult,
    UVVisResult,
)
from ..logging import setup_logger

logger = setup_logger(__name__)


class QmProc(GenericProc):
    """
    QmProc base class
    """

    @abstractmethod
    @GenericProc._run
    def sp(
        self, job: ParallelJob, jobdir: Path | str, config: SPJobConfig, **kwargs
    ) -> tuple[SPResult, MetaData]:
        """
        Perform single-point calculation.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: SP configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (SP result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def gsolv(
        self, job: ParallelJob, jobdir: Path | str, config: SPJobConfig, **kwargs
    ) -> tuple[GsolvResult, MetaData]:
        """
        Perform solvation energy calculation.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: SP configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (Gsolv result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def xtb_opt(
        self, job: ParallelJob, jobdir: Path | str, config: XTBOptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        """
        Perform xTB geometry optimization.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: XTB optimization configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (optimization result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def opt(
        self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        """
        Perform geometry optimization.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: Optimization configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (optimization result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def nmr(
        self, job: ParallelJob, jobdir: Path | str, config: NMRJobConfig, **kwargs
    ) -> tuple[NMRResult, MetaData]:
        """
        Perform NMR calculation.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: NMR configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (NMR result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def rot(self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs):
        """
        Perform rotational calculation.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: Optimization configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (rotational result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def uvvis(
        self, job: ParallelJob, jobdir: Path | str, config: UVVisJobConfig, **kwargs
    ) -> tuple[UVVisResult, MetaData]:
        """
        Perform UV-Vis calculation.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: UV-Vis configuration.
        :param kwargs: Additional arguments.
        :return: Tuple of (UV-Vis result, metadata).
        """
        raise NotImplementedError
