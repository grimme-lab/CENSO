"""
Base class for all QM-based processors, e.g. ORCA, TM (or QChem, ...).
"""

from abc import abstractmethod


from .processor import GenericProc
from ..config.job_config import (
    NMRJobConfig,
    SPJobConfig,
    UVVisJobConfig,
    RotJobConfig,
    XTBOptJobConfig,
)
from .job import (
    JobContext,
)
from .results import (
    GsolvResult,
    NMRResult,
    OptResult,
    MetaData,
    RotResult,
    UVVisResult,
)
from ..logging import setup_logger

logger = setup_logger(__name__)


class QmProc(GenericProc):
    """
    QmProc base class
    """

    @abstractmethod
    def gsolv(
        self, job: JobContext, config: SPJobConfig
    ) -> tuple[GsolvResult, MetaData]:
        """
        Perform solvation energy calculation.

        :param job: Parallel job.
        :param config: SP configuration.
        :return: Tuple of (Gsolv result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    def xtb_opt(
        self, job: JobContext, config: XTBOptJobConfig
    ) -> tuple[OptResult, MetaData]:
        """
        Perform xTB geometry optimization.

        :param job: Parallel job.
        :param config: XTB optimization configuration.
        :return: Tuple of (optimization result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    def nmr(self, job: JobContext, config: NMRJobConfig) -> tuple[NMRResult, MetaData]:
        """
        Perform NMR calculation.

        :param job: Parallel job.
        :param config: NMR configuration.
        :return: Tuple of (NMR result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    def rot(self, job: JobContext, config: RotJobConfig) -> tuple[RotResult, MetaData]:
        """
        Perform rotational calculation.

        :param job: Parallel job.
        :param config: Optimization configuration.
        :return: Tuple of (rotational result, metadata).
        """
        raise NotImplementedError

    @abstractmethod
    def uvvis(
        self, job: JobContext, config: UVVisJobConfig
    ) -> tuple[UVVisResult, MetaData]:
        """
        Perform UV-Vis calculation.

        :param job: Parallel job.
        :param jobdir: Job directory.
        :param config: UV-Vis configuration.
        :return: Tuple of (UV-Vis result, metadata).
        """
        raise NotImplementedError
