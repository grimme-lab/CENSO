"""
Contains QmProc base class,
Additionally contains functions which should be present irrespective of the QM
code. (xTB always available)
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
    GsolvResult,
    NMRResult,
    OptResult,
    MetaData,
    ParallelJob,
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
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def gsolv(
        self, job: ParallelJob, jobdir: Path | str, config: SPJobConfig, **kwargs
    ) -> tuple[GsolvResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def xtb_opt(
        self, job: ParallelJob, jobdir: Path | str, config: XTBOptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def opt(
        self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def nmr(
        self, job: ParallelJob, jobdir: Path | str, config: NMRJobConfig, **kwargs
    ) -> tuple[NMRResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @GenericProc._run
    def uvvis(
        self, job: ParallelJob, jobdir: Path | str, config: UVVisJobConfig, **kwargs
    ) -> tuple[UVVisResult, MetaData]:
        raise NotImplementedError
