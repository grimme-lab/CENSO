from pathlib import Path
from typing import Union

from ..molecules import GeometryData


class JobContext:
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
        self.mo_guess: Union[
            Path, str, tuple[Union[str, Path], Union[str, Path]], None
        ] = None

        self.from_part: str = ""

        self.charge: int = charge
        self.unpaired: int = unpaired

        # stores all flags for the jobtypes
        self.flags: dict[str, str] = {}
