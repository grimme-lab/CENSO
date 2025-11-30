"""
Contains Psi4Proc class for calculating psi4 related properties of conformers.
"""

import os
import pathlib
import typing

from ..config.job_config import (
    OptJobConfig,
    SPJobConfig,
)
from .job import (
    JobContext,
)
from .results import (
    OptResult,
    SPResult,
    MetaData,
)
from ..utilities import Factory
from ..logging import setup_logger
from ..params import (
    Prog,
)
from ..assets import FUNCTIONALS
from .qm_processor import QmProc

logger = setup_logger(__name__)


@typing.final
class Psi4Proc(QmProc):

    def __prep(
        self,
        job: JobContext,
        config: SPJobConfig,
        jobtype: str,
        no_solv: bool = False,
        xyzfile: str | None = None,
    ) -> list[str]:
        """
        Prepares an list of str to be written an input file for jobtype 'jobtype'
        (e.g. sp).

        TODO: use a template file from user assets folder.

        :param job: JobContext object containing the job information
        :type job: JobContext
        :param config: Configuration for the job
        :type config: SPJobConfig
        :param jobtype: jobtype to prepare the input for
        :type jobtype: str
        :param no_solv: if True, no solvent model is used
        :type no_solv: bool
        :param xyzfile: if not None, the geometry is read from this file instead of the job object
        :type xyzfile: str | None
        :returns: List of strings for the input file
        :rtype: list[str]

        NOTE: the xyzfile has to already exist, this function just bluntly writes the name into the input.
        """

        inp: list[str] = []

        template = "template" in config.model_fields and bool(config.template)
        if template:
            logger.warn("template not implemented for psi4")

        # basis set
        inp.append(f"set basis {config.basis}")
        inp.append('')

        # molecule
        inp.append("molecule {")
        for line in job.conf.toxyz()[2:]:
            inp.append(f"\t{line.strip()}")
        inp.append("}")
        inp.append('')

        # functional
        match jobtype:
            case "sp":
                inp.append(f'print("%f" % energy("{FUNCTIONALS[config.func][Prog.PSI4.value]}"))')
            case "opt":
                inp.append(f'print("%f" % optimize("{FUNCTIONALS[config.func][Prog.PSI4.value]}"))')
            case unknown:
                logger.warn(f"{unknown} is not implemented for psi4")
        inp.append('')

        logger.warn("Grid is not implemnented for psi4")

        return inp

    @typing.final
    def __sp(
        self,
        job: JobContext,
        config: SPJobConfig,
        jobdir: str | pathlib.Path | None = None,
        filename: str = "sp",
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[SPResult, MetaData]:
        """
        PSI4 single-point calculation.
        Unwrapped function to call from other methods.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :type job: JobContext
        :param config: Configuration for the job
        :type config: SPJobConfig
        :param jobdir: path to the job directory
        :type jobdir: str | Path | None
        :param filename: name of the input file
        :type filename: str
        :param no_solv: if True, no solvent model is used
        :type no_solv: bool
        :param prep: if True, a new input file is generated
        :type prep: bool
        :returns: Tuple of (SPResult, MetaData)
        :rtype: tuple[SPResult, MetaData]
        """
        if jobdir is None:
            jobdir = self._setup(job, "sp")

        # set results
        result = SPResult()
        meta = MetaData(job.conf.name)

        # set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # check for unsupported config
        solvation = not (config.gas_phase or no_solv)
        if solvation:
            logger.warning("Solvation is not implemented for psi4")
        if job.omp > 1:
            logger.warning("omp is not implemented for psi4")
        if config.copy_mo:
            logger.warning("copy_mo is not implemented for psi4")

        # prepare input
        if prep:
            inp = self.__prep(job, config, "sp", solvation)

            # write input into file "{filename}.inp" in a subdir
            # created for the conformer
            with open(inputpath, "w") as f:
                f.write('\n'.join(inp))

        # call psi4
        call = [config.paths.psi4, f"{filename}.inp"]
        returncode, _ = self._make_call(call, outputpath, jobdir)
        meta.success = returncode == 0

        # read output
        with open(outputpath) as out:
            lines = out.readlines()

        # Get final energy
        try:
            result.energy = next(
                (
                    float(line.split()[3])
                    for line in lines
                    if "Total Energy = " in line
                ),
            )
        except StopIteration:
            meta.success = False
            meta.error = "Could not parse final energy"

        return result, meta

    @typing.override
    def sp(self, *args, **kwargs):
        """
        Perform single-point calculation.

        :param args: Arguments.
        :param kwargs: Keyword arguments.
        :return: Tuple of (SP result, metadata).
        """

        logger.warn("sp is not fully implemented for psi4")

        return self.__sp(*args, **kwargs)

    @typing.override
    def opt(
        self,
        job: JobContext,
        config: OptJobConfig,
    ) -> tuple[OptResult, MetaData]:
        """
        Geometry optimization using psi4 optimizer.
        Note that solvation in handled here always implicitly.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :param config: Optimization configuration
        :return: Tuple of (OptResult, MetaData)
        """

        logger.warn("opt is not fully implemented for psi4")

        # prepare result
        result = OptResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "opt")
        filename = "opt"

        # set orca input/output paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        inp = self.__prep(job, config, "opt", no_solv=config.gas_phase)

        # write input into file "{filename}.inp" in a subdir
        # created for the conformer
        with open(inputpath, "w") as f:
            f.write('\n'.join(inp))

        # read output
        with open(outputpath) as out:
            lines = out.readlines()

        # Get final energy
        try:
            result.energy = next(
                (
                    float(line.split()[3])
                    for line in lines
                    if "Total Energy = " in line
                ),
            )
        except StopIteration:
            meta.success = False
            meta.error = "Could not parse final energy"

        return result, meta


Factory.register_builder(Prog.PSI4, Psi4Proc)
