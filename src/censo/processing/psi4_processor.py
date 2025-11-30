"""
Contains Psi4Proc class for calculating psi4 related properties of conformers.
"""

import os
import pathlib
import typing

from ..config.job_config import (
    SPJobConfig,
)
from .job import (
    JobContext,
)
from .results import (
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
        Prepares an OrderedDict to be fed into the parser in order to write an input file for jobtype 'jobtype'
        (e.g. sp).

        Can load up a template file from user assets folder.

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

        inp[:0] = self.__prep_main(jobtype, config, no_solv)

        # prepare all options that are supposed to be placed before the
        # geometry definition
        inp[1:1] = self.__prep_pregeom(jobtype, config, job.omp)

        # prepare the geometry
        if template:
            geom_line = next(inp.index(line) for line in inp if "{geom}" in line)
            inp.pop(geom_line)
            inp[geom_line:geom_line] = self.__prep_geom(
                job.conf,
                xyzfile,
                job.charge,
                job.unpaired,
            )
        else:
            inp.extend(
                self.__prep_geom(
                    job.conf,
                    xyzfile,
                    job.charge,
                    job.unpaired,
                )
            )

        inp.extend(self.__prep_postgeom(config, jobtype))

        # Finally append newlines where needed
        return [line + "\n" if not line.endswith("\n") else line for line in inp]

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
            inp = []

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
            inp.append(f'print("%f" % energy("{FUNCTIONALS[config.func][Prog.PSI4.value]}"))')
            inp.append('')

            logger.warn("Grid is not implemnented for psi4")

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

        logger.warn("psi4 is not fully implemented")

        return self.__sp(*args, **kwargs)


Factory.register_builder(Prog.PSI4, Psi4Proc)
