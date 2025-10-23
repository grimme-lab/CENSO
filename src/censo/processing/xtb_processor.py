"""
Processor for xtb calculations.
"""

import json
import os
from pathlib import Path
from typing import final, override


from ..config.job_config import (
    RRHOJobConfig,
    XTBJobConfig,
)
from .job import (
    JobContext,
)
from .results import (
    GsolvResult,
    MetaData,
    RRHOResult,
    SPResult,
)
from ..utilities import Factory
from ..logging import setup_logger
from ..assets import SOLVENTS
from .processor import GenericProc
from ..params import Prog

logger = setup_logger(__name__)


@final
class XtbProc(GenericProc):
    """
    XtbProcessor
    """

    progname = Prog.XTB

    def _sp(
        self,
        job: JobContext,
        config: XTBJobConfig,
        jobdir: str | Path | None = None,
        filename: str = "xtb_sp",
        no_solv: bool = False,
    ) -> tuple[SPResult, MetaData]:
        """
        Calculates the single-point energy with GFNn-xTB or GFN-FF.
        Unwrapped function to call from other methods.

        :param job: job to run
        :type job: JobContext
        :param config: XTB configuration
        :type config: XTBJobConfig
        :param jobdir: path to the jobdir
        :type jobdir: str | Path | None
        :param filename: filename to use for the coord file. Defaults to "xtb_sp".
        :type filename: str
        :param no_solv: whether to run the sp in gas-phase. Defaults to False.
        :type no_solv: bool
        :returns: Tuple of (SPResult, MetaData)
        :rtype: tuple[SPResult, MetaData]
        """
        if jobdir is None:
            jobdir = self._setup(job, "xtb_sp")

        result = SPResult()
        meta = MetaData(job.conf.name)

        # set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.coord")
        outputpath = os.path.join(jobdir, f"{filename}.out")
        xcontrolname = "xtb_sp-xcontrol-inp"

        # cleanup
        files = [
            "xtbrestart",
            "xtbtopo.mol",
            xcontrolname,
            "wbo",
            "charges",
            "gfnff_topo",
            f"{filename}.out",
        ]

        # remove potentially preexisting files to avoid confusion
        for file in files:
            if os.path.isfile(os.path.join(jobdir, file)):
                os.remove(os.path.join(jobdir, file))

        # generate coord file for xtb
        with open(inputpath, "w", newline=None) as f:
            f.writelines(job.conf.tocoord())

        # setup call for xtb single-point
        call: list[str] = [
            config.paths.xtb,
            f"{filename}.coord",
            "--" + config.gfnv,
            "--sp",
            "--chrg",
            f"{job.charge}",
            "--norestart",
            "--parallel",
            f"{job.omp}",
        ]

        # add solvent to xtb call if not a gas-phase sp
        if not no_solv and not config.gas_phase:
            assert config.solvent
            assert config.sm_rrho
            solvent_key = SOLVENTS[config.solvent][config.sm_rrho]
            call.extend(
                [
                    "--" + config.sm_rrho,
                    solvent_key,
                    "reference",
                    "-I",
                    xcontrolname,
                ]
            )

            # set gbsa grid
            with open(os.path.join(jobdir, xcontrolname), "w", newline=None) as xcout:
                xcout.writelines(["$gbsa\n", "  gbsagrid=tight\n", "$end\n"])

        # call xtb
        returncode, _ = self._make_call(call, outputpath, jobdir)

        # if returncode != 0 then some error happened in xtb
        # TODO: returncodes
        if returncode != 0:
            meta.success = False
            meta.error = "unknown_error"
            return result, meta

        # read energy from outputfile
        with open(outputpath) as outputfile:
            for line in outputfile.readlines():
                if "| TOTAL ENERGY" in line:
                    result.energy = float(line.split()[3])
                    meta.success = True
                    break
                if "convergence criteria cannot be satisfied" in line:
                    meta.success = False
                    meta.error = "scf_not_converged"
                    break

        return result, meta

    @override
    def sp(self, *args, **kwargs):
        """
        Wrapped version of _sp.

        :param args: Arguments.
        :param kwargs: Keyword arguments.
        :return: Tuple of (SPResult, MetaData).
        """
        return self._sp(*args, **kwargs)

    def gsolv(
        self, job: JobContext, config: XTBJobConfig
    ) -> tuple[GsolvResult, MetaData]:
        """
        Calculate additive GBSA or ALPB solvation using GFNn-xTB or GFN-FF.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :type job: JobContext
        :param config: XTB configuration
        :type config: XTBJobConfig
        :returns: Tuple of (GsolvResult, MetaData)
        :rtype: tuple[GsolvResult, MetaData]
        """
        result = GsolvResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "xtb_gsolv")

        # run gas-phase GFN single-point
        spres, spmeta = self._sp(
            job, config, jobdir=jobdir, filename="gas", no_solv=True
        )
        if spmeta.success:
            result.energy_gas = spres.energy
        else:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # run single-point in solution:
        spres, spmeta = self._sp(job, config, jobdir=jobdir, filename="solv")
        if spmeta.success:
            result.energy_solv = spres.energy
        else:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # only reached if both gas-phase and solvated sp succeeded
        result.gsolv = result.energy_solv - result.energy_gas
        meta.success = True

        return result, meta

    # TODO: break this down
    def xtb_rrho(
        self,
        job: JobContext,
        config: RRHOJobConfig,
    ) -> tuple[RRHOResult, MetaData]:
        """
        Calculates the mRRHO contribution to the free enthalpy of a conformer with GFNn-xTB/GFN-FF.

        :param job: job to run
        :type job: JobContext
        :param config: RRHO configuration
        :type config: RRHOJobConfig
        :returns: result of the rrho calculation and metadata
        :rtype: tuple[RRHOResult, MetaData]
        """
        # what is returned in the end
        result = RRHOResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "xtb_rrho")
        filename = "xtb_rrho"

        # set in/out path
        outputpath = os.path.join(jobdir, f"{filename}.out")
        xcontrolname = "rrho-xcontrol-inp"
        xcontrolpath = os.path.join(jobdir, xcontrolname)

        # TODO - is this list complete?
        files = [
            "xtbrestart",
            "xtbtopo.mol",
            xcontrolname,
            "wbo",
            "charges",
            "gfnff_topo",
        ]

        # remove potentially preexisting files to avoid confusion
        for file in files:
            if os.path.isfile(os.path.join(jobdir, file)):
                os.remove(os.path.join(jobdir, file))

        # setup xcontrol
        with open(xcontrolpath, "w", newline=None) as xcout:
            xcout.write("$thermo\n")
            # if config.multitemp:
            #     trange = frange(
            #         config.trange[0],
            #         config.trange[1],
            #         step=config.trange[2],
            #     )
            #
            #     # Always append the fixed temperature to the trange so that it is the last value
            #     # (important since --enso will make xtb give the G(T) value for this temperature)
            #     assert config.temperature
            #     trange.append(config.temperature)
            #
            #     # Write trange to the xcontrol file
            #     xcout.write(f"    temp=" + ",".join([str(i) for i in trange]) + "\n")
            #     if config.temperature in trange[:-1]:
            #         trange = trange[:-1]
            # else:
            xcout.write(f"    temp={config.temperature}\n")

            xcout.write(f"    sthr={config.sthr}\n")

            xcout.write(f"    imagthr={config.imagthr}\n")

            # TODO - when do you actually want to set the scale manually?
            # don't set scale manually for now
            """ if self.instructions.get("scale", "automatic") != "automatic":
                # for automatic --> is method dependant leave it to xTB e.g. GFNFF has a
                # different scaling factor than GFN2
                xcout.write(f"    scale={self.instructions['scale']}\n") """

            xcout.write("$symmetry\n")
            xcout.write("     maxat=1000\n")
            # always consider symmetry
            # xcout.write("    desy=0.1\n") # taken from xtb defaults
            # xcout.write("    desy=0.0\n")

            # set gbsa grid
            if not config.gas_phase:
                xcout.write("$gbsa\n")
                xcout.write("  gbsagrid=tight\n")

            xcout.write("$end\n")

        # if config.bhess:
        #     # set ohess or bhess
        #     dohess = "--bhess"
        # else:
        #     dohess = "--ohess"

        # generate coord file for xtb
        with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as f:
            f.writelines(job.conf.tocoord())

        call: list[str] = [
            config.paths.xtb,
            f"{filename}.coord",
            "--" + config.gfnv,
            "--bhess",
            "vtight",
            "--chrg",
            f"{job.charge}",
            "--enso",
            "--norestart",
            "-I",
            xcontrolname,
            "--parallel",
            f"{job.omp}",
        ]

        # add solvent to xtb call if necessary
        if not config.gas_phase:
            assert config.solvent
            assert config.sm_rrho
            solvent_key = SOLVENTS[config.solvent][config.sm_rrho]
            call.extend(
                [
                    "--" + config.sm_rrho,
                    solvent_key,
                ]
            )

        # TODO: if rmsd bias is used (should be defined in censo workdir (cwd))
        # if config.rmsdbias:
        #     cwd = os.getcwd()
        #
        #     call.extend(
        #         [
        #             "--bias-input",
        #             os.path.join(cwd, "rmsdpot.xyz"),
        #         ]
        #     )

        # call xtb
        returncode, errors = self._make_call(call, outputpath, jobdir)

        # check if converged:
        if returncode != 0:
            meta.success = False
            meta.error = "unknown_error"
            return result, meta

        # read output and store lines
        with open(outputpath) as outputfile:
            lines = outputfile.readlines()

        # if config.multitemp:
        #     # get gibbs energy, enthalpy and entropy for given temperature range
        #     # gibbs energy
        #     gt: dict[float, float] = {}
        #
        #     # enthalpy
        #     ht: dict[float, float] = {}
        #
        #     # Get Gibbs energy and enthalpy
        #     for i, line in enumerate(lines):
        #         if "T/K" in line:
        #             for line2 in lines[i + 2 :]:
        #                 if "----------------------------------" in line2:
        #                     break
        #
        #                 temp = float(line2.split()[0])
        #                 gt[temp] = float(line2.split()[4])
        #                 ht[temp] = float(line2.split()[2])
        #             break
        #
        #     # rotational entropy
        #     rotS: dict[float, float] = {}
        #
        #     # Get rotational entropy
        #     entropy_lines = [
        #         (line, lines[i + 1]) for i, line in enumerate(lines) if "VIB" in line
        #     ]
        #     for line in entropy_lines:
        #         temp = float(line[0].split()[0])
        #         rotS[temp] = float(line[1].split()[4])
        #
        #     # check if xtb calculated the temperature range correctly
        #     if not (
        #         len(trange) == len(gt)
        #         and len(trange) == len(ht)
        #         and len(trange) == len(rotS)
        #     ):
        #         meta.success = False
        #         meta.error = "xtb_trange_inconsistent"
        #         return result, meta
        #     else:
        #         result.gibbs = gt
        #         result.enthalpy = ht
        #         result.entropy = rotS

        # Extract symmetry
        result.linear = next(
            {"true": True, "false": False}[line.split()[2]]
            for line in lines
            if ":  linear? " in line
        )

        # Extract rmsd
        # if config.bhess:
        #     result.rmsd = next(
        #         (float(line.split()[3]) for line in lines if "final rmsd / " in line)
        #     )

        # xtb_enso.json is generated by xtb by using the '--enso' argument *only* when using --bhess or --ohess
        # (when a hessian is calculated)
        # contains output from xtb in json format to be more easily digestible by CENSO
        with open(
            os.path.join(jobdir, "xtb_enso.json"),
        ) as f:
            data = json.load(f)

        # read number of imaginary frequencies and print warning
        if "number of imags" in data:
            if data["number of imags"] > 0:
                logger.warning(
                    f"Found {data['number of imags']} significant"
                    + " imaginary frequencies for "
                    + f"{job.conf.name}."
                )

        # get gibbs energy
        if "G(T)" in data:
            meta.success = True

            temp = float(config.temperature)
            if config.temperature == 0.0:
                result.energy = data.get("ZPVE", 0.0)
                result.gibbs[temp] = data.get("ZPVE", 0.0)
                result.enthalpy[temp] = data.get("ZPVE", 0.0)
                # if not config.multitemp:
                result.entropy[temp] = 0.0  # set this to None for predictability (?)
            else:
                result.energy = data.get("G(T)", 0.0)
                result.gibbs[temp] = data.get("G(T)", 0.0)
                # if not config.multitemp:
                result.enthalpy[temp] = 0.0  # set this to None for predictability (?)
                result.entropy[temp] = 0.0  # set this to None for predictability (?)

            # only determine symmetry if all the needed information is there
            if "point group" and "linear" in data:
                result.symmetry = data["point group"]
                result.linear = data["linear"]
            else:
                # could not determine symmetry correctly
                result.symmetry = "c1"
                result.linear = False

            # calculate symnum
            result.symnum = self._get_sym_num(sym=result.symmetry, linear=result.linear)
        else:
            meta.success = False
            meta.error = "Could not read xtb_enso.json"

        return result, meta

    @override
    def opt(self, *args, **kwargs):
        raise NotImplementedError


Factory.register_builder(Prog.XTB, XtbProc)
