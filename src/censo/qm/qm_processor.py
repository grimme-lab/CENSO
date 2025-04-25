"""
Contains QmProc base class,
Additionally contains functions which should be present irrespective of the QM
code. (xTB always available)
"""

from abc import abstractmethod
import functools
import json
import os
from pathlib import Path
import subprocess
from collections.abc import Callable
from typing import final


from ..config.job_config import (
    RRHOJobConfig,
    SPJobConfig,
    XTBJobConfig,
    OptJobConfig,
    XTBOptJobConfig,
)
from ..parallel import (
    GsolvResult,
    NMRResult,
    OptResult,
    QmResult,
    MetaData,
    RRHOResult,
    ResourceMonitor,
    ParallelJob,
    SPResult,
    UVVisResult,
)
from ..params import (
    PLENGTH,
    DIGILEN,
    WARNLEN,
    ENVIRON,
)
from ..utilities import printf, frange, Factory
from ..logging import setup_logger
from ..assets import SOLVENTS

logger = setup_logger(__name__)


class QmProc:
    """
    QmProc base class
    """

    paths = {
        "orca": "",
        "orcaversion": "",
        "xtb": "",
        "cosmorssetup": "",
        "cosmotherm": "",
        "dbpath": "",
        "cosmothermversion": "",
        "mpshift": "",
        "escf": "",
    }

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
                jobdir = Path(self.workdir) / job.conf.name / jobtype
                jobdir.mkdir(exist_ok=True, parents=True)
                logger.info(
                    f"{f'worker{os.getpid()}:':{WARNLEN}}Running "
                    + f"{jobtype} calculation using {self.__class__.__name__} in {jobdir} on {job.omp} cores."
                )
                printf(f"Running {jobtype} calculation for {job.conf.name}.")
                result, meta = f(self, job, jobdir, job_config, **kwargs)
            return result, meta

        return wrapper

    @classmethod
    def print_paths(cls) -> None:
        """
        Print out the paths of all external QM programs.
        """
        # Create an empty list to store the lines of the output.
        lines = []

        # Append a separator line to the output.
        lines.append("\n" + "".ljust(PLENGTH, "-") + "\n")

        # Append the title of the section to the output, centered.
        lines.append("PATHS of external QM programs".center(PLENGTH, " ") + "\n")

        # Append a separator line to the output.
        lines.append("".ljust(PLENGTH, "-") + "\n")

        # Iterate over each program and its path in the settings.
        for program, path in cls.paths.items():
            # Append a line with the program and its path to the output.
            lines.append(f"{program}:".ljust(DIGILEN, " ") + f"{path}\n")

        # Print each line of the output.
        for line in lines:
            printf(line)

    def __init__(self, workdir: Path):
        """QM processor base class containing only xtb-related functions."""
        self.workdir: Path = workdir

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

            # wait for process to finish
            _, errors = sub.communicate()
            errors = errors.decode(errors="replace")
            returncode = sub.returncode

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

    @final
    def _xtb_sp(
        self,
        job: ParallelJob,
        jobdir: str | Path,
        config: XTBJobConfig,
        filename: str = "xtb_sp",
        no_solv: bool = False,
    ) -> tuple[SPResult, MetaData]:
        """
        Calculates the single-point energy with GFNn-xTB or GFN-FF.
        Unwrapped function to call from other methods.

        Args:
            job (ParallelJob): job to run
            jobdir (str): path to the jobdir
            filename (str, optional): filename to use for the coord file. Defaults to "xtb_sp".
            no_solv (bool, optional): whether to run the sp in gas-phase. Defaults to False.

        Returns:
            result (SPResult): result of the sp calculation
            meta (MetaData): metadata about the job

        """
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
        with open(inputpath, "w", newline=None) as file:
            file.writelines(job.conf.tocoord())

        # setup call for xtb single-point
        call: list[str] = [
            self.paths["xtb"],
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
        returncode, errors = self._make_call(call, outputpath, jobdir)

        # if returncode != 0 then some error happened in xtb
        # TODO: returncodes
        if returncode != 0:
            meta.success = False
            meta.error = "unknown_error"
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            return result, meta

        # read energy from outputfile
        with open(outputpath, "r") as outputfile:
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

    @final
    @_run
    def xtb_sp(self, *args, **kwargs):
        """Wrapped version of xtb_sp."""
        return self._xtb_sp(*args, **kwargs)

    @final
    @_run
    def xtb_gsolv(
        self, job: ParallelJob, jobdir: str | Path, config: XTBJobConfig
    ) -> tuple[GsolvResult, MetaData]:
        """
        Calculate additive GBSA or ALPB solvation using GFNn-xTB or GFN-FF.

        Args:
            job (ParallelJob): job to run
            jobdir (str): path to the jobdir

        Returns:
            result (GsolvResult): result of the gsolv calculation
            meta (MetaData): metadata about the calculation

        result = {
            "gsolv": None,
            "energy_xtb_gas": None,
            "energy_xtb_solv": None,
        }
        """
        result = GsolvResult()
        meta = MetaData(job.conf.name)

        # run gas-phase GFN single-point
        spres, spmeta = self._xtb_sp(job, jobdir, config, filename="gas", no_solv=True)
        if spmeta.success:
            result.energy_gas = spres.energy
        else:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # run single-point in solution:
        spres, spmeta = self._xtb_sp(job, jobdir, config, filename="solv")
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
    @final
    @_run
    def xtb_rrho(
        self,
        job: ParallelJob,
        jobdir: str | Path,
        config: RRHOJobConfig,
        filename: str = "xtb_rrho",
    ) -> tuple[RRHOResult, MetaData]:
        """
        Calculates the mRRHO contribution to the free enthalpy of a conformer with GFNn-xTB/GFN-FF.

        Args:
            job (ParallelJob): job to run
            jobdir (str): path to the jobdir
            filename (str, optional): filename to use for the coord file. Defaults to "xtb_rrho".

        Returns:
            result (dict[str, any]): result of the rrho calculation
            meta (dict[str, str | bool | None]): metadata about the calculation

        result = {
            "energy": None, # contains the gibbs energy at given temperature (might be ZPVE if T = 0K)
            "rmsd": None,
            "gibbs": None,
            "enthalpy": None,
            "entropy": None,
            "symmetry": None,
            "symnum": None,
        }
        """
        # what is returned in the end
        result = RRHOResult()
        meta = MetaData(job.conf.name)

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
            if config.multitemp:
                trange = frange(
                    config.trange[0],
                    config.trange[1],
                    step=config.trange[2],
                )

                # Always append the fixed temperature to the trange so that it is the last value
                # (important since --enso will make xtb give the G(T) value for this temperature)
                assert config.temperature
                trange.append(config.temperature)

                # Write trange to the xcontrol file
                xcout.write(f"    temp=" + ",".join([str(i) for i in trange]) + "\n")
                if config.temperature in trange[:-1]:
                    trange = trange[:-1]
            else:
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

        if config.bhess:
            # set ohess or bhess
            dohess = "--bhess"
        else:
            dohess = "--ohess"

        # generate coord file for xtb
        with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as file:
            file.writelines(job.conf.tocoord())

        call: list[str] = [
            self.paths["xtb"],
            f"{filename}.coord",
            "--" + config.gfnv,
            dohess,
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
        if config.rmsdbias:
            cwd = os.getcwd()

            call.extend(
                [
                    "--bias-input",
                    os.path.join(cwd, "rmsdpot.xyz"),
                ]
            )

        # call xtb
        returncode, errors = self._make_call(call, outputpath, jobdir)

        # check if converged:
        if returncode != 0:
            meta.success = False
            meta.error = "unknown_error"
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            return result, meta

        # read output and store lines
        with open(outputpath, "r") as outputfile:
            lines = outputfile.readlines()

        if config.multitemp:
            # get gibbs energy, enthalpy and entropy for given temperature range
            # gibbs energy
            gt: dict[float, float] = {}

            # enthalpy
            ht: dict[float, float] = {}

            # Get Gibbs energy and enthalpy
            for line in lines:
                if "T/K" in line:
                    for line2 in lines[lines.index(line) + 2 :]:
                        if "----------------------------------" in line2:
                            break

                        temp = float(line2.split()[0])
                        gt[temp] = float(line2.split()[4])
                        ht[temp] = float(line2.split()[2])
                break

            # rotational entropy
            rotS: dict[float, float] = {}

            # Get rotational entropy
            entropy_lines = [
                (line, lines[i + 1]) for i, line in enumerate(lines) if "VIB" in line
            ]
            for line in entropy_lines:
                temp = float(line[0].split()[0])
                rotS[temp] = float(line[1].split()[4])

            # check if xtb calculated the temperature range correctly
            if not (
                len(trange) == len(gt)
                and len(trange) == len(ht)
                and len(trange) == len(rotS)
            ):
                meta.success = False
                meta.error = "xtb_trange_inconsistent"
                return result, meta
            else:
                result.gibbs = gt
                result.enthalpy = ht
                result.entropy = rotS

        # Extract symmetry
        result.linear = next(
            (
                {"true": True, "false": False}[line.split()[2]]
                for line in lines
                if ":  linear? " in line
            )
        )

        # Extract rmsd
        if config.bhess:
            result.rmsd = next(
                (float(line.split()[3]) for line in lines if "final rmsd / " in line)
            )

        # xtb_enso.json is generated by xtb by using the '--enso' argument *only* when using --bhess or --ohess
        # (when a hessian is calculated)
        # contains output from xtb in json format to be more easily digestible by CENSO
        with open(
            os.path.join(jobdir, "xtb_enso.json"),
            "r",
        ) as f:
            data = json.load(f)

        # read number of imaginary frequencies and print warning
        if "number of imags" in data:
            if data["number of imags"] > 0:
                logger.warning(
                    f"Found {data['number of imags']} significant"
                    + f" imaginary frequencies for "
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
                if not config.multitemp:
                    result.entropy[temp] = (
                        0.0  # set this to None for predictability (?)
                    )
            else:
                result.energy = data.get("G(T)", 0.0)
                result.gibbs[temp] = data.get("G(T)", 0.0)
                if not config.multitemp:
                    result.enthalpy[temp] = (
                        0.0  # set this to None for predictability (?)
                    )
                    result.entropy[temp] = (
                        0.0  # set this to None for predictability (?)
                    )

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

    @abstractmethod
    @_run
    def sp(
        self, job: ParallelJob, jobdir: Path | str, config: SPJobConfig, **kwargs
    ) -> tuple[SPResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @_run
    def gsolv(
        self, job: ParallelJob, jobdir: Path | str, config: SPJobConfig, **kwargs
    ) -> tuple[GsolvResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @_run
    def xtb_opt(
        self, job: ParallelJob, jobdir: Path | str, config: XTBOptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @_run
    def opt(
        self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs
    ) -> tuple[OptResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @_run
    def nmr(
        self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs
    ) -> tuple[NMRResult, MetaData]:
        raise NotImplementedError

    @abstractmethod
    @_run
    def uvvis(
        self, job: ParallelJob, jobdir: Path | str, config: OptJobConfig, **kwargs
    ) -> tuple[UVVisResult, MetaData]:
        raise NotImplementedError


Factory.register_builder("xtb", QmProc)
