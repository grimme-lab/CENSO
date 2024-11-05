"""
Contains QmProc base class,
Additionally contains functions which should be present irrespective of the QM
code. (xTB always available)
"""

import json
import os
import signal
import subprocess
from time import perf_counter
from collections.abc import Callable

from .datastructure import ParallelJob
from .params import (
    Config,
    PLENGTH,
    DIGILEN,
    WARNLEN,
)
from .utilities import print, frange, Factory
from .logging import setup_logger

logger = setup_logger(__name__)


def handle_sigterm(signum, frame, sub):
    logger.critical(
        f"{f'worker{os.getpid()}:':{WARNLEN}}Received SIGTERM. Terminating."
    )
    sub.send_signal(signal.SIGTERM)


class QmProc:
    """
    QmProc base class
    """

    _progname = "generic"

    _paths = {
        "orcapath": "",
        "orcaversion": "",
        "xtbpath": "",
        "crestpath": "",
        "cefinepath": "",
        "cosmorssetup": "",
        "cosmothermpath": "",
        "dbpath": "",
        "cosmothermversion": "",
        "mpshiftpath": "",
        "escfpath": "",
    }

    _req_settings_xtb = {
        "xtb_sp": [
            "gfnv",
        ],
        "xtb_gsolv": [
            # no special requirements
        ],
        "xtb_rrho": [
            "gfnv",
        ],
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
        for program, path in cls._paths.items():
            # Append a line with the program and its path to the output.
            lines.append(f"{program}:".ljust(DIGILEN, " ") + f"{path}\n")

        # Print each line of the output.
        for line in lines:
            print(line)

    def __init__(self, workdir: str):
        # dict to map the jobtypes to their respective methods
        self._jobtypes: dict[str, Callable] = {
            "xtb_sp": self._xtb_sp,
            "xtb_gsolv": self._xtb_gsolv,
            "xtb_rrho": self._xtb_rrho,
        }

        self.workdir = workdir

    def run(self, job: ParallelJob) -> ParallelJob:
        """
        Run methods depending on jobtype.
        DO NOT OVERRIDE OR OVERLOAD! this will break e.g. censo.parallel.execute

        Args:
            job (ParallelJob): job to run

        Returns:
            job (ParallelJob): job with results
        """
        logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running on {job.omp} cores.")
        # jobtype is basically an ordered (!!!) (important e.g. if sp is required before the next step)
        # list containing the types of computations to do
        if not all(t in self._jobtypes for t in job.jobtype):
            raise RuntimeError(
                f"At least one jobtype of {job.jobtype} is not available for "
                + f"{self.__class__.__name__}.\nAvailable "
                + f"jobtypes are: {list(self._jobtypes)}"
            )

        # run all the computations
        for j in job.jobtype:
            # Create jobdir
            jobdir = self._create_jobdir(job.conf.name, j)

            # Time execution
            start = perf_counter()

            logger.info(
                f"{f'worker{os.getpid()}:':{WARNLEN}}Running "
                + f"{j} calculation in {jobdir}."
            )
            print(f"Running {j} calculation for {job.conf.name}.")
            job.results[j], job.meta[j] = self._jobtypes[j](job, jobdir)

            # Copy mo path if possible to be used for further calculations
            # (processors grab the mo path from the 'mo_path' key that is always present in the meta dict)
            if job.meta[j].get("mo_path", None) is not None:
                job.meta["mo_path"] = job.meta[j]["mo_path"]

            end = perf_counter()

            job.meta[j]["time"] = end - start

            # if a calculation failed all following calculations will not be executed
            if not job.meta[j]["success"]:
                for j2 in job.jobtype[job.jobtype.index(j) + 1 :]:
                    job.results[j2] = None
                    job.meta[j2]["success"] = False
                    job.meta[j2]["error"] = "Previous calculation failed"
                break

        job.meta["total_time"] = sum(
            job.meta[j]["time"]
            for j in job.jobtype
            if job.meta[j]["error"] != "Previous calculation failed"
        )

        # returns modified job object with result dict e.g.: {"sp": ..., "gsolv": ..., etc.}
        return job

    def _make_call(
        self, prog: str, call: list, outputpath: str, jobdir: str
    ) -> tuple[int, str]:
        """
        Make a call to an external program and write output into outputfile.

        Args:
            prog (str): program to call
            call (list): list containing the call args to the external program
            outputpath (str): path to the outputfile
            jobdir (str): path to the jobdir

        Returns:
            returncode (int): returncode of the external program
        """
        # make sure program path is not empty
        if prog != "tm":
            pathmap = {
                "xtb": "xtbpath",
                "orca": "orcapath",
            }
            try:
                assert self._paths[pathmap[prog]].strip() != ""
                # NOTE: turbomole does not need this step since you can just call a binary w/o path
                # Also the binaries are called differently for different purposes
                # (e.g. ridft for single-points, but not for NMR, instead the binary names are passed in the call)
                call.insert(0, self._paths[pathmap[prog]])
            except AssertionError as exc:
                raise AssertionError(
                    f"Path for {prog} not found. Please set up {pathmap[prog]} in the rcfile."
                ) from exc

        # call external program and write output into outputfile
        with open(outputpath, "w", newline=None) as outputfile:
            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running {call}...")

            # create subprocess for external program
            sub = subprocess.Popen(
                call,
                shell=False,
                text=True,
                stdin=None,
                stderr=subprocess.PIPE,
                cwd=jobdir,
                stdout=outputfile,
                env=Config.ENVIRON,
            )

            logger.debug(
                f"{f'worker{os.getpid()}:':{WARNLEN}}Started (PID: {sub.pid})."
            )

            # make sure to send SIGTERM to subprocess if program is quit
            signal.signal(
                signal.SIGTERM, lambda signum, frame: handle_sigterm(signum, frame, sub)
            )

            # wait for process to finish
            _, errors = sub.communicate()
            returncode = sub.returncode

            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Done.")

        return returncode, errors

    def _create_jobdir(self, confname: str, job: str) -> str:
        """
        Creates a subdir in confdir for the job.
        """

        jobdir = os.path.join(self.workdir, confname, job)
        try:
            # Create the directory
            os.makedirs(jobdir)
        except FileExistsError:
            # logger.warning(
            #    f"{f'worker{os.getpid()}:':{WARNLEN}}Jobdir {jobdir} already exists!"
            #    " Files will be overwritten."
            # )
            pass

        return jobdir

    def _get_sym_num(self, sym=None, linear=None):
        """Get rotational symmetry number from Schoenfließ symbol"""
        if sym is None:
            sym = "c1"
        if linear is None:
            linear = False
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

    def _xtb_sp(
        self,
        job: ParallelJob,
        jobdir: str,
        filename: str = "xtb_sp",
        no_solv: bool = False,
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        Calculates the single-point energy with GFNn-xTB or GFN-FF.

        Args:
            job (ParallelJob): job to run
            jobdir (str): path to the jobdir
            filename (str, optional): filename to use for the coord file. Defaults to "xtb_sp".
            no_solv (bool, optional): whether to run the sp in gas-phase. Defaults to False.

        Returns:
            result (dict[str, float | None]): result of the sp calculation
            meta (dict[str, any]): metadata about the job

        result = {
            "energy": None,
        }
        """
        # set results
        result = {
            "energy": None,
        }

        # set metadata
        meta = {
            "success": None,
            "error": None,
        }

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
        call = [
            f"{filename}.coord",
            "--" + job.prepinfo["xtb_sp"]["gfnv"],
            "--sp",
            "--chrg",
            f"{job.prepinfo['charge']}",
            "--norestart",
            "--parallel",
            f"{job.omp}",
        ]

        # add solvent to xtb call if not a gas-phase sp
        # (set either through run settings or by call kwarg e.g. for _xtb_gsolv)
        # NOTE on solvents_dict (or rather censo_solvents.json):
        # [0] is the normal name of the solvent, if it is available, [1] is the replacement
        if not (job.prepinfo["general"].get("gas-phase", False) or no_solv):
            call.extend(
                [
                    "--" + job.prepinfo["general"]["sm_rrho"],
                    job.prepinfo["xtb_sp"]["solvent_key_xtb"],
                    "reference",
                    "-I",
                    xcontrolname,
                ]
            )

            # set gbsa grid
            with open(os.path.join(jobdir, xcontrolname), "w", newline=None) as xcout:
                xcout.write("$gbsa\n")
                xcout.write("  gbsagrid=tight\n")
                xcout.write("$end\n")

        # call xtb
        returncode, errors = self._make_call("xtb", call, outputpath, jobdir)

        # if returncode != 0 then some error happened in xtb
        # TODO - returncodes
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "unknown_error"
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            return result, meta

        # read energy from outputfile
        with open(outputpath, "r", encoding=Config.CODING, newline=None) as outputfile:
            for line in outputfile.readlines():
                if "| TOTAL ENERGY" in line:
                    result["energy"] = float(line.split()[3])
                    meta["success"] = True
                # TODO - important - what to do if calculation not converged?

        # FIXME - right now the case meta["success"] = None might appear if "TOTAL ENERGY" is not found in outputfile
        return result, meta

    def _xtb_gsolv(
        self, job: ParallelJob, jobdir: str
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        Calculate additive GBSA or ALPB solvation using GFNn-xTB or GFN-FF.

        Args:
            job (ParallelJob): job to run
            jobdir (str): path to the jobdir

        Returns:
            result (dict[str, any]): result of the gsolv calculation

        result = {
            "gsolv": None,
            "energy_xtb_gas": None,
            "energy_xtb_solv": None,
        }
        """
        # what is returned in the end
        result = {
            "gsolv": None,
            "energy_xtb_gas": None,
            "energy_xtb_solv": None,
        }

        meta = {
            "success": None,
            "error": None,
        }

        # run gas-phase GFN single-point
        spres, spmeta = self._xtb_sp(job, jobdir, filename="gas", no_solv=True)
        if spmeta["success"]:
            result["energy_xtb_gas"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = spmeta["error"]
            return result, meta

        # run single-point in solution:
        # ''reference'' corresponds to 1\;bar of ideal gas and 1\;mol/L of liquid
        #   solution at infinite dilution,
        spres, spmeta = self._xtb_sp(job, jobdir, filename="solv")
        if spmeta["success"]:
            result["energy_xtb_solv"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = spmeta["error"]
            return result, meta

        # only reached if both gas-phase and solvated sp succeeded
        result["gsolv"] = result["energy_xtb_solv"] - result["energy_xtb_gas"]
        meta["success"] = True

        return result, meta

    # TODO - break this down
    def _xtb_rrho(
        self, job: ParallelJob, jobdir: str, filename: str = "xtb_rrho"
    ) -> tuple[dict[str, any], dict[str, any]]:
        """
        Calculates the mRRHO contribution to the free enthalpy of a conformer with GFNn-xTB/GFN-FF.

        Args:
            job (ParallelJob): job to run
            jobdir (str): path to the jobdir
            filename (str, optional): filename to use for the coord file. Defaults to "xtb_rrho".

        Returns:
            result (dict[str, any]): result of the rrho calculation

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
        result = {
            "energy": None,
            "rmsd": None,
            "gibbs": None,
            "enthalpy": None,
            "entropy": None,
            "symmetry": None,
            "symnum": None,
        }

        meta = {
            "success": None,
            "error": None,
        }

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
            if job.prepinfo["general"]["multitemp"]:
                trange = frange(
                    job.prepinfo["general"]["trange"][0],
                    job.prepinfo["general"]["trange"][1],
                    step=job.prepinfo["general"]["trange"][2],
                )

                # Always append the fixed temperature to the trange so that it is the last value
                # (important since --enso will make xtb give the G(T) value for this temperature)
                trange.append(job.prepinfo["general"]["temperature"])

                # Write trange to the xcontrol file
                xcout.write(f"    temp=" f"{','.join([str(i) for i in trange])}\n")
            else:
                xcout.write(f"    temp={job.prepinfo['general']['temperature']}\n")

            xcout.write(f"    sthr={job.prepinfo['general']['sthr']}\n")

            xcout.write(f"    imagthr={job.prepinfo['general']['imagthr']}\n")

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
            if not job.prepinfo["general"]["gas-phase"]:
                xcout.write("$gbsa\n")
                xcout.write("  gbsagrid=tight\n")

            xcout.write("$end\n")

        if job.prepinfo["general"]["bhess"]:
            # set ohess or bhess
            dohess = "--bhess"
        else:
            dohess = "--ohess"

        # generate coord file for xtb
        with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as file:
            file.writelines(job.conf.tocoord())

        call = [
            f"{filename}.coord",
            "--" + job.prepinfo["xtb_rrho"]["gfnv"],
            dohess,
            "vtight",
            "--chrg",
            f"{job.prepinfo['charge']}",
            "--enso",
            "--norestart",
            "-I",
            xcontrolname,
            "--parallel",
            f"{job.omp}",
        ]

        # add solvent to xtb call if necessary
        if not job.prepinfo["general"]["gas-phase"]:
            call.extend(
                [
                    "--" + job.prepinfo["general"]["sm_rrho"],
                    job.prepinfo["xtb_rrho"]["solvent_key_xtb"],
                ]
            )

        # if rmsd bias is used (should be defined in censo workdir (cwd)) TODO
        if job.prepinfo["general"]["rmsdbias"]:
            # move one dir up to get to cwd (FIXME)
            cwd = os.path.join(os.path.split(self.workdir)[::-1][1:][::-1])

            call.extend(
                [
                    "--bias-input",
                    os.path.join(cwd, "rmsdpot.xyz"),
                ]
            )

        # call xtb
        returncode, errors = self._make_call("xtb", call, outputpath, jobdir)

        # check if converged:
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "unknown_error"
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            return result, meta

        # read output and store lines
        with open(outputpath, "r", encoding=Config.CODING, newline=None) as outputfile:
            lines = outputfile.readlines()

        if job.prepinfo["general"]["multitemp"]:
            # get gibbs energy, enthalpy and entropy for given temperature range
            trange = frange(
                job.prepinfo["general"]["trange"][0],
                job.prepinfo["general"]["trange"][1],
                step=job.prepinfo["general"]["trange"][2],
            )

            # gibbs energy
            gt = {}

            # enthalpy
            ht = {}

            # Get Gibbs energy and enthalpy
            for line in lines:
                if "T/K" in line:
                    for line2 in lines[lines.index(line) + 2 :]:
                        if "----------------------------------" in line2:
                            break

                        T = float(line2.split()[0])
                        gt[T] = float(line2.split()[4])
                        ht[T] = float(line2.split()[2])

            # rotational entropy
            rotS = {}

            # Get rotational entropy
            entropy_lines = (
                (line, lines[i + 1]) for i, line in enumerate(lines) if "VIB" in line
            )
            for line in entropy_lines:
                T = float(line[0].split()[0])
                rotS[T] = float(line[1].split()[4])

        # Extract symmetry
        result["linear"] = next(
            (
                {"true": True, "false": False}[line.split()[2]]
                for line in lines
                if ":  linear? " in line
            ),
            None,
        )

        # Extract rmsd
        result["rmsd"] = next(
            (
                float(line.split()[3])
                for line in lines
                if "final rmsd / " in line and job.prepinfo["general"]["bhess"]
            ),
            None,
        )

        # check if xtb calculated the temperature range correctly
        if job.prepinfo["general"]["multitemp"] and not (
            len(trange) == len(gt)
            and len(trange) == len(ht)
            and len(trange) == len(rotS)
        ):
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_rrho"
            return result, meta
        elif job.prepinfo["general"]["multitemp"]:
            result["gibbs"] = gt
            result["enthalpy"] = ht
            result["entropy"] = rotS
        else:
            result["gibbs"] = {}
            result["enthalpy"] = {}
            result["entropy"] = {}

        # xtb_enso.json is generated by xtb by using the '--enso' argument *only* when using --bhess or --ohess
        # (when a hessian is calculated)
        # contains output from xtb in json format to be more easily digestible by CENSO
        with open(
            os.path.join(jobdir, "xtb_enso.json"),
            "r",
            encoding=Config.CODING,
            newline=None,
        ) as f:
            data = json.load(f)

        # read number of imaginary frequencies and print warning
        if "number of imags" in data:
            if data["number of imags"] > 0:
                logger.warning(
                    f"Found {data['number of imags']} significant"
                    f" imaginary frequencies for "
                    f"{job.conf.name}."
                )

        # get gibbs energy
        if "G(T)" in data:
            meta["success"] = True

            if job.prepinfo["general"]["temperature"] == 0.0:
                result["energy"] = data.get("ZPVE", 0.0)
                result["gibbs"][job.prepinfo["general"]["temperature"]] = data.get(
                    "ZPVE", 0.0
                )
                result["enthalpy"][job.prepinfo["general"]["temperature"]] = data.get(
                    "ZPVE", 0.0
                )
                if not job.prepinfo["general"]["multitemp"]:
                    result["entropy"][
                        job.prepinfo["general"]["temperature"]
                    ] = None  # set this to None for predictability
            else:
                result["energy"] = data.get("G(T)", 0.0)
                result["gibbs"][job.prepinfo["general"]["temperature"]] = data.get(
                    "G(T)", 0.0
                )
                if not job.prepinfo["general"]["multitemp"]:
                    result["enthalpy"][
                        job.prepinfo["general"]["temperature"]
                    ] = None  # set this to None for predictability
                    result["entropy"][
                        job.prepinfo["general"]["temperature"]
                    ] = None  # set this to None for predictability
                    # FIXME - why though?

            # only determine symmetry if all the needed information is there
            if "point group" and "linear" in data.keys():
                result["symmetry"] = data["point group"]
                result["linear"] = data.get("linear", result["linear"])
            else:
                # could not determine symmetry correctly
                result["symmetry"] = "c1"
                result["linear"] = False

            # calculate symnum
            result["symnum"] = self._get_sym_num(
                sym=result["symmetry"], linear=result["linear"]
            )
        else:
            meta["success"] = False
            meta["error"] = "Could not read xtb_enso.json"

        return result, meta


Factory.register_builder("xtb", QmProc)
