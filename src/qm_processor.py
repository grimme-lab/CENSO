"""
Contains QmProc base class, 
Additionally contains functions which should be present irrespective of the QM
code. (xTB always available)
"""
import functools
import json
import os
import signal
import subprocess
from time import perf_counter
from typing import Any, Callable, Dict, List

from src.datastructure import ParallelJob
from src.params import (
    ENVIRON,
    CODING,
    rot_sym_num,
    PLENGTH, DIGILEN, WARNLEN,
)
from src.utilities import print, frange, setup_logger

logger = setup_logger(__name__)


def handle_sigterm(signum, frame, sub):
    global logger
    logger.critical(f"{f'worker{os.getpid()}:':{WARNLEN}}Received SIGTERM. Terminating.")
    sub.send_signal(signal.SIGTERM)


class QmProc:
    """
    QmProc base class
    """

    _paths = {
        "orcapath": "",
        "orcaversion": "",
        "xtbpath": "",
        "crestpath": "",
        "cosmorssetup": "",
        "dbpath": "",
        "cosmothermversion": "",
        "mpshiftpath": "",
        "escfpath": "",
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

    # FIXME - for gsolv (calls the sp methods) this still tries to create a sp subdir, even though this is not correct
    @staticmethod
    def _create_jobdir(f: Callable) -> Callable:
        """
        Creates a subdir in confdir for the job.
        
        This method needs to be defined as @staticmethod to be accessible from within the class via the @_create_jobdir decorator.
        The wrapper function within will be able to access the instance variables of the class.
        To access this method from child classes, the decorator must be called like: @QmProc._create_jobdir.
        """

        @functools.wraps(f)
        def wrapper(self, job, *args, **kwargs):
            """
            A function that wraps the given `f` function and performs some operations before and after calling it.

            Parameters:
                self (object): The instance of the class that the function is a method of.
                job (object): The job object that contains the necessary information.
                *args (tuple): The positional arguments passed to the `f` function.
                **kwargs (dict): The keyword arguments passed to the `f` function.

            Returns:
                The return value of the `f` function.
            """
            # getting the name starting from 1: to strip the _
            jobdir = os.path.join(self.workdir, job.conf.name, f.__name__[1:])
            try:
                # Create the directory
                os.makedirs(jobdir)
            except FileExistsError:
                global logger
                logger.warning(f"{f'worker{os.getpid()}:':{WARNLEN}}Jobdir {jobdir} already exists!"
                               " Files will be overwritten.")

            return f(self, job, *args, **kwargs)

        return wrapper

    def __init__(self, instructions: Dict[str, Any], workdir: str):
        # stores instructions, meaning the settings that should be applied for all jobs
        # e.g. 'gfnv' (GFN version for xtb_sp/xtb_rrho/xtb_gsolv)
        self.instructions: Dict[str, Any] = instructions

        # dict to map the jobtypes to their respective methods
        self._jobtypes: Dict[str, Callable] = {
            "sp": self._sp,
            "gsolv": self._gsolv,
            "opt": self._opt,
            "genericoutput": self._genericoutput,
            "xtb_sp": self._xtb_sp,
            "xtb_gsolv": self._xtb_gsolv,
            "xtb_opt": self._xtb_opt,
            "xtb_rrho": self._xtb_rrho,
        }

        self.workdir = workdir

    def run(self, job: ParallelJob) -> ParallelJob:
        """
        run methods depending on jobtype
        DO NOT OVERRIDE OR OVERLOAD! this will break e.g. ProcessHandler.execute
        """
        global logger
        logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running on {job.omp} cores.")
        # jobtype is basically an ordered (!!!) (important e.g. if sp is required before the next step)
        # list containing the types of computations to do
        if not all(t in self._jobtypes.keys() for t in job.jobtype):
            raise RuntimeError(
                f"At least one jobtype of {self._jobtypes} is not available for {self.__class__.__name__}.\nAvailable "
                + f"jobtypes are: {self._jobtypes.keys()}")

        # run all the computations
        for j in job.jobtype:
            # Time execution
            start = perf_counter()

            job.results[j] = self._jobtypes[j](job)

            end = perf_counter()

            job.meta[j]["time"] = end - start

            # if a calculation failed all following calculations will not be executed
            if not job.meta[j]["success"]:
                for j2 in job.jobtype[job.jobtype.index(j) + 1:]:
                    job.meta[j2]["success"] = False
                    job.meta[j2]["error"] = "Previous calculation failed"
                break

        job.meta["total_time"] = sum(
            job.meta[j]["time"] for j in job.jobtype
            if job.meta[j]["error"] != "Previous calculation failed"
        )

        # returns modified job object with result dict e.g.: {"sp": ..., "gsolv": ..., etc.}
        return job

    @staticmethod
    def _make_call(call: List, outputpath: str, jobdir: str) -> int:
        # call external program and write output into outputfile
        with open(outputpath, "w", newline=None) as outputfile:
            global logger
            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running {call}...")

            # create subprocess for external program
            sub = subprocess.Popen(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=jobdir,
                stdout=outputfile,
                env=ENVIRON,
            )

            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Started (PID: {sub.pid}).")

            # make sure to send SIGTERM to subprocess if program is quit
            signal.signal(signal.SIGTERM, lambda signum, frame: handle_sigterm(signum, frame, sub))

            # wait for process to finish
            returncode = sub.wait()

            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Done.")

        return returncode

    def _sp(self):
        """
        single-point calculation
        """
        pass

    def _opt(self):
        """
        geometry ensembleopt
        """
        pass

    def _gsolv(self):
        """
        gsolv calculation using the respective solvent model
        """
        pass

    def _genericoutput(self):
        """
        Read shielding and coupling constants and write them to plain output
        The first natom lines contain the shielding constants, and from
        line natom +1 the coupling constants are written.
        """
        pass

    @staticmethod
    def _get_sym_num(sym=None, linear=None):
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
        for key in rot_sym_num.keys():
            if key in sym.lower():
                symnum = rot_sym_num.get(key, 1)
                break
        return symnum

    @_create_jobdir
    def _xtb_sp(self, job: ParallelJob, filename: str = "xtb_sp", no_solv: bool = False, silent: bool = False,
                jobtype: str = "xtb_sp") -> dict[str, float | None]:
        """
        Get single-point energy from xtb
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
        jobdir = os.path.join(self.workdir, job.conf.name, jobtype)
        inputpath = os.path.join(jobdir, f"{filename}.coord")
        outputpath = os.path.join(jobdir, f"{filename}.out")
        xcontrolname = "xtb_sp-xcontrol-inp"

        global logger
        if not silent:
            logger.info(f"{f'worker{os.getpid()}:':{WARNLEN}}Running xtb_sp calculation in {jobdir}.")
        else:
            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running xtb_sp calculation in {jobdir}.")

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
            self._paths["xtbpath"],
            f"{filename}.coord",
            "--" + self.instructions["gfnv"],
            "--sp",
            "--chrg",
            f"{self.instructions['charge']}",
            "--norestart",
            "--parallel",
            f"{job.omp}",
        ]

        # add solvent to xtb call if not a gas-phase sp 
        # (set either through run settings or by call kwarg e.g. for _xtb_gsolv)
        # NOTE on solvents_dict (or rather censo_solvents.json): 
        # [0] is the normal name of the solvent, if it is available, [1] is the replacement
        if not (self.instructions.get("gas-phase", False) or no_solv):
            call.extend(
                [
                    "--" + self.instructions["sm_rrho"],
                    self.instructions["solvent_key_xtb"],
                    "reference",
                    "-I",
                    xcontrolname
                ]
            )

            # set gbsa grid
            with open(
                    os.path.join(jobdir, xcontrolname), "w", newline=None
            ) as xcout:
                xcout.write("$gbsa\n")
                xcout.write("  gbsagrid=tight\n")
                xcout.write("$end\n")

        # call xtb
        returncode = self._make_call(call, outputpath, jobdir)

        # if returncode != 0 then some error happened in xtb
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_sp"
            job.meta[jobtype].update(meta)
            return result

        # read energy from outputfile
        with open(outputpath, "r", encoding=CODING, newline=None) as outputfile:
            for line in outputfile.readlines():
                if "| TOTAL ENERGY" in line:
                    result["energy"] = float(line.split()[3])
                    meta["success"] = True
                # TODO - what to do if calculation not converged?

        # FIXME - right now the case meta["success"] = None might appear if "TOTAL ENERGY" is not found in outputfile
        job.meta[jobtype].update(meta)
        return result

    @_create_jobdir
    def _xtb_gsolv(self, job: ParallelJob) -> dict[str, Any | None]:
        """
        Calculate additive GBSA or ALPB solvation contribution by
        Gsolv = Esolv - Egas, using GFNn-xTB or GFN-FF
        result = {
            "gsolv": None,
            "energy_xtb_gas": None,
            "energy_xtb_solv": None,
        }
        """
        jobdir = os.path.join(self.workdir, job.conf.name, "xtb_gsolv")

        global logger
        logger.info(
            f"{f'worker{os.getpid()}:':{WARNLEN}}Running xtb_gsolv calculation in {jobdir}."
        )

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
        # TODO - does this need it's own folder?
        spres = self._xtb_sp(job, filename="gas", silent=True, no_solv=True, jobtype="xtb_gsolv")
        if job.meta["xtb_gsolv"]["success"]:
            result["energy_xtb_gas"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_gsolv"
            job.meta["xtb_gsolv"].update(meta)
            return result

        # run single-point in solution:
        # ''reference'' corresponds to 1\;bar of ideal gas and 1\;mol/L of liquid
        #   solution at infinite dilution,
        spres = self._xtb_sp(job, filename="solv", silent=True, jobtype="xtb_gsolv")
        if job.meta["xtb_gsolv"]["success"]:
            result["energy_xtb_solv"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_gsolv"
            job.meta["xtb_gsolv"].update(meta)
            return result

        # only reached if both gas-phase and solvated sp succeeded   
        result["gsolv"] = result["energy_xtb_solv"] - result["energy_xtb_gas"]
        meta["success"] = True

        job.meta["xtb_gsolv"].update(meta)
        return result

    def _xtb_opt(self):
        """
        geometry ensembleopt using xtb as driver, has to be implemented for each qm code
        """
        pass

    # TODO - break this down
    @_create_jobdir
    def _xtb_rrho(self, job: ParallelJob, filename: str = "xtb_rrho") -> dict[str, Any]:
        """
        mRRHO contribution with GFNn-xTB/GFN-FF
        
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

        # if not self.job["onlyread"]: # TODO ???

        # set in/out path
        jobdir = os.path.join(self.workdir, job.conf.name, "xtb_rrho")
        outputpath = os.path.join(jobdir, f"{filename}.out")
        xcontrolname = "rrho-xcontrol-inp"
        xcontrolpath = os.path.join(jobdir, xcontrolname)

        global logger
        logger.info(
            f"{f'worker{os.getpid()}:':{WARNLEN}}Running {str(self.instructions['gfnv']).upper()} mRRHO in {jobdir}."
        )

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

        # if a temperature in the trange is close (+-0.6) to the fixed temperature then replace this value in trange
        # don't do this for now, to make the program more predictable
        """ for t in trange:
            if isclose(self.instructions["temperature"], t, abs_tol=0.6):
                trange[trange.index(t)] = self.instructions["temperature"]
                break """

        # setup xcontrol
        with open(xcontrolpath, "w", newline=None) as xcout:
            xcout.write("$thermo\n")
            if self.instructions.get("trange", False):
                trange = frange(self.instructions['trange'][0], self.instructions['trange'][1],
                                self.instructions['trange'][2])
                xcout.write(
                    f"    temp="
                    f"{','.join([str(i) for i in trange])}\n")
            else:
                xcout.write(f"    temp={self.instructions['temperature']}\n")

            xcout.write(f"    sthr={self.instructions['sthr']}\n")

            xcout.write(f"    imagthr={self.instructions['imagthr']}\n")

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
            if not self.instructions["gas-phase"]:
                xcout.write("$gbsa\n")
                xcout.write("  gbsagrid=tight\n")

            xcout.write("$end\n")

        if self.instructions["bhess"]:
            # set ohess or bhess
            dohess = "--bhess"
            olevel = "vtight"
        else:
            dohess = "--ohess"
            olevel = "vtight"

        # generate coord file for xtb
        with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as file:
            file.writelines(job.conf.tocoord())

        call = [
            self._paths["xtbpath"],
            f"{filename}.coord",
            "--" + self.instructions["gfnv"],
            dohess,
            olevel,
            "--chrg",
            f"{self.instructions['charge']}",
            "--enso",
            "--norestart",
            "-I",
            xcontrolname,
            "--parallel",
            f"{job.omp}",
        ]

        # add solvent to xtb call if necessary
        if not self.instructions["gas-phase"]:
            call.extend(
                [
                    "--" + self.instructions["sm_rrho"],
                    self.instructions["solvent_key_xtb"],
                ]
            )

        # if rmsd bias is used (should be defined in censo workdir (cwd)) TODO
        if self.instructions["rmsdbias"]:
            # move one dir up to get to cwd (FIXME)
            cwd = os.path.join(os.path.split(self.workdir)[::-1][1:][::-1])

            call.extend(
                [
                    "--bias-input",
                    os.path.join(cwd, "rmsdpot.xyz"),
                ]
            )

        # call xtb
        returncode = self._make_call(call, outputpath, jobdir)

        # check if converged:
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_rrho"
            job.meta["xtb_rrho"].update(meta)
            return result

        # read output and store lines
        with open(outputpath, "r", encoding=CODING, newline=None) as outputfile:
            lines = outputfile.readlines()

        # get gibbs energy, enthalpy and entropy for given temperature range
        trange = frange(self.instructions["trange"][0], self.instructions["trange"][1], self.instructions["trange"][2])

        # gibbs energy
        gt = {}

        # enthalpy
        ht = {}

        # rotational entropy
        rotS = {}

        # Extract symmetry
        result["linear"] = next(
            ({"true": True, "false": False}[line.split()[2]] for line in lines if ":  linear? " in line), None)

        # Extract rmsd
        result["rmsd"] = next(
            (float(line.split()[3]) for line in lines if "final rmsd / " in line and self.instructions["bhess"]), None)

        # Get rotational entropy
        entropy_lines = ((line, lines[i + 1]) for i, line in enumerate(lines) if "VIB" in line)
        for line in entropy_lines:
            T = float(line[0].split()[0])
            rotS[T] = float(line[1].split()[4])

        # Get Gibbs energy and enthalpy
        for line in lines:
            if "T/K" in line:
                for line2 in lines[lines.index(line) + 2:]:
                    if "----------------------------------" in line2:
                        break

                    T = float(line2.split()[0])
                    gt[T] = float(line2.split()[4])
                    ht[T] = float(line2.split()[2])

        # check if xtb calculated the temperature range correctly
        if len(trange) == len(gt) and len(trange) == len(ht) and len(trange) == len(rotS):
            result["gibbs"] = gt
            result["enthalpy"] = ht
            result["entropy"] = rotS
        else:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_rrho"
            job.meta["xtb_rrho"].update(meta)
            return result

        # xtb_enso.json is generated by xtb by using the '--enso' argument *only* when using --bhess or --ohess (when a hessian is calculated)
        # contains output from xtb in json format to be more easily digestible by CENSO
        with open(
                os.path.join(jobdir, "xtb_enso.json"),
                "r",
                encoding=CODING,
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

            if self.instructions["temperature"] == 0:
                result["energy"] = data.get("ZPVE", 0.0)
                result["gibbs"][self.instructions["temperature"]] = data.get("ZPVE", 0.0)
                result["enthalpy"][self.instructions["temperature"]] = data.get("ZPVE", 0.0)
                result["entropy"][self.instructions["temperature"]] = None  # set this to None for predictability
            else:
                result["energy"] = data.get("G(T)", 0.0)
                result["gibbs"][self.instructions["temperature"]] = data.get("G(T)", 0.0)
                result["enthalpy"][self.instructions["temperature"]] = None  # set this to None for predictability
                result["entropy"][self.instructions["temperature"]] = None  # set this to None for predictability

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

        job.meta["xtb_rrho"].update(meta)
        return result
