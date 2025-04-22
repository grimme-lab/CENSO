"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""

import os
import shutil
from collections import OrderedDict
from functools import reduce
from typing import Any

from .utilities import od_insert, Factory
from .logging import setup_logger
from .datastructure import GeometryData, ParallelJob
from .params import (
    Config,
    WARNLEN,
)
from .qm_processor import QmProc

logger = setup_logger(__name__)


class OrcaProc(QmProc):
    """
    Performs calculations with ORCA.
    """

    _progname = "orca"

    # contains grid settings for ORCA 5.0+ (True) and older versions (False)
    # can be chosen by simple keyword (low/low+/high/high+)
    __gridsettings = {
        False: {
            "low": ["grid4", "nofinalgrid", "loosescf"],
            "low+": ["grid4", "nofinalgrid", "scfconv6"],
            "high": ["grid4", "nofinalgrid", "scfconv7"],
            "high+": ["grid5", "nofinalgrid", "scfconv7"],
            "nmr": ["grid5", "nofinalgrid", "scfconv7"],
        },
        True: {
            "low": ["DEFGRID1", "loosescf"],
            "low+": ["DEFGRID2", "scfconv6"],
            "high": ["DEFGRID2", "scfconv7"],
            "high+": ["DEFGRID2", "scfconv7"],
            "nmr": ["DEFGRID2", "scfconv7"],
        },
    }

    # Dict to map orca returncodes to error messages
    __returncode_to_err: dict[int, str] = {
        24: "input_error",
        25: "input_error",
        25: "input_error",
        30: "input_error",
        52: "input_error",
        55: "input_error",
        125: "unknown_error",
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # expand jobtypes with special orca jobtypes
        self._jobtypes = {
            **self._jobtypes,
            **{
                "sp": self._sp,
                "gsolv": self._gsolv,
                "xtb_opt": self._xtb_opt,
                "opt": self._opt,
                "nmr": self._nmr,
                "uvvis": self._uvvis,
            },
        }

        # Stores setting wether to copy MO-files for faster SCFs
        self.copy_mo: bool = False

    def __prep(
        self, job: ParallelJob, jobtype: str, no_solv: bool = False, xyzfile: str = None
    ) -> list[str]:
        """
        Prepares an OrderedDict to be fed into the parser in order to write an input file for jobtype 'jobtype'
        (e.g. sp).

        Can load up a template file from user assets folder.

        Args:
            job: ParallelJob object containing the job information
            jobtype: jobtype to prepare the input for
            no_solv: if True, no solvent model is used
            xyzfile: if not None, the geometry is read from this file instead of the job object

        Returns:
            list[str]: the orca input.

        NOTE: the xyzfile has to already exist, this function just bluntly writes the name into the input.
        """

        # check ORCA version (orca5 = True means at least ORCA version 5)
        orca5 = not self._paths["orcaversion"].startswith("4")

        inp = []

        if job.prepinfo[jobtype]["template"]:
            # NOTE: when using templates we're not going to check for double definitions!
            # if the template is messed up, orca will fail and the user should deal with that
            # load template file
            with open(
                os.path.join(
                    Config.USER_ASSETS_PATH,
                    f"{job.prepinfo['partname']}.orca.template",
                ),
                "r",
            ) as f:
                inp = f.readlines()

            main_line = next(inp.index(l) for l in inp if "{main}" in l)
            inp.pop(main_line)

        # prepare the main line of the orca input
        inp[:0] = self.__prep_main(job.prepinfo, jobtype, orca5)

        # prepare all options that are supposed to be placed before the
        # geometry definition
        inp[1:1] = self.__prep_pregeom(job.prepinfo, orca5, jobtype, no_solv, job.omp)

        # prepare the geometry
        if job.prepinfo[jobtype]["template"]:
            geom_line = next(inp.index(l) for l in inp if "{geom}" in l)
            inp.pop(geom_line)
            inp[geom_line:geom_line] = self.__prep_geom(
                job.conf, xyzfile, job.prepinfo["charge"], job.prepinfo["unpaired"]
            )
        else:
            inp.extend(
                self.__prep_geom(
                    job.conf, xyzfile, job.prepinfo["charge"], job.prepinfo["unpaired"]
                )
            )

        inp.extend(self.__prep_postgeom(job.prepinfo, jobtype))

        return [line + "\n" if not line.endswith("\n") else line for line in inp]

    def __prep_main(
        self, prepinfo: dict[str, Any], jobtype: str, orca5: bool
    ) -> list[str]:
        main = ["!"]
        pregeom = []

        # grab func, basis
        func = prepinfo[jobtype]["func_name"]
        basis = prepinfo[jobtype]["basis"]
        functype = prepinfo[jobtype]["func_type"]
        disp = prepinfo[jobtype]["disp"]

        if func != "dummy":
            main.append(func)

        if "composite" not in functype:
            main.append(basis)

        # set  RI def2/J,   RIJCOSX def2/J

        # Set def2/J in case of def2 basis
        if "def2" in basis.lower():
            main.append("def2/J")
        # Otherwise use autoaux
        else:
            main.append("autoaux")

        # settings for double hybrids
        if "double" in functype:
            main.append("RIJCOSX")

            if "nmr" in jobtype:
                main.append("NOFROZENCORE")
                pregeom.extend(["%mp2", "RI true", "end"])
            else:
                main.append("frozencore")
                pregeom.extend(["%mp2", "RI true", "end"])

            def2cbasis = ("def2-svp", "def2-tzvp", "def2-tzvpp", "def2-qzvpp")
            if basis.lower() in def2cbasis:
                main.append(f"{basis}/C")

            if not orca5:
                main.extend(["GRIDX6", "NOFINALGRIDX"])

        # settings for hybrids
        elif "hybrid" in functype:
            main.append("RIJCOSX")
            if not orca5:
                main.extend(["GRIDX6", "NOFINALGRIDX"])

        # settings for (m)ggas
        elif "gga" in functype:
            main.append("RI")

        # dummy type falls through every case, nothing is done in that case

        # use 'grid' setting from instructions to quickly configure the grid
        main.extend(self.__gridsettings[orca5][prepinfo[jobtype]["grid"]])

        # add dispersion
        # dispersion correction information
        # FIXME - temporary solution (not very nice)
        mapping = {
            "d3bj": "d3bj",
            "d3(0)": "D3ZERO",
            "d4": "d4",
            "nl": "NL",
        }

        if disp != "composite" and disp != "included":
            main.append(mapping.get(disp, ""))

        if disp == "nl" and not orca5:
            main.append("vdwgrid3")

        if "composite" not in functype:
            # try to apply gcp if basis set available
            # TODO - extend this
            gcp_keywords = {
                "minis": "MINIS",
                "sv": "SV",
                "6-31g(d)": "631GD",
                "def2-sv(p)": "SV(P)",
                "def2-svp": "SVP",
                "def2-tzvp": "TZ",
            }
            if prepinfo[jobtype]["gcp"]:
                if basis.lower() in gcp_keywords.keys():
                    main.append(f"GCP(DFT/{gcp_keywords[basis.lower()]})")

        # add job keyword for geometry optimizations
        # with ANCOPT
        if jobtype == "xtb_opt":
            main.extend(["ENGRAD", "tightSCF"])
        # for standard geometry optimization
        elif jobtype == "opt":
            main.extend(["OPT", "tightSCF"])

        # additional print settings
        if jobtype in ["xtb_opt", "opt"]:
            main.append("miniprint")
        else:
            main.append("printgap")

        return [" ".join(main)] + pregeom

    def __prep_pregeom(
        self,
        prepinfo: dict[str, Any],
        orca5: bool,
        jobtype: str,
        no_solv: bool,
        nprocs: int,
    ) -> list[str]:
        pregeom = []

        # Check ORCA version (important for grid keywords)
        if orca5 and nprocs > 1:
            pregeom.extend(["%pal", f"nprocs {nprocs}", "end"])

        # TODO: maybe limit TRAH macro steps (or when TRAH activates) to avoid
        # single jobs clogging everything up

        # set keywords for the selected solvent model
        if (
            not prepinfo["general"]["gas-phase"]
            and not no_solv
            and ("sm" in prepinfo[jobtype].keys())
        ):
            sm = prepinfo[jobtype]["sm"]
            solv_key = prepinfo[jobtype]["solvent_key_prog"]
            solv_key = f'"{solv_key}"'

            if sm == "smd":
                pregeom.extend(["%cpcm", "smd true", f"smdsolvent {solv_key}", "end"])
            elif sm == "cpcm":
                pregeom.insert(0, f"! CPCM({solv_key})")

        if jobtype == "uvvis":
            pregeom.extend(["%tddft", f"nroots {prepinfo['uvvis']['nroots']}", "end"])

        # Additional print settings
        if jobtype not in ["xtb_opt", "opt"]:
            pregeom.extend(["%output", "printlevel normal", "end"])

        if jobtype == "opt":
            if prepinfo[jobtype]["macrocycles"]:
                # Set max number of optimization cycles for ORCA driven optimization
                pregeom.extend(["%geom", f"maxiter {prepinfo['opt']['optcycles']}"])
            else:
                pregeom.extend(["%geom"])

            # Set optlevel
            mapping = {
                "crude": "loose",
                "sloppy": "loose",
                "loose": "loose",
                "lax": "loose",
                "normal": "loose",
                "tight": "normal",
                "vtight": "tight",
                "extreme": "tight",
            }

            # Try to apply literally first
            if prepinfo["opt"]["optlevel"] in ["loose", "normal", "tight"]:
                pregeom.extend([f"convergence {prepinfo['opt']['optlevel']}"])
            # Otherwise map to roughly corresponding orca optlevel
            else:
                pregeom.extend([f"convergence {mapping[prepinfo['opt']['optlevel']]}"])

            # Insert constraints (if provided)
            if prepinfo["opt"]["constraints"] is not None:
                with open(prepinfo["opt"]["constraints"], "r") as f:
                    constraints = f.readlines()
                pregeom.extend(constraints)

            pregeom.append("end")

        return pregeom

    def __prep_geom(
        self,
        conf: GeometryData,
        xyzfile: str,
        charge: int,
        unpaired: int,
    ) -> list[str]:
        lines = []
        # unpaired, charge, and coordinates
        # by default coordinates are written directly into input file
        if xyzfile is None:
            lines.extend([f"* xyz {charge} {unpaired + 1}"] + conf.toorca())
            lines.append("*")
        else:
            lines.extend([f"* xyzfile {charge} {unpaired + 1} {xyzfile}"])

        return lines

    def __prep_postgeom(
        self,
        prepinfo: dict[str, Any],
        jobtype: str,
    ) -> list[str]:
        lines = []
        # Set NMR parameters
        if "nmr" in jobtype:
            # Determine the settings that need to be put into the input file for the NMR calculation
            active_elements_map = {
                "H": prepinfo[jobtype]["h_active"],
                "C": prepinfo[jobtype]["c_active"],
                "F": prepinfo[jobtype]["f_active"],
                "Si": prepinfo[jobtype]["si_active"],
                "P": prepinfo[jobtype]["p_active"],
            }
            todo = [
                element for element, active in active_elements_map.items() if active
            ]

            todo2 = []
            todo3 = {}
            lines.append("%eprnmr")
            if jobtype.endswith("_s") or jobtype == "nmr":
                todo2.append("shift")
                lines.extend(
                    [
                        "origin giao",
                        "giao_2el giao_2el_same_as_scf",
                        "giao_1el giao_1el_analytic",
                    ]
                )
            if jobtype.endswith("_j") or jobtype == "nmr":
                if prepinfo[jobtype]["fc_only"]:
                    todo2.append("ssfc")
                else:
                    todo2.append("ssall")
                lines.append(f"SpinSpinRThresh {prepinfo[jobtype]['ss_cutoff']:.4f}")

            for i, element in enumerate(todo):
                lines.append(
                    f"Nuclei = all "
                    + element
                    + " { "
                    + ",".join(x for x in todo2)
                    + " }"
                )

            lines.append("end")

        return lines

    @staticmethod
    def __check_output(lines: list[str]) -> str | None:
        """
        Checks the lines from the output file for errors and returns them.

        Args:
            lines: list of lines from the output file.

        Returns:
            str | None: error message if an error was found, None otherwise
        """
        # Dict mapping specific messages from the output to error messages
        # TODO - this should be extended later
        out_to_err = {
            "SCF NOT CONVERGED": "scf_not_converged",
        }
        for line in lines:
            if any(key in line for key in out_to_err.keys()):
                # Returns the first error found
                key = next(filter(lambda x: x in line, out_to_err.keys()))
                return out_to_err[key]
        return None

    @staticmethod
    def __copy_mo(jobdir: str, filename: str, guess_file: str | tuple) -> None:
        """
        Copy MO file if possible (should be ORCA .gbw file).

        Args:
            jobdir: path to the job directory
            filename: name of the input file
            guess_file: path to the .gbw file to copy

        Returns:
            None
        """
        if guess_file is not None and type(guess_file) is not tuple:
            if os.path.isfile(guess_file) and ".gbw" in os.path.split(guess_file)[1]:
                if os.path.join(jobdir, f"{filename}.gbw") != guess_file:
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {guess_file}."
                    )
                    shutil.copy(guess_file, os.path.join(jobdir, f"{filename}.gbw"))

    def _sp(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="sp",
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[dict[str, float | None], dict[str, Any]]:
        """
        ORCA single-point calculation.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file
            no_solv: if True, no solvent model is used
            prep: if True, a new input file is generated (you only really want to make use of this for NMR)

        Returns:
            result (dict[str, float | None]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job

        result = {
            "energy": None,
        }
        """
        # set results
        result = {
            "energy": None,
        }

        meta = {
            "success": None,
            "error": None,
            "mo_path": None,
        }

        # set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        if prep:
            inp = self.__prep(job, "sp", no_solv=no_solv)

            # check for flags raised for this jobtype
            # NOTE: all other jobtypes call this function so flags are always checked here
            # except xtb_opt
            inp = self.__apply_flags(job, inp)

            # write input into file "{filename}.inp" in a subdir created for the
            # conformer
            with open(inputpath, "w") as f:
                f.writelines(inp)

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if self.copy_mo:
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # call orca
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("orca", call, outputpath, jobdir)
        # NOTE: using orca returncodes it is not possible to determine wether the calculation converged

        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")

        # read output
        with open(outputpath, "r", encoding=Config.CODING, newline=None) as out:
            lines = out.readlines()

        # Get final energy
        result["energy"] = next(
            (
                float(line.split()[4])
                for line in lines
                if "FINAL SINGLE POINT ENERGY" in line
            ),
            None,
        )

        # Check for errors in the output file in case returncode is 0
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and result["energy"] is not None
        else:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        if self.copy_mo:
            # store the path to the current .gbw file for this conformer if
            # possible
            if os.path.isfile(os.path.join(jobdir, f"{filename}.gbw")):
                meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        return result, meta

    def _gsolv(
        self, job: ParallelJob, jobdir: str
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """
        Calculates the solvation free enthalpy of a conformer using ORCA.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory

        Returns:
            result (dict[str, Any]): dictionary containing the results of the calculation
            meta (dict[str, Any]): metadata about the job

        result = {
            "gsolv": None,
            "energy_gas": None,
            "energy_solv": None,
        }
        """
        # what is returned in the end
        result = {
            "gsolv": None,
            "energy_gas": None,
            "energy_solv": None,
        }

        meta = {
            "success": None,
            "error": None,
            "mo_path": None,
        }

        # calculate gas phase
        spres, spmeta = self._sp(job, jobdir, filename="sp_gas", no_solv=True)

        if spmeta["success"]:
            result["energy_gas"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = spmeta["error"]
            return result, meta

        # calculate in solution
        spres, spmeta = self._sp(job, jobdir, filename="sp_solv")

        if spmeta["success"]:
            result["energy_solv"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = spmeta["error"]
            return result, meta

        if self.copy_mo:
            # store the path to the current .gbw file for this conformer if
            # possible
            if os.path.isfile(os.path.join(jobdir, "sp_solv.gbw")):
                meta["mo_path"] = os.path.join(jobdir, "sp_solv.gbw")

        # calculate solvation free enthalpy
        result["gsolv"] = result["energy_solv"] - result["energy_gas"]
        meta["success"] = True

        return result, meta

    def _opt(
        self, job: ParallelJob, jobdir: str, filename: str = "opt"
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """
        Geometry optimization using ORCA optimizer.
        Note that solvation in handled here always implicitly.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file

        Returns:
            result (dict[str, Any]): dictionary containing the results of the calculation
            meta (dict[str, Any]): metadata about the job
        """
        # prepare result
        # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
        # 'gncyc' contains the gradient norms for all cycles
        # 'energy' contains the final energy of the optimization (converged or unconverged)
        # 'geom' stores the optimized geometry in GeometryData.xyz format
        result = {
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "gncyc": None,
            "geom": None,
        }

        meta = {
            "success": None,
            "error": None,
            "mo_path": None,
        }

        # set orca input/output paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        inp = self.__prep(job, "opt")

        # check for flags raised for this jobtype
        # NOTE: all other jobtypes call this function so flags are always checked here
        # except xtb_opt
        inp = self.__apply_flags(job, inp)

        # write input into file "{filename}.inp" in a subdir created for the
        # conformer
        with open(inputpath, "w") as f:
            f.writelines(inp)

        # Get gbw files for initial guess
        if self.copy_mo:
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # call orca
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("orca", call, outputpath, jobdir)
        # NOTE: using orca returncodes it is not possible to determine wether the calculation converged

        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")

        # read output
        with open(outputpath, "r", encoding=Config.CODING, newline=None) as out:
            lines = out.readlines()

        meta["error"] = self.__check_output(lines)
        meta["success"] = meta["error"] is None

        # Check for errors in the output file in case returncode is 0
        if meta["success"]:
            # Check convergence
            if (
                next((True for x in lines if "OPTIMIZATION HAS CONVERGED" in x), None)
                is True
            ):
                result["converged"] = True
            else:
                result["converged"] = False

            # Get the number of cycles
            if result["converged"] is not None:
                for line in lines:
                    if "GEOMETRY OPTIMIZATION CYCLE" in line:
                        result["cycles"] = int(line.split()[4])

                # Get energies for each cycle
                result["ecyc"] = [
                    float(line.split("....")[-1].split()[0])
                    for line in filter(lambda x: "Current Energy" in x, lines)
                ]

                result["energy"] = result["ecyc"][-1]

                # Get all gradient norms for evaluation
                result["gncyc"] = [
                    float(line.split("....")[-1].split()[0])
                    for line in filter(lambda x: "Current gradient norm" in x, lines)
                ]

                # Get the last gradient norm
                result["grad_norm"] = result["gncyc"][-1]
                meta["success"] = True

                # Read out optimized geometry and update conformer geometry with this
                job.conf.fromxyz(os.path.join(jobdir, f"{filename}.xyz"))
                result["geom"] = job.conf.xyz
        elif meta["error"] is not None:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        if self.copy_mo:
            # store the path to the current .gbw file for this conformer if
            # possible
            if os.path.isfile(os.path.join(jobdir, f"{filename}.gbw")):
                meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        return result, meta

    # TODO - split this up
    def _xtb_opt(
        self, job: ParallelJob, jobdir: str, filename: str = "xtb_opt"
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """
        Geometry optimization using ANCOPT and ORCA gradients.
        Note that solvation is handled here always implicitly.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file

        Returns:
            result (dict[str, Any]): dictionary containing the results of the calculation
            meta (dict[str, Any]): metadata about the job

        Keywords required in prepinfo:
        - optcycles
        - hlow
        - optlevel
        - macrocycles
        - constraints

        result = {
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "geom": None,
        }
        """
        # NOTE: some "intuitivity problems":
        # the geometry of the conformer is written into a coord file and also into a xyz-file to be used by orca
        # xtb then outputs a file with the optimized geometry as 'xtbopt.coord', which is then read into the conformer
        # to update it's geometry

        # prepare result
        # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
        # 'gncyc' contains the gradient norms for all cycles
        # 'energy' contains the final energy of the optimization (converged or unconverged)
        # 'geom' stores the optimized geometry in GeometryData.xyz format
        result = {
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "gncyc": None,
            "geom": None,
        }

        meta = {
            "success": None,
            "error": None,
            "mo_path": None,
        }

        xcontrolname = "xtb_opt-xcontrol-inp"

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

        # write conformer geometry to coord file
        with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as file:
            file.writelines(job.conf.tocoord())

        # write xyz-file for orca
        with open(os.path.join(jobdir, f"{filename}.xyz"), "w", newline=None) as file:
            file.writelines(job.conf.toxyz())

        # set orca input path
        inputpath = os.path.join(jobdir, f"{filename}.inp")

        # prepare input dict
        inp = self.__prep(job, "xtb_opt")

        # check for flags raised for this jobtype
        # NOTE: all other jobtypes call this function so flags are always checked here
        # except xtb_opt
        inp = self.__apply_flags(job, inp)

        # write input into file "{filename}.inp" in a subdir created for the
        # conformer
        with open(inputpath, "w") as f:
            f.writelines(inp)

        # append some additional lines to the coord file for ancopt
        with open(
            os.path.join(jobdir, f"{filename}.coord"), "a", newline=None
        ) as newcoord:
            newcoord.writelines(
                [
                    "$external\n",
                    f"   orca input file= {filename}.inp\n",
                    f"   orca bin= {self._paths['orcapath']}\n",
                    "$end\n",
                ]
            )

        # prepare configuration file for ancopt (xcontrol file)
        with open(os.path.join(jobdir, xcontrolname), "w", newline=None) as out:
            out.write("$opt \n")
            if job.prepinfo["xtb_opt"]["macrocycles"]:
                out.write(f"maxcycle={job.prepinfo['xtb_opt']['optcycles']} \n")
                out.write(f"microcycle={job.prepinfo['xtb_opt']['optcycles']} \n")

            out.writelines(
                [
                    "average conv=true \n",
                    f"hlow={job.prepinfo['xtb_opt']['hlow']} \n",
                    "s6=30.00 \n",
                    "engine=lbfgs\n",
                    "$external\n",
                    f"   orca input file= {filename}.inp\n",
                    f"   orca bin= {self._paths['orcapath']} \n",
                ]
            )

            # Import constraints
            if job.prepinfo["xtb_opt"]["constraints"] is not None:
                with open(job.prepinfo["xtb_opt"]["constraints"], "r") as f:
                    lines = f.readlines()

                out.writelines(lines)

            out.write("$end \n")

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if self.copy_mo:
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # prepare xtb call
        call = [
            f"{filename}.coord",  # name of the coord file generated above
            "--opt",
            job.prepinfo["xtb_opt"]["optlevel"],
            "--orca",
            "-I",
            xcontrolname,
        ]

        # set path to the ancopt output file
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # call xtb
        returncode, errors = self._make_call("xtb", call, outputpath, jobdir)

        # check if optimization finished without errors
        # NOTE: right now, not converging scfs are not handled because returncodes need to be implemented first
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "unknown_error"
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            # TODO - the xtb returncodes should be handled
            return result, meta

        # read output
        with open(outputpath, "r", encoding=Config.CODING, newline=None) as file:
            lines = file.readlines()

        result["ecyc"] = []
        result["cycles"] = 0

        # Substrings indicating error in xtb
        error_ind = [
            "external code error",
            "|grad| > 500, something is totally wrong!",
            "abnormal termination of xtb",
        ]

        # Check if xtb terminated normally (if there are any error indicators
        # in the output)
        meta["success"] = (
            False
            if next((x for x in lines if any(y in x for y in error_ind)), None)
            is not None
            else True
        )
        if not meta["success"]:
            meta["error"] = "unknown_error"
            return result, meta

        # check convergence
        if (
            next((True for x in lines if "GEOMETRY OPTIMIZATION CONVERGED" in x), None)
            is True
        ):
            result["converged"] = True
        elif (
            next((True for x in lines if "FAILED TO CONVERGE GEOMETRY" in x), None)
            is True
        ):
            result["converged"] = False

        # Get the number of cycles
        if result["converged"] is not None:
            # tmp is one of the values from the dict defined below
            tmp = {
                True: ("GEOMETRY OPTIMIZATION CONVERGED", 5),
                False: ("FAILED TO CONVERGE GEOMETRY", 7),
            }
            tmp = tmp[result["converged"]]

            result["cycles"] = int(
                next(x for x in lines if tmp[0] in x).split()[tmp[1]]
            )

            # Get energies for each cycle
            result["ecyc"].extend(
                float(line.split("->")[-1])
                for line in filter(lambda x: "av. E: " in x, lines)
            )

            # Get all gradient norms for evaluation
            result["gncyc"] = [
                float(line.split()[3])
                for line in filter(lambda x: " gradient norm " in x, lines)
            ]

            # Get the last gradient norm
            result["grad_norm"] = result["gncyc"][-1]

            # store the final energy of the optimization in 'energy'
            result["energy"] = result["ecyc"][-1]
            meta["success"] = True

        if self.copy_mo:
            # store the path to the current .gbw file for this conformer
            meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # read out optimized geometry and update conformer geometry with this
        job.conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
        result["geom"] = job.conf.xyz

        # TODO - this might be a case where it would be reasonable to raise an
        # exception
        try:
            assert result["converged"] is not None
        except AssertionError:
            meta["success"] = False
            meta["error"] = "unknown_error"

        return result, meta

    def _nmr(
        self, job: ParallelJob, jobdir: str, filename: str = "nmr"
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """
        Calculate the NMR shieldings and/or couplings for a conformer using ORCA. ORCA gives only the active cores in the output
        so there is not need for more thinking here.
        Formatting:
            'shielding' contains a list of tuples (atom_index, shielding), with atom_index being the index of the atom
            in the internal coordinates of the GeometryData.
            'couplings' contains a list of tuples ((atom_index1, atom_index2), coupling), with the indices of the atoms
            in the internal coordinates of the GeometryData. A set is used to represent an atom pair and then converted
            to tuple to be serializable.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory

        Returns:
            result (dict[str, Any]): dictionary containing the results of the calculation
            meta (dict[str, Any]): metadata about the job
        """
        # Set results
        result = {
            "energy": None,
            "shieldings": None,
            "couplings": None,
        }

        meta = {
            "success": None,
            "error": None,
        }

        # If settings are the same, endings = [""], otherwise it contains "_s"
        # and/or "_j", depending on conditions
        for ending in [x[4:] for x in job.prepinfo.keys() if "nmr" in x]:
            # Set in/out path
            inputpath = os.path.join(jobdir, f"{filename}{ending}.inp")
            outputpath = os.path.join(jobdir, f"{filename}{ending}.out")

            # Prepare an input file for _sp
            inp = self.__prep(job, f"nmr{ending}")
            with open(inputpath, "w") as f:
                f.writelines(inp)

            # Run _sp using the input file generated above
            _, spmeta = self._sp(
                job, jobdir, filename=f"{filename}{ending}", prep=False
            )

            if not spmeta["success"]:
                meta["success"] = False
                meta["error"] = spmeta["error"]
                return result, meta

            # Grab shieldings and energy from the output
            with open(outputpath, "r") as f:
                lines = f.readlines()

            # Get final energy
            result["energy"] = next(
                (
                    float(line.split()[4])
                    for line in lines
                    if "FINAL SINGLE POINT ENERGY" in line
                ),
                None,
            )

            if result["energy"] is None:
                meta["success"] = False
                meta["error"] = "unknown_error"
                return result, meta

            # For shieldings watch out for the line "CHEMICAL SHIELDING SUMMARY
            # (ppm)"
            if ending in ["", "_s"]:
                start = (
                    lines.index(
                        next(x for x in lines if "CHEMICAL SHIELDING SUMMARY" in x)
                    )
                    + 6
                )

                result["shieldings"] = []

                for line in lines[start:]:
                    # Try to extract values from the lines, if that fails
                    # (probably IndexError) stop
                    try:
                        result["shieldings"].append(
                            (int(line.split()[0]), float(line.split()[2]))
                        )
                    except IndexError:
                        break

                # Sort shieldings by atom index
                result["shieldings"].sort(key=lambda x: x[0])

            if ending in ["", "_j"]:
                with open(outputpath, "r") as f:
                    lines = f.readlines()

                start = (
                    lines.index(next(x for x in lines if "SPIN-SPIN COUPLING" in x)) + 6
                )

                lines = lines[start:]

                end = (
                    lines.index(next(x for x in lines if "SUMMARY OF ISOTROPIC" in x))
                    - 3
                )

                lines = lines[:end]

                result["couplings"] = []

                # This goes through all the pairs, even though ORCA gives also non-unique pairs, since it makes life
                # easier (basically every element of a symmetric square matrix)
                # The iterator created here unpacks a tuple consisting of the multiples of three (indices of every
                # third line) and every third line in 'lines'
                for i, line in enumerate(lines):
                    if "NUCLEUS" in line:
                        # pair needs to be a frozenset because normal sets are not hashable and can therefore not be part
                        # of a normal set
                        pair = frozenset((int(line.split()[4]), int(line.split()[9])))
                        for line2 in lines[i:]:
                            if "Total" in line2 and "iso=" in line2:
                                coupling = float(line2.split()[5])
                                break

                        result["couplings"].append((pair, coupling))

                # Convert to set and back to get rid of duplicates
                # ('vectorizing the symmetric matrix')
                result["couplings"] = list(set(result["couplings"]))

                # Convert all the frozensets to a tuple to be serializable
                for i in range(len(result["couplings"])):
                    result["couplings"][i] = (
                        tuple(result["couplings"][i][0]),
                        result["couplings"][i][1],
                    )

                # Sort couplings by pairs
                result["couplings"].sort(key=lambda x: x[0])

        meta["success"] = True

        return result, meta

    def _uvvis(
        self, job: ParallelJob, jobdir: str, filename: str = "uvvis"
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """
        Run a single-point to calculate the oscillator strengths and excitation wavelengths.
        """
        # Set results
        result = {
            "energy": None,
            "excitations": None,
        }

        meta = {
            "success": None,
            "error": None,
        }

        # Set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare an input file for _sp
        inp = self.__prep(job, "uvvis")
        with open(inputpath, "w") as f:
            f.writelines(inp)

        # Run _sp using the input file generated above
        _, spmeta = self._sp(job, jobdir, filename=f"{filename}", prep=False)

        if not spmeta["success"]:
            meta["success"] = False
            meta["error"] = spmeta["error"]
            return result, meta

        # Grab excitations and energy from the output
        with open(outputpath, "r") as f:
            lines = f.readlines()

        # Get final energy
        result["energy"] = next(
            (
                float(line.split()[4])
                for line in lines
                if "FINAL SINGLE POINT ENERGY" in line
            ),
            None,
        )

        if result["energy"] is None:
            meta["success"] = False
            meta["error"] = "unknown_error"
            return result, meta

        # Find the line index of the actual excitations information
        start = lines.index(
            next(filter(lambda line: "ABSORPTION SPECTRUM" in line, lines))
        )

        # Get all the lines where the excitations are listed (+5 because of spacing in output file)
        uvvis_table = lines[start + 5 : start + 5 + job.prepinfo["uvvis"]["nroots"]]

        # Extract excitation wavelengths and oscillator strengths
        result["excitations"] = []
        for row in uvvis_table:
            spl = row.split()
            result["excitations"].append(
                {"wavelength": float(spl[-6]), "osc_str": float(spl[-5])}
            )

        meta["success"] = True

        return result, meta

    @staticmethod
    def __apply_flags(job: ParallelJob, input: list[str]) -> list[str]:
        """
        apply flags to an orca input

        this is very much work in progress
        """
        flag_to_setting = {
            "scf_not_converged": {
                "! veryslowconv",
                "%scf",
                "maxiter 300",
                "AutoTRAH false",
                "AutoTRAHIter 275",
                "CNVSOSCF true",
                "SOSCFStart 0.0002",
                "SOSCFMaxIt 200",
            },
        }
        flags = [flag for flag in job.flags.values() if flag in flag_to_setting]

        for flag in flags:
            # insert all settings dictated by the flags after the main
            input.extend(flag_to_setting[flag])

        return input


Factory.register_builder("orca", OrcaProc)
