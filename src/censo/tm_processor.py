"""
Contains TmProc class for calculating TURBOMOLE related properties of conformers.
"""
import subprocess
import os
import shutil
import math

from .qm_processor import QmProc
from .logging import setup_logger
from .parallel import ParallelJob
from .params import ENVIRON, ASSETS_PATH, WARNLEN, R, AU2KCAL
from .datastructure import GeometryData
from .utilities import frange

logger = setup_logger(__name__)


class TmProc(QmProc):
    """
    Performs calculations using TURBOMOLE.
    """
    _progname = "tm"

    __gridsettings = {
        "low": ["-grid", "m3", "-scfconv", "6"],
        "low+": ["-grid", "m4", "-scfconv", "6"],
        "high": ["-grid", "m4", "-scfconv", "7"],
        "high+": ["-grid", "m5", "-scfconv", "7"],
    }

    __returncode_to_err = {}

    # COSMO dielectric constants stored here
    # (mapped to solvent names under 'cosmo' in the solvent helper)
    __cosmo_dcs = {
        "acetone": 20.7,
        "acetonitrile": 36.6,
        "aniline": 6.9,
        "benzaldehyde": 18.2,
        "benzene": 2.3,
        "ccl4": 2.2,
        "ch2cl2": 9.1,
        "chcl3": 4.8,
        "cs2": 2.6,
        "cyclohexane": 2.0,
        "dichloroethane": 10.125,
        "diethylether": 4.4,
        "dioxane": 2.2,
        "dmf": 38.3,
        "dmso": 47.2,
        "ethanol": 24.6,
        "ethylacetate": 5.9,
        "furan": 3.0,
        "h2o": 80.1,
        "hexadecane": 2.1,
        "hexane": 1.9,
        "methanol": 32.7,
        "nitromethane": 38.2,
        "octane": 1.94,
        "octanol": 9.9,
        "phenol": 8.0,
        "thf": 7.6,
        "toluene": 2.4,
        "woctanol": 8.1
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # expand jobtypes with special turbomole jobtypes
        self._jobtypes = {
            **self._jobtypes,
            **{
                "sp": self._sp,
                "gsolv": self._gsolv,
                "xtb_opt": self._xtb_opt,
                "opt": self._opt,
                "nmr": self._nmr,
                "rot": self._rot,
            },
        }

        # Stores setting wether to copy MO-files for faster SCFs
        self.copy_mo: bool = False

    def __prep(self,
               job: ParallelJob,
               jobtype: str,
               jobdir: str,
               no_solv: bool = False) -> None:
        """
        Prepares TURBOMOLE input files using cefine for a specified jobtype.
        """
        # TODO - copy mo files
        func = job.prepinfo[jobtype]["func_name"]
        func_type = job.prepinfo[jobtype]["func_type"]
        basis = job.prepinfo[jobtype]["basis"]
        disp = job.prepinfo[jobtype]["disp"]

        # Set up basic cefine call
        call = [
            self._paths["cefinepath"], "-func", func, "-bas", basis, "-sym",
            "c1", "-noopt"
        ]

        # Configure grid
        call.extend(self.__gridsettings[job.prepinfo[jobtype]["grid"]])

        # r2scan-3c should use m4 grid
        if func == "r2scan-3c":
            call[call.index("m3")] = "m4"

        # Add dispersion
        # dispersion correction information
        # FIXME - temporary solution (not very nice)
        mapping = {
            "d3bj": "-d3",
            "d3(0)": "-zero",
            "d4": "-d4",
            "novdw": "-novdw",
            "included": "-novdw",
        }

        if disp not in ["composite", "nl"]:
            call.extend(mapping[disp])

        # Add charge and unpaired info
        if job.prepinfo["unpaired"] > 0:
            call.extend(["-uhf", f"{job.prepinfo['unpaired']}"])
        if job.prepinfo["charge"] != 0:
            call.extend(["-chrg", f"{job.prepinfo['charge']}"])

        # Call cefine
        cefine_output = subprocess.check_output(
            call,
            shell=False,
            text=True,
            stdin=None,
            stderr=subprocess.PIPE,
            cwd=jobdir,
            env=ENVIRON,
        ).decode("utf-8").splitlines()

        # TODO - Check output for errors
        for line in cefine_output:
            """
            if "define ended abnormally" in line:
                self.job["success"] = False
                broken = True
            elif "define_huge" in line:
                print(f"{'ERROR:':{WARNLEN}}define_huge: not found!")
                self.job["success"] = False
                broken = True
            elif "Could not find the beginning of the MO-eigenvalue data" in line:
                self.job["success"] = False
                broken = True
            """

        # Write coord file
        with open(os.path.join(jobdir, "coord"), "w") as f:
            f.writelines(job.conf.tocoord())

        # Do further manipulations of input files
        with open(os.path.join(jobdir, "control"), "r+") as f:
            lines = f.readlines()

            self.__prep_main(lines, func, disp, func_type, basis)
            if not no_solv and not job.prepinfo["general"]["gas-phase"]:
                self.__prep_solv(lines, job.prepinfo, jobtype)
            self.__prep_special(lines, job.prepinfo, jobtype)

            f.seek(0)  # Reset cursor to 0
            f.writelines(lines)
            f.truncate()  # Truncate in case the content is shorter than before

    def __prep_main(self, lines: list[str], func: str, disp: str,
                    func_type: str, basis: str):
        # Special treatment for KT1/KT2
        if "kt" in func:
            func_line_index = next(
                lines.index(l) for l in lines if "functional" in l)
            lines[func_line_index] = "   functional xcfun set-gga\n"
            lines.insert(func_line_index + 1,
                         f"   functional xcfun kt{func[2]} 1.0\n")
        # Special treatment for b97-3c
        elif func == "b97-3c":
            # Needs three-body dispersion
            disp_line_index = next(
                lines.index(l) for l in lines if "disp" in l)
            lines[disp_line_index] = "$disp3 -bj -abc\n"

        # Enable non local dispersion
        if disp == "nl":
            lines.insert(-1, "$donl\n")

        # Handle GCP
        if func_type != "composite":
            if basis.lower() == "def2-sv(p)":
                lines.insert(-1, "$gcp dft/sv(p)\n")
            else:
                lines.insert(-1,
                             f"$gcp dft/{basis.lower().replace('-', '')}\n")

    def __prep_solv(self, lines: list[str], prepinfo: dict[str, any],
                    jobtype: str):
        lines.insert(-1, "$cosmo\n")

        # write DC in any case
        lines.insert(-1,
                     f" epsilon= {self.__cosmo_dcs[prepinfo[jobtype]['solvent']]}")

        if prepinfo[jobtype]["sm"] == "dcosmors":
            # if using dcosmors also add the potential file path
            # NOTE: the value for solvent should never be None
            # (should be prevented in setup_prepinfo functions, as e.g. in optimizer.py)
            if prepinfo[jobtype]["solvent"] not in ["woctanol", "hexadecane", "octanol"]:
                lines.insert(-1,
                             f"$dcosmo_rs file={prepinfo[jobtype]['solvent']}_25.pot")
            else:
                # The three solvents above are specifically defined in the assets
                # TODO - this opens the possibility to insert your own potential files
                lines.insert(-1,
                             f"$dcosmo_rs file={os.path.join(ASSETS_PATH, prepinfo[jobtype]['solvent'])}_25.pot")

        if jobtype == "rot":
            lines[-1:-1] = [
                " cavity closed\n", " use_contcav\n", " nspa=272\n",
                " nsph=162\n", "$cosmo_isorad\n"
            ]

    def __prep_special(self, lines: list[str], prepinfo: dict[str, any],
                       jobtype: str):
        # Set NMR parameters
        if "nmr" in jobtype:
            # Determine the settings that need to be put into the input file for the NMR calculation
            active_elements_map = {
                '"H"': prepinfo[jobtype]["h_active"],
                '"C"': prepinfo[jobtype]["c_active"],
                '"F"': prepinfo[jobtype]["f_active"],
                '"Si"': prepinfo[jobtype]["si_active"],
                '"P"': prepinfo[jobtype]["p_active"],
            }

            todo = [
                element for element, active in active_elements_map.items()
                if active
            ]

            rpacor_line_index = next(
                lines.index(l) for l in lines if "rpacor" in l)
            rpacor = float(lines[rpacor_line_index].split()[-1])
            rpacor = rpacor if rpacor > 10000 else 10000
            lines[rpacor_line_index] = f"$rpacor {rpacor}\n"

            lines[-1:-1] = ["$ncoupling\n", " simple\n", " thr=0.0\n"]

            # nucsel only required if not all elements are active
            if not all(element in todo for element in active_elements_map):
                lines[-1:-1] = [
                    "$nucsel " + " ".join(todo) + "\n",
                    "$nucsel2 " + " ".join(todo) + "\n"
                ]

            lines.insert(-1, "$rpaconv 8\n")
        elif jobtype == "rot":
            # TODO
            """
            controlappend.append("$scfinstab dynpol nm")
            for i in self.job["freq_or"]:
                controlappend.append(f" {i}")  # e.g. 589
            controlappend.append("$velocity gauge")
            controlappend.append("$rpaconv 4")

            """
            raise NotImplementedError("Optical rotation not available yet!")

    @staticmethod
    def __copy_mo(jobdir: str, guess_file: str | tuple[str, str]) -> None:
        """
        Copy the MO file(s) for TURBOMOLE (should be TM format).
        """
        if guess_file is not None:
            if isinstance(guess_file, tuple):
                # open shell guess
                if all(os.path.isfile(f)
                       and any(g in f for g in ["alpha", "beta"])
                       and not any(os.path.join(jobdir, g) == guess_file for g in ["alpha", "beta"]) for f in guess_file):
                    # All MO files found and not already in dir
                    # Copy MO files
                    for g in ["alpha", "beta"]:
                        logger.debug(
                            f"{f'worker{os.getpid()}:':{WARNLEN}}Copying {g} file from {
                                guess_file}."
                        )
                        shutil.copy(guess_file,
                                    os.path.join(jobdir, g))
            else:
                # closed shell guess
                if (os.path.isfile(guess_file)
                    and os.path.split(guess_file)[1] == "mos"
                        and os.path.join(jobdir, "mos") != guess_file):
                    # Copy MO file
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying mos file from {
                            guess_file}."
                    )
                    shutil.copy(guess_file,
                                os.path.join(jobdir, "mos"))

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

    def _sp(
        self,
        job: ParallelJob,
        jobdir: str,
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        TURBOMOLE single-point calculation.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            no_solv: if True, no solvent model is used
            prep: if True, a new input file is generated

        Returns:
            result (dict[str, float | None]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job

        result = {i can
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
        outputpath = os.path.join(jobdir, "ridft.out")

        # check, if there is an existing mo/alpha,beta file and copy it if option
        # 'copy_mo' is true
        # mo files: mos/alpha,beta
        if self.copy_mo:
            self.__copy_mo(jobdir, job.mo_guess)

        # call turbomole
        call = ["ridft"]
        returncode, errors = self._make_call("tm", call, outputpath, jobdir)

        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(
                f"Job for {job.conf.name} failed. Stderr output:\n{errors}")

        with open(outputpath, "r") as f:
            lines = f.readlines()

        # Get final energy
        result["energy"] = next(
            (float(line.split()[4])
             for line in lines if "|  total energy      = " in line),
            None,
        )

        # Check for errors in the output file in case returncode is 0
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and result[
                "energy"] is not None
        else:
            meta["error"] = self.__returncode_to_err.get(
                returncode, "unknown_error")

        if self.copy_mo:
            # store the path to the current MO file(s) for this conformer if possible
            if os.path.isfile(os.path.join(jobdir, "mos")):
                meta["mo_path"] = os.path.join(jobdir, "mos")
            elif os.path.isfile(os.path.join(jobdir, "alpha")):
                meta["mo_path"] = (os.path.join(
                    jobdir, "alpha"), os.path.join(jobdir, "beta"))

        return result, meta

    def _gsolv(self, job: ParallelJob,
               jobdir: str) -> tuple[dict[str, any], dict[str, any]]:
        """
        Calculate the solvation contribution to the free enthalpy explicitely using (D)COSMO(RS).

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory

        Returns:
            result (dict[str, any]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job
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

        if job.prepinfo["sp"]["sm"] in ["cosmo", "dcosmors"]:
            # Non-COSMORS procedure:
            # Run gas-phase sp
            spres, spmeta = self._sp(job, jobdir, no_solv=True)

            if spmeta["success"]:
                result["energy_gas"] = spres["energy"]
            else:
                meta["success"] = False
                meta["error"] = spmeta["error"]
                return result, meta

            # Run solution sp
            spres, spmeta = self._sp(job, jobdir)

            if spmeta["success"]:
                result["energy_solv"] = spres["energy"]
            else:
                meta["success"] = False
                meta["error"] = spmeta["error"]
                return result, meta

            # Calculate gsolv from energy difference
            result["gsolv"] = result["energy_solv"] - result["energy_gas"]
            meta["success"] = True
        else:
            # COSMORS procedure:
            # Run gas-phase sp with unaltered settings
            spres, spmeta = self._sp(job, jobdir, no_solv=True)

            if spmeta["success"]:
                result["energy_gas"] = spres["energy"]
            else:
                meta["success"] = False
                meta["error"] = spmeta["error"]
                return result, meta

            # Run gas-phase sp with BP86 and def2-TZVP (normal)/def2-TZVPD (fine)
            job.prepinfo["sp"]["func_name"] = "b-p"
            job.prepinfo["sp"]["disp"] = "novdw"

            if job.prepinfo["sp"]["sm"] == "cosmors-fine":
                job.prepinfo["sp"]["basis"] = "def2-tzvpd"
            else:
                job.prepinfo["sp"]["basis"] = "def2-tzvp"

            spres, spmeta = self._sp(job, jobdir, no_solv=True)

            if not spmeta["success"]:
                meta["success"] = False
                meta["error"] = spmeta["error"]
                return result, meta

            # Run special cosmo sp with BP86 and def2-TZVP (normal)/def2-TZVPD (fine)
            self.__prep(job, "sp", jobdir, no_solv=True)

            # Write special settings for cosmo into control file
            with open(os.path.join(jobdir, "control"), "r+") as f:
                lines = f.readlines()

                lines[-1:-1] = ["$cosmo\n", " epsilon=infinity\n", " use_contcav\n",
                                " cavity closed\n", " nspa=272\n", " nsph=162\n", "$cosmo_out  file=out.cosmo\n"]

                f.seek(0)
                f.writelines(lines)
                f.truncate()

            # Run sp
            spres, spmeta = self._sp(job, jobdir, prep=False)

            if not spmeta["success"]:
                meta["success"] = False
                meta["error"] = spmeta["error"]
                return result, meta

            # Prepare cosmotherm.inp
            lines = [self._paths["cosmorssetup"] +
                     "\n", "EFILE VPFILE\n", "!!\n"]
            if job.prepinfo["sp"]["solvent_key_prog"] == "woctanol":
                lines.extend(
                    [f"f = h2o.cosmo fdir={self._paths['dbpath']} autoc\n",
                     f"f = 1-octanol.cosmo fdir={self._paths['dbpath']} autoc\n"]
                )
                mix = "0.27 0.73"
            else:
                lines.append(
                    f"f = {job.prepinfo['sp']['solvent_key_prog']}.cosmo"
                    + f" fdir={self._paths['dbpath']} autoc\n"
                )
                mix = "1.0 0.0"

            lines.append("f = out.cosmo\n")

            if job.prepinfo["general"]["multitemp"]:
                trange = frange(
                    job.prepinfo["general"]["trange"][0],
                    job.prepinfo["general"]["trange"][1],
                    step=job.prepinfo["general"]["trange"][2],
                )

                # Always append the fixed temperature to the trange so that it is the last value
                trange.append(job.prepinfo["general"]["temperature"])

                # Write trange to the xcontrol file
                for t in trange:
                    lines.append(
                        f"henry xh={{mix}} tc={t - 273.15} Gsolv\n")
            else:
                lines.append(
                    f"henry xh={{mix}} tc={job.prepinfo['general']['temperature'] - 273.15} Gsolv\n")

            # TODO - Run cosmotherm
            returncode, errors = self._make_call()

            meta["success"] = returncode == 0
            if not meta["success"]:
                logger.warning(
                    f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
                meta["error"] = "unknown_error"
                return result, meta

            # Read output
            gsolvt = {}
            videal = 24.789561955 / 298.15  # molar volume for ideal gas at 298.15 K 100.0 kPa

            cosmothermtab = os.path.join(jobdir, "cosmotherm.tab")
            with open(cosmothermtab, "r") as inp:
                lines = inp.readlines()
            for line in lines:
                if "T=" in line:
                    T = float(line.split()[5])

                    # Add volume work
                    vwork = R * T * math.log(videal * T)
                elif " out " in line:
                    # Add volume work
                    gsolvt[T] = float(line.split()[-1]) / \
                        AU2KCAL + vwork / AU2KCAL

            result["gsolvt"] = gsolvt
            result["gsolv"] = gsolvt[job.prepinfo["general"]["temperature"]]
            result["energy_solv"] = result["energy_gas"] + result["gsolv"]

            # cosmothermd
            with open(os.path.join(jobdir, "cosmors.out"), "w",) as out:
                T = job.prepinfo["general"]["temperature"]
                vwork = R * T * math.log(videal * T)

                out.writelines([
                    "This is cosmothermrd (python version in ENSO) (SG,FB,SAW, 06/18)\n",
                    "final thermochemical solvation properties in kcal/mol\n"
                    "----------------------------------------------------------\n",
                    " Gsolv({} K)= {:10.3f}\n".format(
                        T, result["gsolv"] *
                        AU2KCAL - vwork),
                    " VWork({} K)= {:10.3f}\n".format(
                        T, vwork),
                    " Gsolv+VWork({} K)= {:10.3f}\n".format(
                        # volwork already included!
                        T, result["gsolv"] * AU2KCAL)
                ])

        return result, meta

    def _xtb_opt(self,
                 job: ParallelJob,
                 jobdir: str,
                 filename: str = "xtb_opt"
                 ) -> tuple[dict[str, any], dict[str, any]]:
        """
        Geometry optimization using ANCOPT and ORCA gradients.
        Note that solvation is handled here always implicitly.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file

        Returns:
            result (dict[str, any]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job

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
        with open(os.path.join(jobdir, f"{filename}.coord"), "w",
                  newline=None) as file:
            file.writelines(job.conf.tocoord())

        # prepare configuration file for ancopt (xcontrol file)
        with open(os.path.join(jobdir, xcontrolname), "w",
                  newline=None) as out:
            out.write("$opt \n")
            if job.prepinfo["xtb_opt"]["macrocycles"]:
                out.write(
                    f"maxcycle={job.prepinfo['xtb_opt']['optcycles']} \n")
                out.write(
                    f"microcycle={job.prepinfo['xtb_opt']['optcycles']} \n")

            out.writelines([
                "average conv=true \n",
                f"hlow={job.prepinfo['xtb_opt']['hlow']} \n",
                "s6=30.00 \n",
                # remove unnecessary sp/gradient call in xTB
                "engine=lbfgs\n",
                "$external\n",
                f"   orca input file= {filename}.inp\n",
                f"   orca bin= {self._paths['orcapath']} \n",
            ])

            # Import constraints
            if job.prepinfo["xtb_opt"]["constraints"] is not None:
                with open(job.prepinfo["xtb_opt"]["constraints"], "r") as f:
                    lines = f.readlines()

                out.writelines(lines)

            out.write("$end \n")

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if self.copy_mo:
            self.__copy_mo(jobdir, job.mo_guess)

        # prepare xtb call
        call = [
            f"{filename}.coord",  # name of the coord file generated above
            "--opt",
            job.prepinfo["xtb_opt"]["optlevel"],
            "--tm",
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
            logger.warning(
                f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            # TODO - the xtb returncodes should be handled
            return result, meta

        # read output
        with open(outputpath, "r") as file:
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
        meta["success"] = (False if next((x for x in lines if any(
            y in x for y in error_ind)), None) is not None else True)
        if not meta["success"]:
            meta["error"] = "unknown_error"
            return result, meta

        # check convergence
        if (next(
            (True for x in lines if "GEOMETRY OPTIMIZATION CONVERGED" in x),
                None) is True):
            result["converged"] = True
        elif (next((True for x in lines if "FAILED TO CONVERGE GEOMETRY" in x),
                   None) is True):
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
                next(x for x in lines if tmp[0] in x).split()[tmp[1]])

            # Get energies for each cycle
            result["ecyc"].extend(
                float(line.split("->")[-1])
                for line in filter(lambda x: "av. E: " in x, lines))

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
            # store the path to the current MO file(s) for this conformer if possible
            if os.path.isfile(os.path.join(jobdir, "mos")):
                meta["mo_path"] = os.path.join(jobdir, "mos")
            elif os.path.isfile(os.path.join(jobdir, "alpha")):
                meta["mo_path"] = (os.path.join(
                    jobdir, "alpha"), os.path.join(jobdir, "beta"))

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

    def _opt(self):
        raise NotImplementedError(
            "Pure TURBOMOLE geometry optimization not available yet.")

    def _nmr(self,
             job: ParallelJob,
             jobdir: str) -> tuple[dict[str, any], dict[str, any]]:
        """
        Calculate the NMR shieldings and/or couplings for a conformer using TURBOMOLE.
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
            result (dict[str, any]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job
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

        # Run sp first
        spres, spmeta = self._sp(job, jobdir)

        if not spmeta["success"]:
            meta["success"] = False
            meta["error"] = spmeta["error"]
            return result, meta
        
        result["energy"] = spres["energy"]

        # Run shielding calculations if requested
        if "nmr_s" in job.prepinfo or "nmr" in job.prepinfo:
            # Set in/out path
            outputpath = os.path.join(jobdir, "mpshift.out")

            call = [self._paths["mpshiftpath"], "-smpcpus", f"{job.omp}"]

            # Run mpshift for shielding calculation
            returncode, errors = self._make_call("tm", call, outputpath, jobdir)

            meta["success"] = returncode == 0
            if not meta["success"]:
                logger.warning(
                    f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
                meta["error"] = "unknown_error"
                return result, meta

            # Grab shieldings from the output
            with open(outputpath, "r") as f:
                lines = f.readlines()

            start = lines.index(
                next(x for x in lines
                     if ">>>>> DFT MAGNETIC SHIELDINGS <<<<<" in x))

            lines = lines[start:]

            result["shieldings"] = []

            # Get lines with "ATOM" in it
            line_indices = [lines.index(l) for l in lines if "ATOM" in l]

            for i in line_indices:
                split = lines[i].split()
                result["shieldings"].append((int(split[2]), float(split[4])))

            # Sort shieldings by atom index
            result["shieldings"].sort(key=lambda x: x[0])

        # Run couplings calculation if requested
        if "nmr_j" in job.prepinfo or "nmr" in job.prepinfo:
            # Set in/out path
            outputpath = os.path.join(jobdir, "escf.out")

            call = [self._paths["escfpath"], "-smpcpus", f"{job.omp}"]

            # Run escf for couplings calculation
            returncode, errors = self._make_call("tm", call, outputpath, jobdir)

            meta["success"] = returncode == 0
            if not meta["success"]:
                logger.warning(
                    f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
                meta["error"] = "unknown_error"
                return result, meta

            # Grab couplings from the output
            with open(outputpath, "r") as f:
                lines = f.readlines()

            start = lines.index(next(x for x in lines if "Nuclear coupling constants" in x)) + 3

            lines = lines[start:]

            end = lines.index(next(x for x in lines if "-----" in x))

            lines = lines[:end]

            result["couplings"] = []

            line_indices = [lines.index(l) for l in lines if len(l.split()) in [6, 7]]

            for i in line_indices:
                # pair needs to be a frozenset because normal sets are not hashable and can therefore not be part
                # of a normal set
                split = lines[i].split()
                pair = frozenset((int(split[1]), int(split[4].split(":")[0])))
                coupling = float(split[5])
                result["couplings"].append((pair, coupling))

            # Convert to set and back to get rid of duplicates
            # ('vectorizing the symmetric matrix')
            result["couplings"] = list(set(result["couplings"]))

            # Convert all the frozensets to a tuple to be serializable
            for i in range(len(result["couplings"])):
                result["couplings"][i] = (tuple(result["couplings"][i][0]),
                                          result["couplings"][i][1])

            # Sort couplings by pairs
            result["couplings"].sort(key=lambda x: x[0])

        meta["success"] = True

        return result, meta

    def _rot(self):
        pass