"""
Contains TmProc class for calculating TURBOMOLE related properties of conformers.
"""
import subprocess
import os
import shutil

from .qm_processor import QmProc
from .logging import setup_logger
from .parallel import ParallelJob
from .params import ENVIRON, ASSETS_PATH
from .datastructure import GeometryData

logger = setup_logger(__name__)


class TmProc(QmProc):
    """
    Performs calculations using TURBOMOLE.
    """

    __gridsettings = {
        "low": ["-grid", "m3", "-scfconv", "6"],
        "low+": ["-grid", "m4", "-scfconv", "6"],
        "high": ["-grid", "m4", "-scfconv", "7"],
        "high+": ["-grid", "m5", "-scfconv", "7"],
    }

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
        """self._jobtypes = {
            **self._jobtypes, **{
                "opt": self._opt,
                "uvvis": self._uvvis,
            }
        }"""

    def __prep(self,
               job: ParallelJob,
               jobtype: str,
               jobdir: str,
               no_solv: bool = False) -> list:
        """
        Prepares TURBOMOLE input files using cefine for a specified jobtype.
        """
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

        # Check output for errors
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
            if not no_solv:
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
        # NOTE: the latter is basically given due to compatibility checks
        if not prepinfo["general"]["gas-phase"] and prepinfo[jobtype]["sm"] in ("cosmo", "dcosmors"):
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
            """
            controlappend.append("$scfinstab dynpol nm")
            for i in self.job["freq_or"]:
                controlappend.append(f" {i}")  # e.g. 589
            controlappend.append("$velocity gauge")
            controlappend.append("$rpaconv 4")

            """
            raise NotImplementedError("Optical rotation not available yet!")

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

        # check, if there is an existing mo/alpha,beta file and copy it if option
        # 'copy_mo' is true
        # mo files: mos/alpha,beta
        if self.copy_mo:
            if job.mo_guess is not None and os.path.isfile(job.mo_guess):
                if os.path.join(jobdir, f"{filename}.gbw") != job.mo_guess:
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {
                            job.mo_guess}."
                    )
                    shutil.copy(job.mo_guess,
                                os.path.join(jobdir, f"{filename}.gbw"))

        return result, meta

    def _gsolv(self, job: ParallelJob,
               jobdir: str) -> tuple[dict[str, any], dict[str, any]]:
        """
        Calculate the solvation contribution to the free enthalpy explicitely using (D)COSMO(RS).
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

            # Run solution sp
            spres, spmeta = self._sp(job, jobdir)

            # Calculate gsolv from energy difference
        else:
            # COSMORS procedure:
            # Run gas-phase sp with unaltered settings
            # Run gas-phase sp with BP86 and def2-TZVP (normal)/def2-TZVPD (fine)
            # Run special cosmo sp with BP86 and def2-TZVP (normal)/def2-TZVPD (fine)
            """
            out.write("$cosmo \n")
            out.write(" epsilon=infinity \n")
            out.write(" use_contcav \n")
            out.write(" cavity closed \n")
            out.write(" nspa=272 \n")
            out.write(" nsph=162 \n")
            out.write("$cosmo_out file=out.cosmo \n")
            out.write("$end \n")
            """
            # Prepare cosmotherm.inp
            # Run cosmotherm
            # Add volume work

        return result, meta

    def _xtb_opt(self):
        pass

    def _opt(self):
        pass

    def _nmr(self):
        pass

    def _rot(self):
        pass
