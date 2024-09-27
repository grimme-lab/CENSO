"""
Contains TmProc class for calculating TURBOMOLE related properties of conformers.
"""
import subprocess
import os

from .qm_processor import QmProc
from .logging import setup_logger
from .parallel import ParallelJob
from .params import ENVIRON
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
               no_solv: bool = False,
               coordfile: str = None) -> list:
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

        # Do further manipulations of input files
        with open(os.path.join(jobdir, "control"), "r+") as f:
            lines = f.readlines()

            self.__prep_main(lines, func, disp, func_type, basis)
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
        if prepinfo[jobtype]["sm"] in ("cosmors", "dcosmors"):
            lines[-1:-1] = ["$cosmo\n", f" {prepinfo}\n"]  # TODO TODO TODO

            if jobtype == "rot":
                # TODO - disable isorad option?
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
            raise NotImplementedError

    def _sp(self):
        # check, if there is an existing mo/alpha,beta file and copy it if option
        # 'copy_mo' is true
        if self.copy_mo:
            if job.mo_guess is not None and os.path.isfile(job.mo_guess):
                if os.path.join(jobdir, f"{filename}.gbw") != job.mo_guess:
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {job.mo_guess}."
                    )
                    shutil.copy(job.mo_guess,
                                os.path.join(jobdir, f"{filename}.gbw"))

    def _gsolv(self):
        pass

    def _xtb_opt(self):
        pass

    def _opt(self):
        pass

    def _nmr(self):
        pass

    def _rot(self):
        pass
