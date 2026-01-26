"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""

import os
from pathlib import Path
import shutil
from typing import final, cast, override

from ..config.job_config import (
    NMRJobConfig,
    OptJobConfig,
    SPJobConfig,
    UVVisJobConfig,
    XTBOptJobConfig,
)
from .job import (
    JobContext,
)
from .results import (
    GsolvResult,
    NMRResult,
    OptResult,
    SPResult,
    MetaData,
    UVVisResult,
)
from ..utilities import Factory
from ..logging import setup_logger
from ..molecules import GeometryData
from ..params import (
    WARNLEN,
    USER_ASSETS_PATH,
    OrcaSolvMod,
    Prog,
)
from ..assets import FUNCTIONALS, SOLVENTS
from .qm_processor import QmProc

logger = setup_logger(__name__)


@final
class OrcaProc(QmProc):
    """
    Performs calculations with ORCA.
    """

    progname = Prog.ORCA

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
        30: "input_error",
        52: "input_error",
        55: "unknown_error",
        125: "unknown_error",
    }

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

        # check ORCA version (orca5 = True means at least ORCA version 5)
        # orca5 = int(config.paths.orcaversion[0]) > 4

        inp: list[str] = []

        template = "template" in config.model_fields and bool(config.template)
        if template:
            # NOTE: when using templates we're not going to check for double definitions!
            # if the template is messed up, orca will fail and the user should deal with that
            # load template file
            try:
                inp = (
                    (USER_ASSETS_PATH / f"{job.from_part}.orca.template")
                    .read_text()
                    .split("\n")
                )
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not find template file {job.from_part}.orca.template."
                )

        # prepare the main line of the orca input
        if template:
            main_line = next(inp.index(line) for line in inp if "{main}" in line)
            inp.pop(main_line)
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

    def __prep_main(
        self, jobtype: str, config: SPJobConfig, no_solv: bool
    ) -> list[str]:
        orca5 = int(config.paths.orcaversion[0]) > 4

        # grab func, basis
        func = config.func
        basis = config.basis

        # All of this should be validated while setting up the config
        funcname = FUNCTIONALS[func][Prog.ORCA.value]
        functype = FUNCTIONALS[func]["type"]
        disp = FUNCTIONALS[func]["disp"]

        main = ["!", funcname]
        pregeom = []

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
            if "nmr" in jobtype:
                main.append("NOFROZENCORE")
                pregeom.extend(["%mp2", "RI true", "density relaxed", "end"])
            else:
                main.append("frozencore")
                pregeom.extend(["%mp2", "RI true", "end"])

            def2cbasis = ("def2-svp", "def2-tzvp", "def2-tzvpp", "def2-qzvpp")
            if basis.lower() in def2cbasis:
                main.append(f"{basis}/C")
            else:
                main.append("def2-tzvpp/C")

            if not orca5:
                main.extend(["GRIDX6", "NOFINALGRIDX"])
        # settings for hybrids
        elif "hybrid" in functype:
            if not orca5:
                main.extend(["GRIDX6", "NOFINALGRIDX"])

        # dummy type falls through every case, nothing is done in that case

        # use 'grid' setting from instructions to quickly configure the grid
        main.extend(self.__gridsettings[orca5][config.grid])

        # add dispersion
        # dispersion correction information
        # TODO: temporary solution (not very nice)
        mapping = {
            "d3bj": "d3bj",
            "d3(0)": "D3ZERO",
            "d4": "d4",
            "nl": "NL",
            "novdw": "",
        }

        if disp != "composite" and disp != "included":
            main.append(mapping[disp])

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
            if basis.lower() in gcp_keywords.keys():
                main.append(f"GCP(DFT/{gcp_keywords[basis.lower()]})")

        # add job keyword for geometry optimizations
        # with ANCOPT
        if jobtype == "xtb_opt":
            main.extend(["ENGRAD", "tightSCF"])
        # for standard geometry optimization
        elif jobtype == "opt":
            main.extend(["OPT", "tightSCF"])

        # set keywords for the selected solvent model
        if not no_solv:
            assert config.solvent
            assert config.sm
            sm = config.sm
            solv_key = f"{SOLVENTS[config.solvent][sm]}"

            orca6 = int(config.paths.orcaversion[0]) > 5
            if orca6 or sm != OrcaSolvMod.SMD:
                main.append(f"{sm.upper()}({solv_key})")
            else:
                pregeom.extend(["%cpcm", "smd true", f'smdsolvent "{solv_key}"', "end"])

        # additional print settings
        if "opt" in jobtype:
            main.append("miniprint")
        else:
            main.append("printgap")

        return [" ".join(main)] + pregeom

    def __prep_pregeom(
        self,
        jobtype: str,
        config: SPJobConfig,
        nprocs: int,
    ) -> list[str]:
        orca5 = int(config.paths.orcaversion[0]) > 4
        pregeom = []

        # Check ORCA version (important for grid keywords)
        if orca5 and nprocs > 1:
            pregeom.extend(["%pal", f"nprocs {nprocs}", "end"])

        # TODO: maybe limit TRAH macro steps (or when TRAH activates) to avoid
        # single jobs clogging everything up

        if jobtype == "uvvis":
            assert isinstance(config, UVVisJobConfig)
            pregeom.extend(["%tddft", f"nroots {config.nroots}", "end"])

        # Additional print settings
        if "opt" not in jobtype:
            pregeom.extend(["%output", "printlevel normal", "end"])

        if jobtype == "opt":
            assert isinstance(config, OptJobConfig)
            if config.macrocycles:
                # Set max number of optimization cycles for ORCA driven optimization
                pregeom.extend(["%geom", f"maxiter {config.optcycles}"])
            else:
                pregeom.extend(["%geom"])

            # Set optlevel
            # NOTE: this is only approximatively equivalent,
            # convergence criteria differ too much between ANCOPT and ORCA
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
            if config.optlevel in ["loose", "normal", "tight"]:
                pregeom.extend([f"convergence {config.optlevel}", "end"])
            # Otherwise map to roughly corresponding orca optlevel
            else:
                pregeom.extend([f"convergence {mapping[config.optlevel]}", "end"])

            # Insert constraints (if provided)
            # FIXME: without better parser this will not work
            """
            if prepinfo["opt"]["constraints"] is not None:
                parser = OrcaParser()
                constraints = parser.read_input(prepinfo["opt"]["constraints"])
                indict["geom"].update(constraints["geom"])
            """

        return pregeom

    def __prep_geom(
        self,
        conf: GeometryData,
        xyzfile: str | None,
        charge: int,
        unpaired: int,
    ) -> list[str]:
        lines: list[str] = []
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
        config: SPJobConfig,
        jobtype: str,
    ) -> list[str]:
        lines = []
        # Set NMR parameters
        if "nmr" in jobtype:
            assert isinstance(config, NMRJobConfig)
            # Determine the settings that need to be put into the input file for the NMR calculation
            todo = config.active_nuclei.split(",")

            todo2 = []
            lines.append("%eprnmr")
            if config.shieldings:
                todo2.append("shift")
                lines.extend(
                    [
                        "origin giao",
                        "giao_2el giao_2el_same_as_scf",
                        "giao_1el giao_1el_analytic",
                    ]
                )
            if config.couplings:
                if config.fc_only:
                    todo2.append("ssfc")
                else:
                    todo2.append("ssall")
                lines.append(f"SpinSpinRThresh {config.ss_cutoff:.4f}")

            lines.extend(
                [
                    f"Nuclei = all {element.capitalize()} "
                    + "{ "
                    + ",".join(todo2)
                    + " }"
                    for element in todo
                ]
            )
            lines.append("end")

        return lines

    @staticmethod
    def __check_output(lines: list[str]) -> str:
        """
        Checks the lines from the output file for errors and returns them.

        :param lines: list of lines from the output file.
        :type lines: list[str]
        :returns: error message if an error was found, empty string otherwise
        :rtype: str
        """
        # Dict mapping specific messages from the output to error messages
        # TODO: this should be extended later
        out_to_err = {
            "SCF NOT CONVERGED": "scf_not_converged",
        }
        for line in lines:
            if any(key in line for key in out_to_err.keys()):
                # Returns the first error found
                key = next(filter(lambda x: x in line, out_to_err.keys()))
                return out_to_err[key]
        return ""

    @staticmethod
    def __copy_mo(jobdir: str | Path, filename: str, guess_file: str | Path) -> None:
        """
        Copy MO file if possible (should be ORCA .gbw file).

        :param jobdir: path to the job directory
        :type jobdir: str | Path
        :param filename: name of the input file
        :type filename: str
        :param guess_file: path to the .gbw file to copy
        :type guess_file: str | Path
        :returns: None
        :rtype: None
        """
        if os.path.isfile(guess_file) and ".gbw" in os.path.split(guess_file)[1]:
            if os.path.join(jobdir, f"{filename}.gbw") != guess_file:
                logger.debug(
                    f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {guess_file}."
                )
                shutil.copy(guess_file, os.path.join(jobdir, f"{filename}.gbw"))

    @final
    def _sp(
        self,
        job: JobContext,
        config: SPJobConfig,
        jobdir: str | Path | None = None,
        filename: str = "sp",
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[SPResult, MetaData]:
        """
        ORCA single-point calculation.
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
        :param prep: if True, a new input file is generated (you only really want to make use of this for NMR)
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

        # prepare input
        if prep:
            inp = self.__prep(job, config, "sp", no_solv=config.gas_phase or no_solv)

            # write input into file "{filename}.inp" in a subdir created for the
            # conformer
            with open(inputpath, "w") as f:
                f.writelines(inp)

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if (
            config.copy_mo
            and job.mo_guess is not None
            and not isinstance(job.mo_guess, tuple)
        ):
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # call orca
        call = [config.paths.orca, f"{filename}.inp"]
        returncode, _ = self._make_call(call, outputpath, jobdir)
        # NOTE: using orca returncodes it is not possible to determine wether the calculation converged

        meta.success = returncode == 0

        # read output
        with open(outputpath) as out:
            lines = out.readlines()

        # Check for errors in the output file in case returncode is 0
        if meta.success:
            meta.error = self.__check_output(lines)
            if meta.error != "":
                return result, meta
        else:
            meta.error = self.__returncode_to_err.get(returncode, "unknown_error")
            return result, meta

        # Get final energy
        try:
            result.energy = next(
                (
                    float(line.split()[4])
                    for line in lines
                    if "FINAL SINGLE POINT ENERGY" in line
                ),
            )
        except StopIteration:
            meta.success = False
            meta.error = "Could not parse final energy"

        if config.copy_mo:
            # store the path to the current .gbw file for this conformer if possible
            if os.path.isfile(os.path.join(jobdir, f"{filename}.gbw")):
                result.mo_path = os.path.join(jobdir, f"{filename}.gbw")

        return result, meta

    @override
    def sp(self, *args, **kwargs):
        """
        Perform single-point calculation.

        :param args: Arguments.
        :param kwargs: Keyword arguments.
        :return: Tuple of (SP result, metadata).
        """
        return self._sp(*args, **kwargs)

    @override
    def gsolv(
        self, job: JobContext, config: SPJobConfig
    ) -> tuple[GsolvResult, MetaData]:
        """
        Calculates the solvation free enthalpy of a conformer using ORCA.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :param config: SP configuration
        :return: Tuple of (GsolvResult, MetaData)
        """
        # Check required settings

        # what is returned in the end
        result = GsolvResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "gsolv")

        # calculate gas phase
        spres, spmeta = self._sp(
            job, config, jobdir=jobdir, filename="sp_gas", no_solv=True
        )

        if spmeta.success:
            result.energy_gas = spres.energy
        else:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # calculate in solution
        spres, spmeta = self._sp(job, config, jobdir=jobdir, filename="sp_solv")

        if spmeta.success:
            result.energy_solv = spres.energy
        else:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        if config.copy_mo:
            # store the path to the current .gbw file for this conformer if
            # possible
            assert isinstance(spres.mo_path, str | Path)
            if os.path.isfile(spres.mo_path):
                result.mo_path = spres.mo_path

        # calculate solvation free enthalpy
        result.gsolv = result.energy_solv - result.energy_gas
        meta.success = True

        return result, meta

    @override
    def opt(
        self,
        job: JobContext,
        config: OptJobConfig,
    ) -> tuple[OptResult, MetaData]:
        """
        Geometry optimization using ORCA optimizer.
        Note that solvation in handled here always implicitly.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :param config: Optimization configuration
        :return: Tuple of (OptResult, MetaData)
        """
        # prepare result
        result = OptResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "opt")
        filename = "opt"

        # set orca input/output paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        # TODO: add constraints
        inp = self.__prep(job, config, "opt", no_solv=config.gas_phase)

        # write input into file "{filename}.inp" in a subdir created for the
        # conformer
        with open(inputpath, "w") as f:
            f.writelines(inp)

        # Get gbw files for initial guess
        if (
            config.copy_mo
            and job.mo_guess is not None
            and not isinstance(job.mo_guess, tuple)
        ):
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # call orca
        call = [config.paths.orca, f"{filename}.inp"]
        returncode, _ = self._make_call(call, outputpath, jobdir)
        # NOTE: using orca returncodes it is not possible to determine wether the calculation converged

        meta.success = returncode == 0

        # read output
        with open(outputpath) as out:
            lines = out.readlines()

        # Check for errors in the output file in case returncode is 0
        if meta.success:
            meta.error = self.__check_output(lines)
            if meta.error != "":
                return result, meta
        else:
            meta.error = self.__returncode_to_err.get(returncode, "unknown_error")
            return result, meta

        # Check convergence
        if (
            next((True for x in lines if "OPTIMIZATION HAS CONVERGED" in x), None)
            is True
        ):
            result.converged = True
        else:
            result.converged = False

        # Get the number of cycles
        for line in lines:
            if "GEOMETRY OPTIMIZATION CYCLE" in line:
                result.cycles = int(line.split()[4])

        # Get energies for each cycle
        result.ecyc = [
            float(line.split("....")[-1].split()[0])
            for line in filter(lambda x: "Current Energy" in x, lines)
        ]

        # Get all gradient norms for evaluation
        result.gncyc = [
            float(line.split("....")[-1].split()[0])
            for line in filter(lambda x: "Current gradient norm" in x, lines)
        ]

        # Get the last gradient norm and energy
        result.energy = result.ecyc[-1]
        result.grad_norm = result.gncyc[-1]
        meta.success = True

        # Read out optimized geometry and update conformer geometry with this
        job.conf.fromxyz(os.path.join(jobdir, f"{filename}.xyz"))
        result.geom = job.conf.xyz

        if config.copy_mo:
            # store the path to the current .gbw file for this conformer if
            # possible
            if os.path.isfile(os.path.join(jobdir, f"{filename}.gbw")):
                result.mo_path = os.path.join(jobdir, f"{filename}.gbw")

        return result, meta

    @override
    def xtb_opt(
        self,
        job: JobContext,
        config: XTBOptJobConfig,
    ) -> tuple[OptResult, MetaData]:
        """
        Geometry optimization using ANCOPT and ORCA gradients.
        Note that solvation is handled here always implicitly.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :param config: XTB optimization configuration
        :return: Tuple of (OptResult, MetaData)
        """
        # NOTE: some "intuitivity problems":
        # the geometry of the conformer is written into a coord file and also into an xyz-file to be used by orca
        # xtb then outputs a file with the optimized geometry as 'xtbopt.coord', which is then read into the conformer
        # to update it's geometry

        # prepare result
        # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
        # 'gncyc' contains the gradient norms for all cycles
        # 'energy' contains the final energy of the optimization (converged or unconverged)
        # 'geom' stores the optimized geometry in GeometryData.xyz format
        result = OptResult()
        meta = MetaData(job.conf.name)
        xcontrolname = "xtb_opt-xcontrol-inp"

        files = [
            "xtbrestart",
            "xtbtopo.mol",
            xcontrolname,
            "wbo",
            "charges",
            "gfnff_topo",
        ]

        jobdir = self._setup(job, "xtb_opt")
        filename = "xtb_opt"

        # remove potentially preexisting files to avoid confusion
        for file in files:
            if os.path.isfile(os.path.join(jobdir, file)):
                os.remove(os.path.join(jobdir, file))

        # write conformer geometry to coord file
        with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as f:
            f.writelines(job.conf.tocoord())

        # write xyz-file for orca
        with open(os.path.join(jobdir, f"{filename}.xyz"), "w", newline=None) as f:
            f.writelines(job.conf.toxyz())

        # set orca input path
        inputpath = os.path.join(jobdir, f"{filename}.inp")

        # prepare input dict
        inp = self.__prep(
            job, config, "xtb_opt", no_solv=config.gas_phase, xyzfile=f"{filename}.xyz"
        )

        # write orca input into file "xtb_opt.inp" in a subdir created for the
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
                    f"   orca bin= {config.paths.orca}\n",
                    "$end\n",
                ]
            )

        # prepare configuration file for ancopt (xcontrol file)
        with open(os.path.join(jobdir, xcontrolname), "w", newline=None) as out:
            out.write("$opt \n")
            if config.macrocycles:
                out.write(f"maxcycle={config.optcycles} \n")
                out.write(f"microcycle={config.optcycles} \n")

            out.writelines(
                [
                    "average conv=true \n",
                    f"hlow={config.hlow} \n",
                    "s6=30.00 \n",
                    "engine=lbfgs\n",
                    "$external\n",
                    f"   orca input file= {filename}.inp\n",
                    f"   orca bin= {config.paths.orca} \n",
                ]
            )

            # Import constraints
            if config.constraints:
                out.writelines(config.constraints)

            out.write("$end \n")

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if (
            config.copy_mo
            and job.mo_guess is not None
            and not isinstance(job.mo_guess, tuple)
        ):
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # prepare xtb call
        call = [
            config.paths.xtb,
            f"{filename}.coord",  # name of the coord file generated above
            "--opt",
            config.optlevel,
            "--orca",
            "-I",
            xcontrolname,
        ]

        # set path to the ancopt output file
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # call xtb
        returncode, _ = self._make_call(call, outputpath, jobdir)

        # check if optimization finished without errors
        # NOTE: right now, not converging scfs are not handled because returncodes need to be implemented first
        if returncode != 0:
            meta.success = False
            meta.error = "unknown_error"
            # TODO: the xtb returncodes should be handled
            return result, meta

        # read output
        with open(outputpath) as f:
            lines = f.readlines()

        result.ecyc = []
        result.cycles = 0

        # Substrings indicating error in xtb
        error_ind = [
            "external code error",
            "|grad| > 500, something is totally wrong!",
            "abnormal termination of xtb",
        ]

        # Check if xtb terminated normally (if there are any error indicators
        # in the output)
        meta.success = (
            False
            if next((x for x in lines if any(y in x for y in error_ind)), None)
            is not None
            else True
        )
        if not meta.success:
            meta.error = "unknown_error"
            return result, meta

        # check convergence
        if (
            next((True for x in lines if "GEOMETRY OPTIMIZATION CONVERGED" in x), None)
            is True
        ):
            result.converged = True
        elif (
            next((True for x in lines if "FAILED TO CONVERGE GEOMETRY" in x), None)
            is True
        ):
            result.converged = False

        # Get the number of cycles
        # tmp is one of the values from the dict defined below
        tmp_val = {
            True: ("GEOMETRY OPTIMIZATION CONVERGED", 5),
            False: ("FAILED TO CONVERGE GEOMETRY", 7),
        }[result.converged]

        result.cycles = int(
            next(x for x in lines if tmp_val[0] in x).split()[tmp_val[1]]
        )

        # Get energies for each cycle
        result.ecyc.extend(
            float(line.split("->")[-1])
            for line in filter(lambda x: "av. E: " in x, lines)
        )

        # Get all gradient norms for evaluation
        result.gncyc = [
            float(line.split()[3])
            for line in filter(lambda x: " gradient norm " in x, lines)
        ]

        # Get the last gradient norm
        result.grad_norm = result.gncyc[-1]

        # store the final energy of the optimization in 'energy'
        result.energy = result.ecyc[-1]
        meta.success = True

        if config.copy_mo:
            # store the path to the current .gbw file for this conformer
            result.mo_path = os.path.join(jobdir, f"{filename}.gbw")

        # read out optimized geometry and update conformer geometry with this
        job.conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
        result.geom = job.conf.xyz

        return result, meta

    @override
    def nmr(
        self,
        job: JobContext,
        config: NMRJobConfig,
    ) -> tuple[NMRResult, MetaData]:
        """
        Calculate the NMR shieldings and/or couplings for a conformer using ORCA. ORCA gives only the active cores in the output
        so there is not need for more thinking here.
        Formatting:
            'shielding' contains a list of tuples (atom_index, shielding), with atom_index being the index of the atom
            in the internal coordinates of the GeometryData.
            'couplings' contains a list of tuples ((atom_index1, atom_index2), coupling), with the indices of the atoms
            in the internal coordinates of the GeometryData. A set is used to represent an atom pair and then converted
            to tuple to be serializable.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :param config: NMR configuration
        :return: Tuple of (NMRResult, MetaData)
        """
        # Set results
        result = NMRResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "nmr")
        filename = "nmr"

        # Set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare an input file for _sp
        inp = self.__prep(job, config, "nmr", no_solv=config.gas_phase)
        with open(inputpath, "w") as f:
            f.writelines(inp)

        # Run _sp using the input file generated above
        _, spmeta = self._sp(job, config, jobdir=jobdir, filename=filename, prep=False)

        if not spmeta.success:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # Grab shieldings and energy from the output
        with open(outputpath) as f:
            lines = f.readlines()

        # Get final energy
        # result.energy = spres.energy

        # For shieldings watch out for the line "CHEMICAL SHIELDING SUMMARY
        # (ppm)"
        start = (
            lines.index(next(x for x in lines if "CHEMICAL SHIELDING SUMMARY" in x)) + 6
        )

        result.shieldings = []

        for line in lines[start:]:
            # Try to extract values from the lines, if that fails
            # (probably IndexError) stop
            try:
                result.shieldings.append((int(line.split()[0]), float(line.split()[2])))
            except IndexError:
                break

        # Sort shieldings by atom index
        result.shieldings.sort(key=lambda x: x[0])

        with open(outputpath) as f:
            lines = f.readlines()

        start = lines.index(next(x for x in lines if "SPIN-SPIN COUPLING" in x)) + 6

        lines = lines[start:]

        end = lines.index(next(x for x in lines if "SUMMARY OF ISOTROPIC" in x)) - 3

        lines = lines[:end]

        couplings: list[tuple[frozenset[int], float]] = []

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
                    if "(Total)" in line2:
                        coupling = float(line2.split()[5])
                        couplings.append((pair, coupling))
                        break

        # Convert to set and back to get rid of duplicates
        # ('vectorizing the symmetric matrix')
        couplings = list(set(couplings))

        # Convert all the frozensets to a tuple to be serializable
        for i in range(len(couplings)):
            pair_tuple = cast(tuple[int, int], tuple(couplings[i][0]))
            coupling = couplings[i][1]
            result.couplings.append(
                (
                    pair_tuple,
                    couplings[i][1],
                )
            )

        # Sort couplings by pairs
        result.couplings.sort(key=lambda x: x[0])

        meta.success = True

        return result, meta

    @override
    def uvvis(
        self,
        job: JobContext,
        config: UVVisJobConfig,
    ) -> tuple[UVVisResult, MetaData]:
        """
        Run a single-point to calculate the oscillator strengths and excitation wavelengths.

        :param job: JobContext object containing the job information, metadata is stored in job.meta
        :param config: UV-Vis configuration
        :return: Tuple of (UVVisResult, MetaData)
        """
        # Set results
        result = UVVisResult()
        meta = MetaData(job.conf.name)

        jobdir = self._setup(job, "uvvis")
        filename = "uvvis"

        # Set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare an input file for _sp
        inp = self.__prep(job, config, "uvvis", no_solv=config.gas_phase)
        with open(inputpath, "w") as f:
            f.writelines(inp)

        # Run _sp using the input file generated above
        _, spmeta = self._sp(
            job, config, jobdir=jobdir, filename=f"{filename}", prep=False
        )

        if not spmeta.success:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # Grab excitations and energy from the output
        with open(outputpath) as f:
            lines = f.readlines()

        # Get final energy
        # result.energy = spres.energy

        # Find the line index of the actual excitations information
        start = lines.index(
            next(filter(lambda line: "ABSORPTION SPECTRUM" in line, lines))
        )

        # Get all the lines where the excitations are listed (+5 because of spacing in output file)
        uvvis_table = lines[start + 5 : start + 5 + config.nroots]

        # Extract excitation wavelengths and oscillator strengths
        result.excitations = []
        for row in uvvis_table:
            spl = row.split()
            result.excitations.append(
                {"wavelength": float(spl[-6]), "osc_str": float(spl[-5])}
            )

        meta.success = True

        return result, meta

    @override
    def rot(self, *args, **kwargs):
        raise NotImplementedError


Factory.register_builder(Prog.ORCA, OrcaProc)
