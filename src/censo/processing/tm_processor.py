"""
Contains TmProc class for calculating TURBOMOLE related properties of conformers.
"""

import os
import shutil
import math
from pathlib import Path
from typing import cast

from .qm_processor import QmProc
from ..logging import setup_logger
from ..parallel import (
    GsolvResult,
    MetaData,
    NMRResult,
    OptResult,
    SPResult,
    ParallelJob,
)
from ..params import WARNLEN, R, AU2KCAL, TmSolvMod, ASSETS_PATH, Prog
from ..utilities import frange, Factory
from ..config.job_config import NMRJobConfig, SPJobConfig, XTBOptJobConfig
from ..assets import FUNCTIONALS, SOLVENTS

logger = setup_logger(__name__)


class TmProc(QmProc):
    """
    Performs calculations using TURBOMOLE.
    """

    progname = Prog.TM

    __gridsettings = {
        "low": ["    gridsize m3", "$scfconv 6"],
        "low+": ["    gridsize m4", "$scfconv 6"],
        "high": ["    gridsize m4", "$scfconv 7"],
        "high+": ["    gridsize m5", "$scfconv 7"],
        "nmr": ["    gridsize 5", "$scfconv 7"],
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
        "woctanol": 8.1,
    }

    # Contains mapping from lower case to proper spelling of bases for TM
    # (should be available for most elements)
    # TODO: consider adding custom basis mapping in user assets
    __basis_mapping = {
        "dzp": "DZP",
        "svp": "SVP",
        "def-svp": "def-SVP",
        "def2-svp": "def2-SVP",
        "dhf-svp": "dhf-SVP",
        "dhf-svp-2c": "dhf-SVP-2c",
        "tz": "TZ",
        "tzv": "TZV",
        "tzp": "TZP",
        "tzvp": "TZVP",
        "def-tzvp": "def-TZVP",
        "def2-tzvp": "def2-TZVP",
        "dhf-tzvp": "dhf-TZVP",
        "dhf-tzvp-2c": "dhf-TZVP-2c",
        "tzvpp": "TZVPP",
        "def-tzvpp": "def-TZVPP",
        "def2-tzvpp": "def2-TZVPP",
        "def2-mtzvpp": "def2-mTZVPP",
        "dhf-tzvpp": "dhf-TZVPP",
        "dhf-tzvpp-2c": "dhf-TZVPP-2c",
        "tzvppp": "TZVPPP",
        "def-qzv": "def-QZV",
        "def2-qzv": "def2-QZV",
        "qzv": "QZV",
        "qz": "QZ",
        "def-qzvp": "def-QZVP",
        "def2-qzvp": "def2-QZVP",
        "dhf-qzvp": "dhf-QZVP",
        "dhf-qzvp-2c": "dhf-QZVP-2c",
        "qzvp": "QZVP",
        "qzp": "QZP",
        "def-qzvpp": "def-QZVPP",
        "def2-qzvpp": "def2-QZVPP",
        "dhf-qzvpp": "dhf-QZVPP",
        "dhf-qzvpp-2c": "dhf-QZVPP-2c",
        "qzvpp": "QZVPP",
        "qzpp": "QZPP",
        "minix": "minix",
        "sto-3g": "sto-3g",
        "3-21g": "3-21g",
        "4-31g": "4-31g",
        "def2-qzvppd": "def2-QZVPPD",
        "def2-qzvpd": "def2-QZVPD",
        "6-31g": "6-31G",
        "6-31g*": "6-31G*",
        "6-31g**": "6-31G**",
        "6-311g": "6-311G",
        "6-311g*": "6-311G*",
        "6-311g**": "6-311G**",
        "6-311++g**": "6-311++G**",
        "6-311g(2df,2pd)": "6-311G(2df,2pd)",
        "cc-pvdz": "cc-pVDZ",
        "aug-cc-pvdz": "aug-cc-pVDZ",
        "yp-aug-cc-pvdz": "YP-aug-cc-pVDZ",
        "d-aug-cc-pvdz": "d-aug-cc-pVDZ",
        "cc-pvtz": "cc-pVTZ",
        "cc-pvtz-kernel": "cc-pVTZ-kernel",
        "aug-cc-pvtz": "aug-cc-pVTZ",
        "yp-aug-cc-pvtz": "YP-aug-cc-pVTZ",
        "d-aug-cc-pvtz": "d-aug-cc-pVTZ",
        "d-aug-cc-pvtz-oep": "d-aug-cc-pVTZ-oep",
        "cc-pvqz": "cc-pVQZ",
        "aug-cc-pvqz": "aug-cc-pVQZ",
        "yp-aug-cc-pvqz": "YP-aug-cc-pVQZ",
        "d-aug-cc-pvqz": "d-aug-cc-pVQZ",
        "cc-pv5z": "cc-pV5Z",
        "aug-cc-pv5z": "aug-cc-pV5Z",
        "yp-aug-cc-pv5z": "YP-aug-cc-pV5Z",
        "d-aug-cc-pv5z": "d-aug-cc-pV5Z",
        "cc-pv6z": "cc-pV6Z",
        "aug-cc-pv6z": "aug-cc-pV6Z",
        "r12": "r12",
        "cc-pvdz-f12": "cc-pVDZ-F12",
        "cc-pvtz-f12": "cc-pVTZ-F12",
        "cc-pvqz-f12": "cc-pVQZ-F12",
        "def2-svpd": "def2-SVPD",
        "def2-tzvppd": "def2-TZVPPD",
        "def2-tzvpd": "def2-TZVPD",
        "dz": "DZ",
        "sv": "SV",
        "sv(p)": "SV(P)",
        "def-sv(p)": "def-SV(P)",
        "def2-sv(p)": "def2-SV(P)",
        "dhf-sv(p)": "dhf-SV(P)",
        "dhf-sv(p)-2c": "dhf-SV(P)-2c",
    }

    def __prep(
        self,
        job: ParallelJob,
        config: SPJobConfig,
        jobtype: str,
        jobdir: str,
        no_solv: bool = False,
    ) -> None:
        """
        Prepares TURBOMOLE input files for a specified jobtype.
        """
        inp = []

        func = config.func
        func_name = FUNCTIONALS[func][Prog.TM.value]
        func_type = FUNCTIONALS[func]["type"]

        try:
            basis = self.__basis_mapping[config.basis.lower()]
        except KeyError as exc:
            raise KeyError(
                f"Basis {config.basis} could not be found for TURBOMOLE input preparation. "
                f"Available basis sets: {list(self.__basis_mapping.values())}"
            ) from exc

        inp.extend(["$atoms", f"    basis={basis}"])

        inp.extend(["$dft", f"   functional {func_name}"])

        # Configure grid
        inp.extend(self.__gridsettings[config.grid])

        # r2scan-3c should use m4 grid and radsize 10
        if func == "r2scan-3c":
            if "m3" in inp:
                inp[inp.index("m3")] = "m4"
            inp.insert(inp.index("$dft") + 2, "    radsize 10")

        # Add dispersion
        # dispersion correction information
        # TODO: temporary solution (not very nice)
        mapping = {
            "d3bj": "$disp3 -bj",
            "d3(0)": "$disp3 -zero",
            "d4": "$disp4",
            "nl": "$donl",
        }

        disp = FUNCTIONALS[func]["disp"]
        if disp not in ["composite", "novdw", "included"]:
            inp.append(mapping[disp])
        elif func.endswith("-v"):
            inp.append("$doscnl")

        inp.append("$rij")

        # Add charge and unpaired info
        inp.append(f"$charge={job.charge} unpaired={job.unpaired}")

        # Write coord file
        with open(os.path.join(jobdir, "coord"), "w") as f:
            f.writelines(job.conf.tocoord())

        inp.append("$coord file=coord")

        # Do further manipulations of input files
        self.__prep_main(
            inp,
            func,
            disp,
            func_type,
            basis,
        )
        if not no_solv:
            self.__prep_solv(inp, config, jobtype)

        self.__prep_special(inp, config, jobtype)

        # TODO: add template

        inp.append("$end")
        inp = [i if i.endswith("\n") else i + "\n" for i in inp]
        (Path(jobdir) / "control").write_text("".join(inp))

    def __prep_main(
        self,
        inp: list[str],
        func: str,
        disp: str,
        func_type: str,
        basis: str,
    ):
        # Special treatment for KT1/KT2
        if "kt" in func:
            func_line_index = next(inp.index(l) for l in inp if "functional" in l)
            inp[func_line_index] = "   functional xcfun set-gga\n"
            inp.insert(func_line_index + 1, f"   functional xcfun kt{func[2]} 1.0\n")
        # Special treatment for b97-3c
        elif func == "b97-3c":
            # Needs three-body dispersion
            disp_line_index = next(inp.index(l) for l in inp if "disp" in l)
            inp[disp_line_index] = "$disp3 -bj -abc\n"

        # Enable non local dispersion
        if disp == "nl":
            inp.append("$donl\n")

        # Handle GCP
        if func_type != "composite":
            gcp_keywords = {
                "minis": "MINIS",
                "sv": "SV",
                "6-31g(d)": "631GD",
                "def2-sv(p)": "SV(P)",
                "def2-svp": "SVP",
                "def2-tzvp": "TZ",
            }
            if basis.lower() in gcp_keywords:
                if basis.lower() == "def2-sv(p)":
                    inp.append("$gcp dft/sv(p)\n")
                else:
                    inp.append(f"$gcp dft/{basis.lower().replace('-', '')}\n")

    def __prep_solv(self, inp: list[str], config: SPJobConfig, jobtype: str):
        inp.insert(-1, "$cosmo\n")

        assert config.solvent
        assert config.sm

        sm = config.sm
        solv_key = SOLVENTS[config.solvent][sm]

        # write DC in any case
        inp.append(f" epsilon= {self.__cosmo_dcs[solv_key]}\n")

        if sm == TmSolvMod.DCOSMORS:
            # if using dcosmors also add the potential file path
            if solv_key not in [
                "woctanol",
                "hexadecane",
                "octanol",
            ]:
                inp.append(
                    f"$dcosmo_rs file={solv_key}_25.pot\n",
                )
            else:
                # The three solvents above are specifically defined in the assets
                # TODO: this opens the possibility to insert your own potential files
                inp.append(
                    f"$dcosmo_rs file={os.path.join(ASSETS_PATH, solv_key)}_25.pot\n",
                )

        if jobtype == "rot":
            inp.extend(
                [
                    " cavity closed\n",
                    " use_contcav\n",
                    " nspa=272\n",
                    " nsph=162\n",
                    "$cosmo_isorad\n",
                ]
            )

    def __prep_special(self, lines: list[str], config: SPJobConfig, jobtype: str):
        # Set NMR parameters
        if "nmr" in jobtype:
            assert isinstance(config, NMRJobConfig)
            # Determine the settings that need to be put into the input file for the NMR calculation
            # $nucsel does not work properly with capital letters
            todo = [element.lower() for element in config.active_nuclei.split(",")]

            lines.append("$rpacor 10000")

            lines.append("$ncoupling")

            if config.fc_only:
                lines.extend([" simple\n", " thr=0.0\n"])

            # nucsel only required if not all elements are active
            if not all(element in todo for element in ["h", "c", "f", "si", "p"]):
                todo = [f'"{e}"' for e in todo]
                lines.extend(
                    [
                        "$nucsel " + " ".join(todo) + "\n",
                        "$nucsel2 " + " ".join(todo) + "\n",
                    ]
                )

            lines.append("$rpaconv 8\n")
        elif jobtype == "rot":
            # TODO:
            """
            controlappend.append("$scfinstab dynpol nm")
            for i in self.job["freq_or"]:
                controlappend.append(f" {i}")  # e.g. 589
            controlappend.append("$velocity gauge")
            controlappend.append("$rpaconv 4")

            """
            raise NotImplementedError("Optical rotation not available yet!")

    @staticmethod
    def __copy_mo(
        jobdir: str, guess_file: str | Path | tuple[str | Path, str | Path]
    ) -> None:
        """
        Copy the MO file(s) for TURBOMOLE (should be TM format).
        """
        if isinstance(guess_file, tuple):
            # open shell guess
            if all(
                os.path.isfile(f)
                and any(g in str(f) for g in ["alpha", "beta"])
                and not any(
                    os.path.join(jobdir, g) == str(f) for g in ["alpha", "beta"]
                )
                for f in guess_file
            ):
                # All MO files found and not already in dir
                # Copy MO files
                for f in guess_file:
                    g = os.path.split(f)[1] if isinstance(f, str) else f.stem
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying {g} file from {f}."
                    )
                    shutil.copy(f, os.path.join(jobdir, g))
        else:
            # closed shell guess
            if (
                os.path.isfile(guess_file)
                and os.path.split(guess_file)[1] == "mos"
                and os.path.join(jobdir, "mos") != guess_file
            ):
                # Copy MO file
                logger.debug(
                    f"{f'worker{os.getpid()}:':{WARNLEN}}Copying mos file from {guess_file}."
                )
                shutil.copy(guess_file, os.path.join(jobdir, "mos"))

    @staticmethod
    def __check_output(lines: list[str]) -> str:
        """
        Checks the lines from the output file for errors and returns them.

        Args:
            lines: list of lines from the output file.

        Returns:
            str | None: error message if an error was found, None otherwise
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

    def _sp(
        self,
        job: ParallelJob,
        jobdir: str,
        config: SPJobConfig,
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[SPResult, MetaData]:
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
        """
        # set results
        result = SPResult()

        meta = MetaData(job.conf.name)

        # set in/out path
        outputpath = os.path.join(jobdir, "ridft.out")

        # check, if there is an existing mo/alpha,beta file and copy it if option
        # 'copy_mo' is true
        # mo files: mos/alpha,beta
        # NOTE: this HAS TO BE in this order, otherwise ridft fails to read mos
        if config.copy_mo and job.mo_guess is not None:
            self.__copy_mo(jobdir, job.mo_guess)

        if prep:
            self.__prep(job, config, "sp", jobdir, no_solv=config.gas_phase or no_solv)

        # call turbomole
        call = ["ridft"]
        returncode, errors = self._make_call(call, outputpath, jobdir)

        meta.success = returncode == 0
        if not meta.success:
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")

        lines = Path(outputpath).read_text().split("\n")

        # Check for errors in the output file in case returncode is 0
        if meta.success:
            meta.error = self.__check_output(lines)
            meta.success = meta.error == ""
        else:
            meta.error = self.__returncode_to_err.get(returncode, "unknown_error")

        # Get final energy
        try:
            result.energy = next(
                (
                    float(line.split()[4])
                    for line in lines
                    if "|  total energy      = " in line
                )
            )
        except StopIteration:
            meta.success = False
            meta.error = "Could not parse final energy"

        if config.copy_mo:
            # store the path to the current MO file(s) for this conformer if possible
            if os.path.isfile(os.path.join(jobdir, "mos")):
                result.mo_path = os.path.join(jobdir, "mos")
            elif os.path.isfile(os.path.join(jobdir, "alpha")):
                result.mo_path = (
                    os.path.join(jobdir, "alpha"),
                    os.path.join(jobdir, "beta"),
                )

        return result, meta

    @QmProc._run
    def sp(self, *args, **kwargs):
        return self._sp(*args, **kwargs)

    @QmProc._run
    def gsolv(
        self,
        job: ParallelJob,
        jobdir: str,
        config: SPJobConfig,
    ) -> tuple[GsolvResult, MetaData]:
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
        result = GsolvResult()
        meta = MetaData(job.conf.name)

        assert config.sm
        assert config.solvent
        # assert config.multitemp
        assert config.temperature
        # assert config.trange

        if config.sm in [TmSolvMod.COSMO, TmSolvMod.DCOSMORS]:
            # Non-COSMORS procedure:
            # Run gas-phase sp
            spres, spmeta = self._sp(job, jobdir, config, no_solv=True)

            if spmeta.success:
                result.energy_gas = spres.energy
            else:
                meta.success = False
                meta.error = spmeta.error
                return result, meta

            # Run solution sp
            spres, spmeta = self._sp(job, jobdir, config)

            if spmeta.success:
                result.energy_solv = spres.energy
            else:
                meta.success = False
                meta.error = spmeta.error
                return result, meta

            # Calculate gsolv from energy difference
            result.gsolv = result.energy_solv - result.energy_gas
            meta.success = True
        else:
            # COSMORS procedure:
            # Run gas-phase sp with unaltered settings
            spres, spmeta = self._sp(job, jobdir, config, no_solv=True)

            if spmeta.success:
                result.energy_gas = spres.energy
            else:
                meta.success = False
                meta.error = spmeta.error
                return result, meta

            # Run gas-phase sp with BP86 and def2-TZVP (normal)/def2-TZVPD (fine)
            # Turn off copying mos since this will lead to errors otherwise
            # (due to incorrect order of prep/copy_mos)
            spconfig = SPJobConfig(
                copy_mo=False,
                func="bp86",
                basis="def2-tzvp" if config.sm == TmSolvMod.COSMORS else "def2-tzvpd",
                grid=config.grid,
                template=config.template,
                gas_phase=True,
            )

            spres, spmeta = self._sp(job, jobdir, spconfig, no_solv=True)

            if not spmeta.success:
                meta.success = False
                meta.error = spmeta.error
                return result, meta

            # Run special cosmo sp with BP86 and def2-TZVP (normal)/def2-TZVPD (fine)
            # Write special settings for cosmo into control file
            with open(os.path.join(jobdir, "control"), "r+") as f:
                lines = f.readlines()

                lines[-1:-1] = [
                    "$cosmo\n",
                    " epsilon=infinity\n",
                    " use_contcav\n",
                    " cavity closed\n",
                    " nspa=272\n",
                    " nsph=162\n",
                    "$cosmo_out  file=out.cosmo\n",
                ]

                f.seek(0)
                f.writelines(lines)
                f.truncate()

            # Run sp
            spconfig.gas_phase = True
            spconfig.sm = config.sm
            spconfig.solvent = config.solvent
            spres, spmeta = self._sp(job, jobdir, spconfig, prep=False)

            if not spmeta.success:
                meta.success = False
                meta.error = spmeta.error
                return result, meta

            # Prepare cosmotherm.inp
            if config.sm == TmSolvMod.COSMORS:
                setup = (
                    self.paths["cosmorssetup"]
                    if not "FINE" in self.paths["cosmorssetup"]
                    else self.paths["cosmorssetup"].replace("TZVPD_FINE", "TZVP")
                )
            else:
                setup = (
                    self.paths["cosmorssetup"]
                    if "FINE" in self.paths["cosmorssetup"]
                    else self.paths["cosmorssetup"].replace("TZVP", "TZVPD_FINE")
                )
            lines = [
                f"ctd = {setup} cdir = {os.path.join(self.paths['cosmotherm'], 'CTDATA-FILES')}\n",
                "EFILE VPFILE\n",
                "!!\n",
            ]
            db = os.path.join(
                self.paths["cosmotherm"],
                "DATABASE-COSMO",
                (
                    "BP-TZVP-COSMO"
                    if config.sm == TmSolvMod.COSMORS
                    else "BP-TZVPD-FINE"
                ),
            )

            solv_key = SOLVENTS[config.solvent][config.sm]
            if solv_key == "woctanol":
                lines.extend(
                    [
                        f"f = h2o.cosmo fdir={db} autoc\n",
                        f"f = 1-octanol.cosmo fdir={db} autoc\n",
                    ]
                )
                mix = "0.27 0.73"
            else:
                lines.append(f"f = {solv_key}.cosmo" + f" fdir={db} autoc\n")
                mix = "1.0 0.0"

            lines.append("f = out.cosmo\n")

            # TODO: are these temperature dependent terms used somewhere?
            # if config.multitemp:
            #     trange = frange(
            #         config.trange[0],
            #         config.trange[1],
            #         step=config.trange[2],
            #     )
            #
            #     # Always append the fixed temperature to the trange so that it is the last value
            #     trange.append(config.temperature)
            #
            #     # Write trange to the xcontrol file
            #     for t in trange:
            #         lines.append(f"henry xh={{{mix}}} tc={t - 273.15} Gsolv\n")
            # else:
            lines.append(f"henry xh={{{mix}}} tc={config.temperature - 273.15} Gsolv\n")

            with open(os.path.join(jobdir, "cosmotherm.inp"), "w") as f:
                f.writelines(lines)

            # Run cosmotherm
            outputpath = os.path.join(jobdir, "cosmotherm.out")
            call = ["cosmotherm", "cosmotherm.inp"]
            returncode, errors = self._make_call(call, outputpath, jobdir)

            meta.success = returncode == 0
            if not meta.success:
                logger.warning(
                    f"Job for {job.conf.name} failed. Stderr output:\n{errors}"
                )
                meta.error = "unknown_error"
                return result, meta

            # Read output
            gsolvt = {}
            videal = (
                24.789561955 / 298.15
            )  # molar volume for ideal gas at 298.15 K 100.0 kPa

            cosmothermtab = os.path.join(jobdir, "cosmotherm.tab")
            with open(cosmothermtab, "r") as inp:
                lines = inp.readlines()
            for line in lines:
                if "T=" in line:
                    temp = float(line.split()[5])

                    # Add volume work
                    vwork = R * temp * math.log(videal * temp)
                elif " out " in line:
                    # Add volume work
                    gsolvt[temp] = float(line.split()[-1]) / AU2KCAL + vwork / AU2KCAL

            # result.gsolvt = gsolvt
            result.gsolv = gsolvt[config.temperature]
            result.energy_solv = result.energy_gas + result.gsolv

            # cosmothermd
            with open(os.path.join(jobdir, "cosmors.out"), "w") as out:
                temp = config.temperature
                vwork = R * temp * math.log(videal * temp)

                out.writelines(
                    [
                        "This is cosmothermrd (python version in ENSO) (SG,FB,SAW, 06/18)\n",
                        "final thermochemical solvation properties in kcal/mol\n"
                        "----------------------------------------------------------\n",
                        " Gsolv({} K)= {:10.3f}\n".format(
                            temp, result.gsolv * AU2KCAL - vwork
                        ),
                        " VWork({} K)= {:10.3f}\n".format(temp, vwork),
                        " Gsolv+VWork({} K)= {:10.3f}\n".format(
                            # volwork already included!
                            temp,
                            result.gsolv * AU2KCAL,
                        ),
                    ]
                )

        return result, meta

    @QmProc._run
    def xtb_opt(
        self,
        job: ParallelJob,
        jobdir: str,
        config: XTBOptJobConfig,
        filename: str = "xtb_opt",
    ) -> tuple[OptResult, MetaData]:
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
        """
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

        # remove potentially preexisting files to avoid confusion
        for file in files:
            if os.path.isfile(os.path.join(jobdir, file)):
                os.remove(os.path.join(jobdir, file))

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
                ]
            )

            # # Import constraints
            # if job.prepinfo["xtb_opt"]["constraints"] is not None:
            #     with open(job.prepinfo["xtb_opt"]["constraints"], "r") as f:
            #         lines = f.readlines()
            #
            #     out.writelines(lines)

            out.write("$end \n")

        # check, if there is an existing mo/alpha,beta file and copy it if option
        # 'copy_mo' is true
        # mo files: mos/alpha,beta
        if config.copy_mo and job.mo_guess is not None:
            self.__copy_mo(jobdir, job.mo_guess)

        self.__prep(job, config, "xtb_opt", jobdir)

        # prepare xtb call
        call = [
            self.paths[Prog.XTB.value],
            "coord",  # name of the coord file generated above
            "--opt",
            config.optlevel,
            "--tm",
            "-I",
            xcontrolname,
        ]

        # set path to the ancopt output file
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # call xtb
        returncode, errors = self._make_call(call, outputpath, jobdir)

        # check if optimization finished without errors
        # NOTE: right now, not converging scfs are not handled because returncodes need to be implemented first
        if returncode != 0:
            meta.success = False
            meta.error = "unknown_error"
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
            # TODO: the xtb returncodes should be handled
            return result, meta

        # read output
        with open(outputpath, "r") as file:
            lines = file.readlines()

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
        tmp = {
            True: ("GEOMETRY OPTIMIZATION CONVERGED", 5),
            False: ("FAILED TO CONVERGE GEOMETRY", 7),
        }
        tmp = tmp[result.converged]

        result.cycles = int(next(x for x in lines if tmp[0] in x).split()[tmp[1]])

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
            # store the path to the current MO file(s) for this conformer if possible
            if os.path.isfile(os.path.join(jobdir, "mos")):
                result.mo_path = os.path.join(jobdir, "mos")
            elif os.path.isfile(os.path.join(jobdir, "alpha")):
                result.mo_path = (
                    os.path.join(jobdir, "alpha"),
                    os.path.join(jobdir, "beta"),
                )

        # read out optimized geometry and update conformer geometry with this
        job.conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
        result.geom = job.conf.xyz

        return result, meta

    @QmProc._run
    def opt(self, *args, **kwargs):
        raise NotImplementedError(
            "Pure TURBOMOLE geometry optimization not available yet."
        )

    @QmProc._run
    def nmr(
        self,
        job: ParallelJob,
        jobdir: str,
        config: NMRJobConfig,
    ) -> tuple[NMRResult, MetaData]:
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
        result = NMRResult()
        meta = MetaData(job.conf.name)

        # Run sp first
        self.__prep(job, config, "nmr", jobdir)
        _, spmeta = self._sp(job, jobdir, config, prep=False)

        if not spmeta.success:
            meta.success = False
            meta.error = spmeta.error
            return result, meta

        # result.energy = spres.energy

        # Run shielding calculations if requested
        if config.shieldings:
            # Set in/out path
            outputpath = os.path.join(jobdir, "mpshift.out")

            call = ["mpshift", "-smpcpus", f"{job.omp}"]

            # Run mpshift for shielding calculation
            returncode, errors = self._make_call(call, outputpath, jobdir)

            meta.success = returncode == 0
            if not meta.success:
                logger.warning(
                    f"Job for {job.conf.name} failed. Stderr output:\n{errors}"
                )
                meta.error = "unknown_error"
                return result, meta

            # Grab shieldings from the output
            with open(outputpath, "r") as f:
                lines = f.readlines()

            try:
                start = lines.index(
                    next(x for x in lines if ">>>>> DFT MAGNETIC SHIELDINGS <<<<<" in x)
                )
            except StopIteration:
                meta.success = False
                meta.error = "Could not read shieldings"
                return result, meta

            lines = lines[start:]

            result.shieldings = []

            # Get lines with "ATOM" in it
            line_indices = [lines.index(l) for l in lines if "ATOM" in l]

            for i in line_indices:
                split = lines[i].split()
                result.shieldings.append((int(split[2]), float(split[4])))

            # Sort shieldings by atom index
            result.shieldings.sort(key=lambda x: x[0])

        # Run couplings calculation if requested
        if config.couplings:
            # Set in/out path
            outputpath = os.path.join(jobdir, "escf.out")

            call = ["escf", "-smpcpus", f"{job.omp}"]

            # Run escf for couplings calculation
            returncode, errors = self._make_call(call, outputpath, jobdir)

            meta.success = returncode == 0
            if not meta.success:
                logger.warning(
                    f"Job for {job.conf.name} failed. Stderr output:\n{errors}"
                )
                meta.error = "unknown_error"
                return result, meta

            # Grab couplings from the output
            with open(outputpath, "r") as f:
                lines = f.readlines()

            try:
                start = (
                    lines.index(
                        next(x for x in lines if "Nuclear coupling constants" in x)
                    )
                    + 3
                )
            except StopIteration:
                meta.success = False
                meta.error = "Could not read couplings"
                return result, meta

            lines = lines[start:]

            end = lines.index(next(x for x in lines if "-----" in x))

            lines = lines[:end]

            line_indices = [lines.index(l) for l in lines if len(l.split()) in [6, 7]]

            couplings: list[tuple[frozenset[int], float]] = []
            for i in line_indices:
                # pair needs to be a frozenset because normal sets are not hashable and can therefore not be part
                # of a normal set
                split = lines[i].split()
                pair = frozenset((int(split[1]), int(split[4].split(":")[0])))
                coupling = float(split[5])
                couplings.append((pair, coupling))

            # Convert to set and back to get rid of duplicates
            # ('vectorizing the symmetric matrix')
            couplings = list(set(couplings))

            # Convert all the frozensets to a tuple to be serializable
            for i in range(len(couplings)):
                pair = cast(tuple[int, int], tuple(couplings[i][0]))
                coupling = couplings[i][1]
                result.couplings.append(
                    (
                        pair,
                        couplings[i][1],
                    )
                )

            # Sort couplings by pairs
            result.couplings.sort(key=lambda x: x[0])

        meta.success = True

        return result, meta

    @QmProc._run
    def rot(self):
        raise NotImplementedError


Factory.register_builder(Prog.TM, TmProc)
