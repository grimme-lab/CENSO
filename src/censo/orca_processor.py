"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""
from .utilities import od_insert, setup_logger
import os
import shutil
from collections import OrderedDict
from functools import reduce

from .datastructure import GeometryData, ParallelJob
from .params import (
    CODING,
    USER_ASSETS_PATH,
    WARNLEN,
)
from .qm_processor import QmProc

logger = setup_logger(__name__)


class OrcaParser:
    """
    Parser for orca input files. Can read input files and transpile them to ordered dict. Also capable of writing an
    input file from a properly ordered dict (see __todict for format).
    """

    def read_input(self, path: str) -> OrderedDict:
        """
        Read orca input file at 'path' and parse into ordered dict.

        Args:
            path: path to the orca input file

        Returns:
            OrderedDict: OrderedDict containing the parsed input file
        """
        with open(path, "r") as infile:
            lines = infile.readlines()

        # Remove comment lines from input
        lines = self.__remove_comments(lines)

        # Convert lines to OrderedDict
        converted = self.__todict(lines)

        return converted

    def write_input(self, path: str, indict: OrderedDict) -> None:
        """
        Writes an orca input file at 'path' from an OrderedDict.

        Args:
            path: path to the output file
            indict: OrderedDict containing the input file

        Returns:
            None
        """
        with open(path, "w") as outfile:
            outfile.writelines(self.__tolines(indict))

    def __todict(self, lines: list[str]) -> OrderedDict:
        """
        Converts lines from ORCA input into OrderedDict

        Can also be used with an ORCA input template, the {main} and {postgeom} strings indicate the ordering.
        Of the settigs read from the template ({main} should come first, {postgeom} should come last) all settings
        below {postgeom} are put after the geometry definition block.

        Args:
            lines: list of lines from an ORCA input file

        Returns:
            OrderedDict: OrderedDict containing the parsed input file

        Format of the result:
        "main": [...], (contains all main input line options)

        "some setting": {"some option": [...], ...}, (contains a settings that is started with a % with the respective keywords and values)
        (the list under "some option" contains just all the following substrings, that, in the input file,
        are separated by a whitespace, nothing fancy)

        "geom": {"def": ["xyz/xyzfile", ...], "coord": [...]} (contains the geometry information block definition in 'def', 'coord' only if you use the option 'xyz')
        'coord' for e.g. H2: "coord": [["H", "0.0", "0.0", "0.0"], ["H", "0.0", "0.0", "0.7"]]

        "some setting after geom": {"some option": [...], ...}

        Comments get removed from the lines and therefore lost.

        TODO - cannot deal with %coord or internal cinpoordinates
        TODO - cannot deal with settings blocks with more 'end's than '%'s
        """

        converted = OrderedDict()

        for i, line in enumerate(lines):
            # strip leading whitespace
            line = line.lstrip()

            # main input line found
            if line.startswith("!"):
                # get rid of '!'
                line = line[1:]

                # orca input may have multiple lines beginning with '!'
                # account for that and read options from line
                if converted.get("main", False):
                    converted["main"].extend(line.split())
                else:
                    converted["main"] = line.split()

            # configuration line found
            elif line.startswith("%"):
                # get rid of '%'
                line = line[1:]
                split = line.split()

                # new configuration entry
                # keep in mind that converted remembers the order
                # in which key/value pairs are inserted
                setting = split[0]
                converted[setting] = {}

                # check it there is already some option in the same line as the
                # setting declaration
                if len(split) > 1:
                    option = split[1]
                    converted[setting][option] = split[2:]

                # find end of definition block
                end = i + self.__eob(lines[i:])

                # in case the eob is found within the line itself remove the
                # 'end' substring
                if end == i:
                    converted[setting][option].remove("end")
                # consume remaining definitions
                else:
                    for line2 in lines[i + 1:end]:
                        split = line2.split()
                        option = split[0]
                        converted[setting][option] = split[1:]
            # geometry input line found
            elif line.startswith("*") and "geom" not in converted.keys():
                if "geom" in converted.keys():
                    raise RuntimeError(
                        "Error parsing ORCA input file (double geom definition)."
                    )

                converted["geom"] = {}

                if "xyzfile" in line:
                    converted["geom"]["def"] = ["xyzfile"]
                # the 'xyz' keyword should be one of the first two substrings
                elif "xyz" in line.split()[0] or "xyz" in line.split()[1]:
                    converted["geom"]["def"] = ["xyz"]
                else:
                    raise RuntimeError("Error parsing ORCA input file.")

                # add the remaining line to the dict
                # get rid of '*'
                line = line[1:]
                converted["geom"]["def"].extend(line.split()[1:])

                # consume the remaining geometry information
                if converted["geom"]["def"][0] == "xyz":
                    converted["geom"]["coord"] = []
                    # find end of definition block
                    # start search from next line since geometry definition
                    # starts with an '*'
                    end = i + self.__eob(lines[i + 1:], endchar="*") + 1

                    for line2 in lines[i + 1: end]:
                        converted["geom"]["coord"].append(line2.split())
            # check for template lines
            # NOTE: only these two need to be checked since they're the only
            # ones that are sensitive to ordering
            elif "{main}" in line:
                # create the main key for ordering if it does not already exist
                if "main" not in converted.keys():
                    converted["main"] = []

                # main needs to be defined before everything
                if "main" != list(converted.keys())[0]:
                    raise RuntimeError(
                        "Error parsing ORCA input file (main not defined first)."
                    )
            elif "{postgeom}" in line:
                # there should be no geometry definition block in a template
                # file
                if "geom" in converted.keys():
                    raise RuntimeError(
                        "Error parsing ORCA input file (double geom definition)."
                    )
                # also main needs to be defined before geom
                elif "main" not in converted.keys():
                    raise RuntimeError(
                        "Error parsing ORCA input file (missing main definition)."
                    )

                # if we reach this point, the geom key should be created
                converted["geom"] = {}

        return converted

    @staticmethod
    def __eob(lines: list[str], endchar: str = "end") -> int:
        """
        find end of definition block
        """
        end = None
        for index, line in enumerate(lines):
            if endchar in line or ("%" in line and index > 0):
                end = index
                break

        if end is None:
            raise RuntimeError("Error parsing ORCA input file.")
        else:
            return end

    @staticmethod
    def __remove_comments(inlist: list[str]) -> list[str]:
        """
        remove all comments from an orca input
        """
        for i in range(len(inlist)):
            if "#" in inlist[i]:
                index = inlist[i].index("#")
                inlist[i] = inlist[i][:index]

        return inlist

    @staticmethod
    def __tolines(indict: OrderedDict) -> list[str]:
        """
        convert ordered dict to lines for orca input
        """
        # write lines into this list (include \n at the end of each line)
        lines = []

        # write main input line
        # start with an '!'
        lines.append("! ")

        # then append all keywords with whitespace inbetween
        lines[0] += reduce(lambda x, y: x + " " + y, indict["main"])
        lines[0] += "\n"

        # next, write all keywords and options that come between the main input
        # line and the geom input
        allkeys = list(indict.keys())

        # skip first key ('main')
        for key in allkeys[1: allkeys.index("geom")]:
            lines.append(f"%{key}\n")
            # FIXME - temporary workaround for definition blocks that have no
            # 'end', this code smells immensely
            try:
                for option in indict[key].keys():
                    lines.append(
                        f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n"
                    )
                lines.append("end\n")
            except TypeError:
                lines[-1] = f"%{key} {list(indict[key].keys())[0]}\n"

        # next, write the geometry input lines
        # geometry definition line (e.g. "* xyzfile 0 1 input.xyz" / "* xyz 0
        # 1")
        lines.append(
            f"* {reduce(lambda x, y: f'{x} {y}', indict['geom']['def'])}\n")

        # write coordinates if "xyz" keyword is used
        # if "xyzfile" is used, nothing more has to be done
        if indict["geom"].get("coord", False):
            for coord in indict["geom"]["coord"]:
                lines.append(f"{reduce(lambda x, y: f'{x} {y}', coord)}\n")
            lines.append("*\n")

        # lastly, write all the keywords and options that should be placed
        # after the geometry input (e.g. NMR settings)
        for key in allkeys[allkeys.index("geom") + 1:]:
            lines.append(f"%{key}\n")
            for option in indict[key].keys():
                lines.append(
                    f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n"
                )
            lines.append("end\n")

        return lines


class OrcaProc(QmProc):
    """
    Performs calculations with ORCA.
    """

    # contains grid settings for ORCA 5.0+ (True) and older versions (False)
    # can be chosen by simple keyword (low/low+/high/high+)
    __gridsettings = {
        False: {
            "low": ["grid4", "nofinalgrid", "loosescf"],
            "low+": ["grid4", "nofinalgrid", "scfconv6"],
            "high": ["grid4", "nofinalgrid", "scfconv7"],
            "high+": ["grid5", "nofinalgrid", "scfconv7"],
        },
        True: {
            "low": ["DEFGRID2", "loosescf"],
            "low+": ["DEFGRID2", "scfconv6"],
            "high": ["DEFGRID2", "scfconv7"],
            "high+": ["DEFGRID2", "scfconv7"],
        },
    }

    __req_settings = {
        **{
            "sp": [
                "func",
                "basis",
                "grid",
                "template",
                "func_name",
                "func_type",
                "disp",
                "gcp",
                # in principal there is the soft requirement of having sm and solvent_key_prog, but this is only used
                # in the case of implicit solvation (this is handled in the _sp method)
            ],
            "gsolv": [
                "sm",
                "solvent_key_prog",
            ],
            "xtb_opt": [
                "macrocycles",
                "optcycles",
                "hlow",
                "optlevel",
            ],
            "nmr": [
                "shieldings",
                "couplings",
                "h_active",
                "c_active",
                "f_active",
                "si_active",
                "p_active",
            ],
        },
        **QmProc._req_settings_xtb
    }

    @classmethod
    def check_requirements(cls, jobs: list[ParallelJob]) -> None:
        """
        Check, if the required settings are implemented in the jobs' prepinfo attributes for all the jobtypes.
        Checks requirements for single-point always, except for xtb-only jobtypes.

        Args:
            jobs(list[ParallelJob]): List of jobs that are to be checked.

        Returns:
            None
        """
        failed = False
        for job in jobs:
            for jt in job.jobtype:
                try:
                    # Check requirements for sp always except for xtb_sp or xtb_gsolv
                    if jt not in ["xtb_sp", "xtb_gsolv"]:
                        assert all(s in job.prepinfo[jt].keys()
                                   for s in cls.__req_settings["sp"])

                    # Check specific requirements
                    assert all(s in job.prepinfo[jt].keys()
                               for s in cls.__req_settings[jt])
                except AssertionError:
                    failed = True
                    logger.debug(
                        "The following settings are missing for implementation:")
                    if jt not in ["xtb_sp", "xtb_gsolv"]:
                        logger.debug(list(s for s in filter(
                            lambda x: x not in job.prepinfo[jt].keys(), cls.__req_settings["sp"])))

                    logger.debug(list(s for s in filter(
                        lambda x: x not in job.prepinfo[jt].keys(), cls.__req_settings[jt])))

        if failed:
            raise RuntimeError(
                "For at least one jobtype at least one setting is missing from job.prepinfo. "
                "Set __loglevel to logging.DEBUG and check log file for more info.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # expand jobtypes with special orca jobtypes
        self._jobtypes = {
            **self._jobtypes,
            **{
                "sp": self._sp,
                "gsolv": self._gsolv,
                "xtb_opt": self._xtb_opt,
                "nmr": self._nmr,
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

    def __prep(
            self,
            job: ParallelJob,
            jobtype: str,
            no_solv: bool = False,
            xyzfile: str = None) -> OrderedDict:
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
            OrderedDict: OrderedDict containing the input file

        NOTE: the xyzfile has to already exist, this function just bluntly writes the name into the input.

        TODO - OR/UVVis etc preparation steps
        """

        # check ORCA version
        orca5 = True if self._paths["orcaversion"].startswith("5") else False

        indict = OrderedDict()

        if self.instructions["template"]:
            # NOTE: when using templates we're not going to check for double definitions!
            # if the template is messed up, orca will fail and the user should deal with that
            # load template file
            try:
                indict = OrcaParser().read_input(
                    os.path.join(
                        USER_ASSETS_PATH,
                        f"{self.instructions['part_name']}.orca.template",
                    )
                )
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not find template file {self.instructions['part_name']}.orca.template."
                )

        # prepare the main line of the orca input
        indict = self.__prep_main(indict, jobtype, orca5)

        # prepare all options that are supposed to be placed before the
        # geometry definition
        indict = self.__prep_pregeom(indict, orca5, jobtype, no_solv, job.omp)

        # prepare the geometry
        indict = self.__prep_geom(indict, job.conf, xyzfile)

        # indict = self.__prep_postgeom(indict, jobtype, orca5)

        return indict

    def __prep_main(
        self, indict: OrderedDict, jobtype: str, orca5: bool
    ) -> OrderedDict:
        if "main" not in indict:
            indict["main"] = []

        # grab func, basis
        # first case catches the shielding & coupling all-in-one
        if jobtype == "nmr_s" or jobtype == "nmr":
            func = self.instructions["func_name_s"]
            basis = self.instructions["basis_s"]
            functype = self.instructions["func_type_s"]
            disp = self.instructions["disp_s"]
        elif jobtype == "nmr_j":
            func = self.instructions["func_name_j"]
            basis = self.instructions["basis_j"]
            functype = self.instructions["func_type_j"]
            disp = self.instructions["disp_j"]
        else:
            func = self.instructions["func_name"]
            basis = self.instructions["basis"]
            functype = self.instructions["func_type"]
            disp = self.instructions["disp"]

        indict["main"].append(func)

        # For non-composite methods, parameters for RI, RIJCOSX, and RIJK are
        # set
        if "composite" not in functype:
            indict["main"].append(basis)
            # set  RI def2/J,   RIJCOSX def2/J
            # this is only set if no composite DFA is used
            # settings for double hybrids
            if functype == "double":
                indict["main"].extend(["def2/J", "RIJCOSX"])

                if "nmr" in jobtype:
                    indict["main"].append("NOFROZENCORE")
                    indict = od_insert(
                        indict,
                        "mp2",
                        {"RI": ["true"], "density": ["relaxed"]},
                        list(indict.keys()).index("main") + 1,
                    )
                else:
                    indict["main"].append("frozencore")
                    indict = od_insert(
                        indict,
                        "mp2",
                        {"RI": ["true"]},
                        list(indict.keys()).index("main") + 1,
                    )

                def2cbasis = ("def2-svp", "def2-tzvp",
                              "def2-tzvpp", "def2-qzvpp")
                if basis.lower() in def2cbasis:
                    indict["main"].append(f"{basis}/C")
                    if not orca5:
                        indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])
                else:
                    indict["main"].append("def2-TZVPP/C")
                    if not orca5:
                        indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])

            # settings for hybrids
            elif "hybrid" in functype:
                indict["main"].append("RIJCOSX")
                if not orca5:
                    indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])

            # settings for (m)ggas
            elif "gga" in functype:
                indict["main"].append("RI")

        # use 'grid' setting from instructions to quickly configure the grid
        indict["main"].extend(self.__gridsettings[orca5]
                              [self.instructions["grid"]])

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
            indict["main"].append(mapping.get(disp, ""))

        if disp == "nl" and not orca5:
            indict["main"].append("vdwgrid3")

        # try to apply gcp if basis set available
        gcp_keywords = {
            "minis": "MINIS",
            "sv": "SV",
            "6-31g(d)": "631GD",
            "def2-sv(p)": "SV(P)",
            "def2-svp": "SVP",
            "def2-tzvp": "TZ",
        }
        if (
            self.instructions["gcp"]
            and basis.lower() in gcp_keywords.keys()
        ):
            indict["main"].append(
                f"GCP(DFT/{gcp_keywords[basis.lower()]})"
            )
        elif self.instructions["gcp"]:
            # TODO - error handling
            logger.warning(
                f"{f'worker{os.getpid()}:':{WARNLEN}}Selected basis not available for GCP."
            )

        # add job keyword for geometry optimizations
        # with ANCOPT
        if jobtype == "xtb_opt":
            indict["main"].append("ENGRAD")
        # for standard geometry optimization
        elif jobtype == "opt":
            indict["main"].append("OPT")

        return indict

    def __prep_pregeom(
            self, indict: OrderedDict, orca5: bool, jobtype: str, no_solv: bool, nprocs: int
    ) -> OrderedDict:
        # Check ORCA version (important for grid keywords)
        if orca5 and nprocs > 1:
            indict = od_insert(
                indict,
                "pal",
                {"nprocs": [nprocs]},
                list(indict.keys()).index("main") + 1,
            )

        # TODO - maybe limit TRAH macro steps (or when TRAH activates) to avoid
        # single jobs clogging everything up

        # set keywords for the selected solvent model
        # TODO - this is not good, reduce cyclomatic complexity
        if (
            not self.instructions["gas-phase"]
            and not no_solv
            and (
                "sm" in self.instructions.keys()
                or "sm_s" in self.instructions.keys()
            )
            or "sm_j" in self.instructions.keys()
        ):
            if jobtype == "nmr_s" or jobtype == "nmr":
                sm = self.instructions["sm_s"]
                solv_key = self.instructions["solvent_key_prog_s"]
            elif jobtype == "nmr_j":
                sm = self.instructions["sm_j"]
                solv_key = self.instructions["solvent_key_prog_j"]
            else:
                sm = self.instructions["sm"]
                solv_key = self.instructions["solvent_key_prog"]

            if sm == "smd":
                indict = od_insert(
                    indict,
                    "cpcm",
                    {
                        "smd": ["true"],
                        "smdsolvent": [f"\"{solv_key}\""],
                    },
                    list(indict.keys()).index("main") + 1,
                )
            elif sm == "cpcm":
                indict["main"].append(
                    f"CPCM({solv_key})")

        return indict

    def __prep_geom(
        self, indict: OrderedDict, conf: GeometryData, xyzfile: str
    ) -> OrderedDict:
        # unpaired, charge, and coordinates
        # by default coordinates are written directly into input file
        if xyzfile is None:
            indict["geom"] = {
                "def": [
                    "xyz",
                    self.instructions["charge"],
                    self.instructions["unpaired"] + 1,
                ],
                "coord": conf.toorca(),
            }
        else:
            indict["geom"] = {
                "def": [
                    "xyzfile",
                    self.instructions["charge"],
                    self.instructions["unpaired"] + 1,
                    xyzfile,
                ],
            }

        return indict

    def __prep_postgeom(
        self, indict: OrderedDict, jobtype: str, orca5: bool
    ) -> OrderedDict:
        # Set NMR parameters
        # TODO - this is not very nice, maybe make a list setting that contains
        # all the active nuclei
        todo = []
        if self.instructions["h_active"]:
            todo.append("H")
        if self.instructions["c_active"]:
            todo.append("C")
        if self.instructions["f_active"]:
            todo.append("F")
        if self.instructions["si_active"]:
            todo.append("Si")
        if self.instructions["p_active"]:
            todo.append("P")

        todo2 = []
        if self.instructions["shieldings"]:
            todo2.append("shift")
        if self.instructions["couplings"]:
            todo2.append("ssfc")

        if jobtype == "nmr":
            indict = od_insert(
                indict,
                "eprnmr",
                {
                    "Nuclei": [
                        "=",
                        ",".join(element for element in todo),
                        "{",
                        ",".join(x for x in todo2),
                        "}",
                    ],
                },
                list(indict.keys()).index("geom") + 1,
            )
            if self.instructions["shieldings"]:
                indict["eprnmr"]["origin"] = ["giao"]
                indict["eprnmr"]["giao_2el", "giao_2el_same_as_scf"]
                indict["eprnmr"]["giao_1el", "giao_1el_analytic"]

            if self.instructions["couplings"]:
                indict["eprnmr"]["SpinSpinRThresh"] = ["8.0"]

        return indict

    def _sp(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="sp",
        no_solv: bool = False,
        prep: bool = True,
    ) -> dict[str, float | None]:
        """
        ORCA single-point calculation.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file
            no_solv: if True, no solvent model is used
            prep: if True, a new input file is generated (you only really want to make use of this for NMR)

        Returns:
            dict[str, float | None]: dictionary containing the results of the calculation

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
        }
        job.meta.setdefault("sp", {})

        # set in/out path
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        indict = self.__prep(job, "sp", no_solv=no_solv)

        # check for flags raised for this jobtype
        if "sp" in job.flags or "gsolv" in job.flags:
            if (
                job.flags["sp"] == "scf_not_converged"
                or job.flags["gsolv"] == "scf_not_converged"
            ):
                indict = self.__apply_flags(indict, "scf_not_converged")

        # write input into file "{filename}.inp" in a subdir created for the
        # conformer
        parser = OrcaParser()
        parser.write_input(inputpath, indict)

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if self.copy_mo:
            if job.mo_guess is not None and os.path.isfile(job.mo_guess):
                if os.path.join(jobdir, f"{filename}.gbw") != job.mo_guess:
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {job.mo_guess}."
                    )
                    shutil.copy(job.mo_guess, os.path.join(
                        jobdir, f"{filename}.gbw"))

        # call orca
        call = [self._paths["orcapath"], f"{filename}.inp"]
        self._make_call(call, outputpath, jobdir)

        # do not check returncode, since orca doesn't give meaningful returncodes
        # check the outputfile instead, basically if the file doesn't say "SCF CONVERGED" or "SCF NOT CONVERGED"
        # somewhere, something went wrong
        # easily remedied:
        # "SCF NOT CONVERGED"
        # "Error encountered when trying to calculate the atomic fitting density!"

        # read output
        with open(outputpath, "r", encoding=CODING, newline=None) as out:
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

            # Check if scf is converged, this results to None if it cannot be
            # determined if the SCF converged or not
            meta["success"] = next(
                (True for line in lines if "SCF CONVERGED" in line), None) or next(
                (False for line in lines if "SCF NOT CONVERGED" in line), None)

        # TODO - this is not really correct, might not mean that the scf didn't
        # converge
        if meta["success"] is not None and not meta["success"]:
            meta["error"] = "SCF not converged"
        elif meta["success"] is None:
            meta["success"] = False
            meta["error"] = "Unknown error"

        if self.copy_mo:
            # store the path to the current .gbw file for this conformer if
            # possible
            if os.path.isfile(os.path.join(jobdir, f"{filename}.gbw")):
                job.meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # TODO - clean up?
        job.meta["sp"].update(meta)

        return result

    def _gsolv(self, job: ParallelJob, jobdir: str) -> dict[str, any]:
        """
        Calculates the solvation free enthalpy of a conformer using ORCA.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory

        Returns:
            dict[str, any | None]: dictionary containing the results of the calculation

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

        job.meta.setdefault("gsolv", {})
        meta = {
            "success": None,
            "error": None,
        }

        # calculate gas phase
        # TODO - does this need it's own folder? this is a bit messy
        spres = self._sp(job, jobdir, filename="sp_gas", no_solv=True)

        if job.meta["sp"]["success"]:
            result["energy_gas"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = "SCF not converged"
            job.meta["gsolv"].update(meta)
            return result

        # calculate in solution
        spres = self._sp(job, jobdir, filename="sp_solv")

        if job.meta["sp"]["success"]:
            result["energy_solv"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = "SCF not converged"
            job.meta["gsolv"].update(meta)
            return result

        # calculate solvation free enthalpy
        result["gsolv"] = result["energy_solv"] - result["energy_gas"]
        meta["success"] = True
        job.meta["gsolv"].update(meta)

        return result

    # TODO - split this up
    def _xtb_opt(self, job: ParallelJob, jobdir: str, filename: str = "xtb_opt") -> dict[str, any]:
        """
        ORCA geometry optimization using ANCOPT.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file

        Returns:
            dict[str, any]: dictionary containing the results of the calculation

        Keywords required in instructions:
        - optcycles
        - hlow
        - optlevel

        result = {
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "geom": None,
        }
        """
        # TODO - check if orca scf converges
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

        job.meta.setdefault("xtb_opt", {})
        meta = {
            "success": None,
            "error": None,
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
        parser = OrcaParser()
        indict = self.__prep(job, filename, xyzfile=f"{filename}.xyz")

        # write orca input into file "xtb_opt.inp" in a subdir created for the
        # conformer
        parser.write_input(inputpath, indict)

        # append some additional lines to the coord file for ancopt
        with open(os.path.join(jobdir, f"{filename}.coord"), "a", newline=None) as newcoord:
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
            if self.instructions["macrocycles"]:
                out.write(f"maxcycle={self.instructions['optcycles']} \n")
                out.write(f"microcycle={self.instructions['optcycles']} \n")

            out.writelines(
                [
                    "average conv=true \n",
                    f"hlow={self.instructions['hlow']} \n",
                    "s6=30.00 \n",
                    # remove unnecessary sp/gradient call in xTB
                    "engine=lbfgs\n",
                    "$external\n",
                    f"   orca input file= {filename}.inp\n",
                    f"   orca bin= {self._paths['orcapath']} \n",
                    "$end \n",
                ]
            )

        # check, if there is an existing .gbw file and copy it if option
        # 'copy_mo' is true
        if self.copy_mo:
            if job.mo_guess is not None and os.path.isfile(job.mo_guess):
                if os.path.join(jobdir, f"{filename}.gbw") != job.mo_guess:
                    logger.debug(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {job.mo_guess}."
                    )
                    shutil.copy(job.mo_guess, os.path.join(
                        jobdir, f"{filename}.gbw"))

        # prepare xtb call
        call = [
            self._paths["xtbpath"],
            f"{filename}.coord",  # name of the coord file generated above
            "--opt",
            self.instructions["optlevel"],
            "--orca",
            "-I",
            xcontrolname,
        ]

        # set path to the ancopt output file
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # call xtb
        returncode = self._make_call(call, outputpath, jobdir)

        # check if optimization finished without errors
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_opt"
            job.meta["xtb_opt"].update(meta)
            return result

        # read output
        with open(outputpath, "r", encoding=CODING, newline=None) as file:
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
            meta["error"] = "what went wrong in xtb_opt"
            job.meta["xtb_opt"].update(meta)
            return result

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
                next(x for x in lines if tmp[0] in x).split()[tmp[1]])

            # Get energies for each cycle
            result["ecyc"].extend(
                float(line.split("->")[-1])
                for line in filter(lambda x: "av. E: " in x, lines)
            )

            # Get the gradient norm (lines reversed so it takes the last
            # gradient norm value)
            result["grad_norm"] = float(
                next((x for x in lines[::-1] if " :: gradient norm      " in x), None).split()[3])

            # Get all other gradient norms for evaluation
            result["gncyc"] = [
                float(line.split()[3])
                for line in filter(lambda x: " :: gradient norm      " in x, lines)
            ]

            # store the final energy of the optimization in 'energy'
            result["energy"] = result["ecyc"][-1]
            meta["success"] = True

        if self.copy_mo:
            # store the path to the current .gbw file for this conformer
            job.meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # read out optimized geometry and update conformer geometry with this
        job.conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
        result["geom"] = job.conf.xyz

        # TODO - this might be a case where it would be reasonable to raise an
        # exception
        try:
            assert result["converged"] is not None
        except AssertionError:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_opt"

        job.meta["xtb_opt"].update(meta)

        return result

    def _nmr(self, job: ParallelJob, jobdir: str, filename: str = "nmr") -> dict[str, any]:
        """
        Calculate the NMR shieldings and/or couplings for a conformer. ORCA gives only the active cores in the output
        so there is not need for more thinking here.
        Formatting:
            'shielding' contains a list of tuples (atom_index, shielding), with atom_index being the index of the atom
            in the internal coordinates of the GeometryData.
            'couplings' contains a list of tuples (set(atom_index1, atom_index2), coupling), with the indices of the atoms
            in the internal coordinates of the GeometryData. A set is used to represent an atom pair.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory

        Returns:
            dict[str, any]: dictionary containing the results of the calculation
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

        # The following is a workaround to avoid a complicated nested if-else
        conds = (
            self.instructions["func_s"] == self.instructions["func_j"]
            and self.instructions["disp_s"] == self.instructions["disp_j"]
            and self.instructions["basis_s"] == self.instructions["basis_j"],
            self.instructions["shieldings"],
            self.instructions["couplings"],
        )

        conds_to_endings = {
            (True, False): ["_s"],
            (False, True): ["_j"],
            (True, True): ["_s", "_j"],
        }

        # Basically this selects the proper jobs to run, so run two calculations if the settings for shieldings and
        # couplings calculations differ even in one parameter, otherwise both
        # can be calculated in one go
        endings = []
        if all(cond for cond in conds):
            endings.append("")
        else:
            endings.extend(conds_to_endings[conds[1:]])

        # If settings are the same, endings = [""], otherwise it contains "_s"
        # and/or "_j", depending on conditions
        for ending in endings:
            job.meta.setdefault("nmr", {})

            # Set in/out path
            inputpath = os.path.join(jobdir, f"{filename}{ending}.inp")
            outputpath = os.path.join(jobdir, f"{filename}{ending}.out")

            # Prepare an input file for _sp
            indict = self.__prep(job, f"nmr{ending}")
            parser = OrcaParser()
            parser.write_input(inputpath, indict)

            # Run _sp using the input file generated above
            self._sp(job, jobdir, filename=f"{filename}{ending}", prep=False)

            if not job.meta["sp"]["success"]:
                meta["success"] = False
                meta["error"] = "sp failed"
                job.meta["nmr"].update(meta)
                return result

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

            # For shieldings watch out for the line "CHEMICAL SHIELDING SUMMARY
            # (ppm)"
            if ending in ["", "_s"]:
                start = lines.index(
                    next(x for x in lines if "CHEMICAL SHIELDING SUMMARY" in x)) + 6

                result["shieldings"] = []

                for line in lines[start:]:
                    # Try to extract values from the lines, if that fails
                    # (probably IndexError) stop
                    try:
                        result["shieldings"].append(
                            (int(line.split()[0]), float(line.split()[2])))
                    except IndexError:
                        break

            if ending in ["", "_j"]:
                # Read couplings from *_properties.txt for easier parsing
                with open(os.path.join(jobdir, f"{filename}{ending}_properties.txt"), "r") as f:
                    lines = f.readlines()

                start = lines.index(
                    next(x for x in lines if "EPRNMR_SSCoupling" in x)) + 13

                lines = lines[start:]

                end = lines.index(next(x for x in lines if "-----" in x))

                lines = lines[:end]

                result["couplings"] = []

                # This goes through all the pairs, even though ORCA gives also non-unique pairs, since it makes life
                # easier (basically every element of a symmetric square matrix)
                # Only every third line a new pair is defined
                # The iterator created here unpacks a tuple consisting of the multiples of three (indices of every
                # third line) and every third line in 'lines'
                for i, line in zip(range(0, len(lines), 3), lines[::3]):
                    pair = {int(line.split()[4]), int(line.split()[11])}
                    coupling = float(lines[i + 2].split()[4])
                    result["couplings"].append((pair, coupling))

                # Convert to set and back to get rid of duplicates
                # ('vectorizing the symmetric matrix')
                result["couplings"] = list(set(result["couplings"]))

        meta["success"] = True

        job.meta["nmr"].update(meta)

        return result

    @staticmethod
    def __apply_flags(indict: OrderedDict[str, any], *args) -> OrderedDict[str, any]:
        """
        apply flags to an orca input

        this is very much work in progress
        """
        flag_to_setting = {
            "scf_not_converged": {
                "scf": {
                    "maxiter": ["300"],
                    "AutoTRAH": ["false"],
                    "AutoTRAHIter": ["275"],
                    "CNVSOSCF": ["true"],
                    "SOSCFStart": ["0.0002"],
                    "SOSCFMaxIt": ["200"],
                },
                "main": ["veryslowconv"],
            },
        }

        for flag in args:
            if flag in flag_to_setting:
                # insert all settings dictated by the flags after the main
                # input line (TODO)
                for key, value in flag_to_setting[flag].items():
                    if key == "main":
                        indict[key].extend(value)
                    else:
                        indict = od_insert(indict, key, value, 1)

        return indict
