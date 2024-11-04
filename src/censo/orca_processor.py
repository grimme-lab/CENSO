"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""

import os
import shutil
from collections import OrderedDict
from functools import reduce

from .utilities import od_insert, Factory
from .logging import setup_logger
from .datastructure import GeometryData, ParallelJob
from .params import (
    Config,
    WARNLEN,
)
from .qm_processor import QmProc

logger = setup_logger(__name__)


class OrcaParser:
    """
    Parser for orca input files. Can read input files and transpile them to ordered dict. Also capable of writing an
    input file from a properly ordered dict (see __todict for format).
    """

    __exceptions = [
        "maxcore",
    ]

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

        "some setting": {"some option": [...], ...}, (contains a setting that is started with a % with the respective keywords and values)
        (the list under "some option" contains just all the following substrings, that, in the input file,
        are separated by a whitespace, nothing fancy)

        "coords": {"def": ["xyz/xyzfile", ...], "coord": [...]} (contains the geometry information block definition in 'def', 'coord' only if you use the option 'xyz')
        'coord' for e.g. H2: "coord": [["H", "0.0", "0.0", "0.0"], ["H", "0.0", "0.0", "0.7"]]

        "some setting after coords": {"some option": [...], ...}

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

                # check if there is already some option in the same line as the
                # setting declaration
                if len(split) > 1:
                    option = split[1]
                    converted[setting][option] = split[2:]

                if setting.lower() not in self.__exceptions:
                    # find end of definition block
                    end = i + self.__eob(lines[i:])

                    # in case the eob is found within the line itself remove the
                    # 'end' substring
                    if end == i:
                        converted[setting][option].remove("end")
                    # consume remaining definitions
                    else:
                        for line2 in lines[i + 1 : end]:
                            split = line2.split()
                            option = split[0]
                            converted[setting][option] = split[1:]
            # geometry input line found
            elif line.startswith("*") and "coords" not in converted.keys():
                converted["coords"] = {}

                if "xyzfile" in line:
                    converted["coords"]["def"] = ["xyzfile"]
                # the 'xyz' keyword should be one of the first two substrings
                elif "xyz" in line.split()[0] or "xyz" in line.split()[1]:
                    converted["coords"]["def"] = ["xyz"]
                else:
                    raise RuntimeError("Error parsing ORCA input file.")

                # add the remaining line to the dict
                # get rid of '*'
                line = line[1:]
                converted["coords"]["def"].extend(line.split()[1:])

                # consume the remaining coordsetry information
                if converted["coords"]["def"][0] == "xyz":
                    converted["coords"]["coord"] = []
                    # find end of definition block
                    # start search from next line since geometry definition
                    # starts with an '*'
                    end = i + self.__eob(lines[i + 1 :], endchar="*") + 1

                    for line2 in lines[i + 1 : end]:
                        converted["coords"]["coord"].append(line2.split())
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
            elif "{geom}" in line:
                # there should be no geometry definition block in a template
                # file
                if "coords" in converted.keys():
                    raise RuntimeError(
                        "Error parsing ORCA input file (double coords definition)."
                    )
                # also main needs to be defined before coords
                elif "main" not in converted.keys():
                    raise RuntimeError(
                        "Error parsing ORCA input file (missing main definition)."
                    )

                # if we reach this point, the coords key should be created
                converted["coords"] = {}

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
        for key in allkeys[1 : allkeys.index("coords")]:
            lines.append(f"%{key}\n")
            # FIXME - temporary workaround for definition blocks that have no
            # 'end', this code smells immensely
            try:
                for option in indict[key].keys():
                    lines.append(
                        f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n"
                    )
                    # NOTE: TypeError is raised if reduce fails due to indict[key][option] not being iterable
                lines.append("end\n")
            except TypeError:
                lines[-1] = f"%{key} {list(indict[key].keys())[0]}\n"

        # next, write the geometry input lines
        # geometry definition line (e.g. "* xyzfile 0 1 input.xyz" / "* xyz 0
        # 1")
        lines.append(f"* {reduce(lambda x, y: f'{x} {y}', indict['coords']['def'])}\n")

        # write coordinates if "xyz" keyword is used
        # if "xyzfile" is used, nothing more has to be done
        if indict["coords"].get("coord", False):
            for coord in indict["coords"]["coord"]:
                lines.append(f"{reduce(lambda x, y: f'{x} {y}', coord)}\n")
            lines.append("*\n")

        # lastly, write all the keywords and options that should be placed
        # after the geometry input (e.g. NMR settings)
        for key in allkeys[allkeys.index("coords") + 1 :]:
            lines.append(f"%{key}\n")
            for option in indict[key].keys():
                if "Nuclei" in option:
                    # Special treatment for the "Nuclei" option in %eprnmr
                    lines.append(
                        f"    {option[:len(option)-1]} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n"
                    )
                else:
                    lines.append(
                        f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n"
                    )
            lines.append("end\n")

        return lines


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

    # NOTE: this is currently not functionally used but could be used as a guideline
    __req_settings = {
        **{
            "sp": [
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
                "h_active",
                "c_active",
                "f_active",
                "si_active",
                "p_active",
            ],
        },
        **QmProc._req_settings_xtb,
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
    ) -> OrderedDict:
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

        # check ORCA version (orca5 = True means at least ORCA version 5)
        orca5 = not self._paths["orcaversion"].startswith("4")

        indict = OrderedDict()

        if job.prepinfo[jobtype]["template"]:
            # NOTE: when using templates we're not going to check for double definitions!
            # if the template is messed up, orca will fail and the user should deal with that
            # load template file
            try:
                indict = OrcaParser().read_input(
                    os.path.join(
                        Config.USER_ASSETS_PATH,
                        f"{job.prepinfo['partname']}.orca.template",
                    )
                )
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not find template file {job.prepinfo['partname']}.orca.template."
                )

        # prepare the main line of the orca input
        indict = self.__prep_main(job.prepinfo, indict, jobtype, orca5)

        # prepare all options that are supposed to be placed before the
        # geometry definition
        indict = self.__prep_pregeom(
            job.prepinfo, indict, orca5, jobtype, no_solv, job.omp
        )

        # prepare the geometry
        indict = self.__prep_geom(
            indict, job.conf, xyzfile, job.prepinfo["charge"], job.prepinfo["unpaired"]
        )

        indict = self.__prep_postgeom(job.prepinfo, indict, job.conf, jobtype, orca5)

        return indict

    def __prep_main(
        self, prepinfo: dict[str, any], indict: OrderedDict, jobtype: str, orca5: bool
    ) -> OrderedDict:
        if "main" not in indict:
            indict["main"] = []

        # grab func, basis
        func = prepinfo[jobtype]["func_name"]
        basis = prepinfo[jobtype]["basis"]
        functype = prepinfo[jobtype]["func_type"]
        disp = prepinfo[jobtype]["disp"]

        indict["main"].append(func)

        if "composite" not in functype:
            indict["main"].append(basis)

        # set  RI def2/J,   RIJCOSX def2/J

        # Set def2/J in case of def2 basis
        if "def2" in basis.lower():
            indict["main"].append("def2/J")
        # Otherwise use autoaux
        else:
            indict["main"].append("autoaux")

        # settings for double hybrids
        if "double" in functype:
            indict["main"].extend(["RIJCOSX"])

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

            def2cbasis = ("def2-svp", "def2-tzvp", "def2-tzvpp", "def2-qzvpp")
            if basis.lower() in def2cbasis:
                indict["main"].append(f"{basis}/C")

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

        # dummy type falls through every case, nothing is done in that case

        # use 'grid' setting from instructions to quickly configure the grid
        indict["main"].extend(self.__gridsettings[orca5][prepinfo[jobtype]["grid"]])

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
                    indict["main"].append(f"GCP(DFT/{gcp_keywords[basis.lower()]})")
                else:
                    logger.warning(
                        f"{f'worker{os.getpid()}:':{WARNLEN}}Selected basis not available for GCP. GCP not employed."
                    )

        # add job keyword for geometry optimizations
        # with ANCOPT
        if jobtype == "xtb_opt":
            indict["main"].extend(["ENGRAD", "tightSCF"])
        # for standard geometry optimization
        elif jobtype == "opt":
            indict["main"].extend(["OPT", "tightSCF"])

        # additional print settings
        if jobtype in ["xtb_opt", "opt"]:
            indict["main"].extend(["miniprint"])
        else:
            indict["main"].extend(["printgap"])

        return indict

    def __prep_pregeom(
        self,
        prepinfo: dict[str, any],
        indict: OrderedDict,
        orca5: bool,
        jobtype: str,
        no_solv: bool,
        nprocs: int,
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
        if (
            not prepinfo["general"]["gas-phase"]
            and not no_solv
            and ("sm" in prepinfo[jobtype].keys())
        ):
            sm = prepinfo[jobtype]["sm"]
            solv_key = prepinfo[jobtype]["solvent_key_prog"]

            if sm == "smd":
                indict = od_insert(
                    indict,
                    "cpcm",
                    {
                        "smd": ["true"],
                        "smdsolvent": [f'"{solv_key}"'],
                    },
                    list(indict.keys()).index("main") + 1,
                )
            elif sm == "cpcm":
                indict["main"].append(f"CPCM({solv_key})")

        if jobtype == "uvvis":
            indict = od_insert(
                indict,
                "tddft",
                {"nroots": [f"{prepinfo['uvvis']['nroots']}"]},
                list(indict.keys()).index("main") + 1,
            )

        # Additional print settings
        if jobtype not in ["xtb_opt", "opt"]:
            indict = od_insert(
                indict,
                "output",
                {"printlevel": ["normal"]},
                list(indict.keys()).index("main") + 1,
            )

        if jobtype == "opt":
            if prepinfo[jobtype]["macrocycles"]:
                # Set max number of optimization cycles for ORCA driven optimization
                indict = od_insert(
                    indict,
                    "geom",
                    {"maxiter": [prepinfo["opt"]["optcycles"]]},
                    list(indict.keys()).index("main") + 1,
                )
            else:
                indict = od_insert(
                    indict, "geom", {}, list(indict.keys()).index("main") + 1
                )

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
                indict["geom"]["convergence"] = [prepinfo["opt"]["optlevel"]]
            # Otherwise map to roughly corresponding orca optlevel
            else:
                indict["geom"]["convergence"] = [mapping[prepinfo["opt"]["optlevel"]]]

            # Insert constraints (if provided)
            # FIXME - without better parser this will not work
            """
            if prepinfo["opt"]["constraints"] is not None:
                parser = OrcaParser()
                constraints = parser.read_input(prepinfo["opt"]["constraints"])
                indict["geom"].update(constraints["geom"])
            """

        return indict

    def __prep_geom(
        self,
        indict: OrderedDict,
        conf: GeometryData,
        xyzfile: str,
        charge: int,
        unpaired: int,
    ) -> OrderedDict:
        # unpaired, charge, and coordinates
        # by default coordinates are written directly into input file
        if xyzfile is None:
            indict["coords"] = {
                "def": [
                    "xyz",
                    charge,
                    unpaired + 1,
                ],
                "coord": conf.toorca(),
            }
        else:
            indict["coords"] = {
                "def": [
                    "xyzfile",
                    charge,
                    unpaired + 1,
                    xyzfile,
                ],
            }

        return indict

    def __prep_postgeom(
        self,
        prepinfo: dict[str, any],
        indict: OrderedDict,
        conf: GeometryData,
        jobtype: str,
        orca5: bool,
    ) -> OrderedDict:
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
            indict.setdefault("eprnmr", {})
            if jobtype.endswith("_s") or jobtype == "nmr":
                todo2.append("shift")
                todo3["origin"] = ["giao"]
                todo3["giao_2el"] = ["giao_2el_same_as_scf"]
                todo3["giao_1el"] = ["giao_1el_analytic"]
            if jobtype.endswith("_j") or jobtype == "nmr":
                if prepinfo[jobtype]["fc_only"]:
                    todo2.append("ssfc")
                else:
                    todo2.append("ssall")
                todo3["SpinSpinRThresh"] = [f"{prepinfo[jobtype]['ss_cutoff']:.4f}"]

            # Workaround for weird orca quirk:
            # Creates a "Nuclei" entry for every active nucleus, which is then parsed
            # into multiple lines with the "Nuclei" keyword by the parser
            # The numbers (Nuclei{i}) are stripped by the parser and necessary because
            # in a dict all keys are unique so we cannot have multiple "Nuclei"
            nuclei = {}
            for i, element in enumerate(todo):
                nuclei[f"Nuclei{i}"] = (
                    ["="] + ["all", element] + ["{", ",".join(x for x in todo2), "}"]
                )

            compiled = {**nuclei, **todo3}

            # Insert the settings for NMR calculation
            indict["eprnmr"] = compiled

        return indict

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
    ) -> tuple[dict[str, float | None], dict[str, any]]:
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
            indict = self.__prep(job, "sp", no_solv=no_solv)

            # check for flags raised for this jobtype
            # NOTE: all other jobtypes call this function so flags are always checked here
            # except xtb_opt
            indict = self.__apply_flags(job, indict)

            # write input into file "{filename}.inp" in a subdir created for the
            # conformer
            parser = OrcaParser()
            parser.write_input(inputpath, indict)

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
    ) -> tuple[dict[str, any], dict[str, any]]:
        """
        Calculates the solvation free enthalpy of a conformer using ORCA.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory

        Returns:
            result (dict[str, any]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job

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
    ) -> tuple[dict[str, any], dict[str, any]]:
        """
        Geometry optimization using ORCA optimizer.
        Note that solvation in handled here always implicitly.

        Args:
            job: ParallelJob object containing the job information, metadata is stored in job.meta
            jobdir: path to the job directory
            filename: name of the input file

        Returns:
            result (dict[str, any]): dictionary containing the results of the calculation
            meta (dict[str, any]): metadata about the job
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
        # TODO - add constraints
        parser = OrcaParser()
        indict = self.__prep(job, "opt")

        # apply flags
        indict = self.__apply_flags(job, indict)

        # write orca input in a subdir created for the conformer
        parser.write_input(inputpath, indict)

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

        # Get final energy
        result["energy"] = next(
            (
                float(line.split()[4])
                for line in lines
                if "FINAL SINGLE POINT ENERGY" in line
            ),
            None,
        )

        meta["error"] = self.__check_output(lines)
        meta["success"] = meta["error"] is None and result["energy"] is not None

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
        parser = OrcaParser()
        indict = self.__prep(job, "xtb_opt", xyzfile=f"{filename}.xyz")

        # apply flags
        indict = self.__apply_flags(job, indict)

        # write orca input into file "xtb_opt.inp" in a subdir created for the
        # conformer
        parser.write_input(inputpath, indict)

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
    ) -> tuple[dict[str, any], dict[str, any]]:
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

        # If settings are the same, endings = [""], otherwise it contains "_s"
        # and/or "_j", depending on conditions
        for ending in [x[4:] for x in job.prepinfo.keys() if "nmr" in x]:
            # Set in/out path
            inputpath = os.path.join(jobdir, f"{filename}{ending}.inp")
            outputpath = os.path.join(jobdir, f"{filename}{ending}.out")

            # Prepare an input file for _sp
            indict = self.__prep(job, f"nmr{ending}")
            parser = OrcaParser()
            parser.write_input(inputpath, indict)

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
                # Read couplings from *_properties.txt for easier parsing
                with open(
                    os.path.join(jobdir, f"{filename}{ending}_property.txt"), "r"
                ) as f:
                    lines = f.readlines()

                start = (
                    lines.index(next(x for x in lines if "EPRNMR_SSCoupling" in x)) + 13
                )

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
                    # pair needs to be a frozenset because normal sets are not hashable and can therefore not be part
                    # of a normal set
                    pair = frozenset((int(line.split()[4]), int(line.split()[11])))
                    coupling = float(lines[i + 2].split()[4])
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
    ) -> tuple[dict[str, any], dict[str, any]]:
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
        indict = self.__prep(job, "uvvis")
        parser = OrcaParser()
        parser.write_input(inputpath, indict)

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
                {"wavelength": float(spl[2]), "osc_str": float(spl[3])}
            )

        meta["success"] = True

        return result, meta

    @staticmethod
    def __apply_flags(
        job: ParallelJob, indict: OrderedDict[str, any]
    ) -> OrderedDict[str, any]:
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
        flags = [flag for flag in job.flags.values() if flag in flag_to_setting]

        for flag in flags:
            # insert all settings dictated by the flags after the main
            for key, value in flag_to_setting[flag].items():
                if key == "main":
                    indict[key].extend(value)
                else:
                    indict = od_insert(indict, key, value, 1)

        return indict


Factory.register_builder("orca", OrcaProc)
