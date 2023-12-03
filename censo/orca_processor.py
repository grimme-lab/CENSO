"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""
import os
import shutil
from collections import OrderedDict
from functools import reduce
from typing import Any, List

from censo.datastructure import GeometryData, ParallelJob
from censo.params import (
    CODING,
    USER_ASSETS_PATH, WARNLEN,
)
from censo.qm_processor import QmProc
from censo.utilities import od_insert, setup_logger

logger = setup_logger(__name__)


class OrcaParser:
    """
    Parser for orca input files
    can read input files and transpile them to ordered dict
    also capable of writing an input file from a properly ordered dict (see __todict for format)
    """

    def __init__(self):
        pass

    def read_input(self, path: str) -> OrderedDict:
        """
        read orca input file at 'path' and parse into ordered dict
        """
        with open(path, "r") as infile:
            lines = infile.readlines()

        lines = self.__remove_comments(lines)
        converted = self.__todict(lines)

        return converted

    def write_input(self, path: str, indict: OrderedDict) -> None:
        """
        write an orca input file at 'path' from an ordered dict
        """
        with open(path, "w") as outfile:
            outfile.writelines(self.__tolines(indict))

    def __todict(self, lines: List[str]) -> OrderedDict:
        """
        convert lines from orca input into ordered dict

        can also be used with an orca input template, the {main} and {postgeom} strings indicate the ordering
        of the settigs read from the template ({main} should come first, {postgeom} should come last)
        all settings below {postgeom} are put after the geometry definition block

        format of the result:
        "main": [...], (contains all main input line options)
        
        "some setting": {"some option": [...], ...}, (contains a settings that is started with a % with the respective keywords and values) 
        (the list under "some option" contains just all the following substrings, nothing fancy)
        
        "geom": {"def": ["xyz/xyzfile", ...], "coord": [...]} (contains the geometry information block definition in 'def', 'coord' only if you use the option 'xyz')
        'coord' for e.g. H2: "coord": [["H", "0.0", "0.0", "0.0"], ["H", "0.0", "0.0", "0.7"]]

        "some setting after geom": {"some option": [...], ...}

        comments get removed from the lines and therefore lost

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

                # check it there is already some option in the same line as the setting declaration
                if len(split) > 1:
                    option = split[1]
                    converted[setting][option] = split[2:]

                # find end of definition block
                end = i + self.__eob(lines[i:])

                # in case the eob is found within the line itself remove the 'end' substring
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
                    raise RuntimeError("Error parsing ORCA input file (double geom definition).")

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
                    # start search from next line since geometry definition starts with an '*'
                    end = i + self.__eob(lines[i + 1:], endchar="*") + 1

                    for line2 in lines[i + 1:end]:
                        converted["geom"]["coord"].append(line2.split())
            # check for template lines
            # NOTE: only these two need to be checked since they're the only ones that are sensitive to ordering
            elif "{main}" in line:
                # create the main key for ordering if it does not already exist
                if "main" not in converted.keys():
                    converted["main"] = []

                # main needs to be defined before everything
                if "main" != list(converted.keys())[0]:
                    raise RuntimeError("Error parsing ORCA input file (main not defined first).")
            elif "{postgeom}" in line:
                # there should be no geometry definition block in a template file
                if "geom" in converted.keys():
                    raise RuntimeError("Error parsing ORCA input file (double geom definition).")
                # also main needs to be defined before geom
                elif "main" not in converted.keys():
                    raise RuntimeError("Error parsing ORCA input file (missing main definition).")

                # if we reach this point, the geom key should be created
                converted["geom"] = {}

        return converted

    @staticmethod
    def __eob(lines: List[str], endchar: str = "end") -> int:
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
    def __remove_comments(inlist: List[str]) -> List[str]:
        """
        remove all comments from an orca input
        """
        for i in range(len(inlist)):
            if "#" in inlist[i]:
                index = inlist[i].index("#")
                inlist[i] = inlist[i][:index]

        return inlist

    @staticmethod
    def __tolines(indict: OrderedDict) -> List[str]:
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

        # next, write all keywords and options that come between the main input line and the geom input
        allkeys = list(indict.keys())

        # skip first key ('main')
        for key in allkeys[1:allkeys.index("geom")]:
            lines.append(f"%{key}\n")
            # FIXME - temporary workaround for definition blocks that have no 'end', this code smells immensely
            try:
                for option in indict[key].keys():
                    lines.append(f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n")
                lines.append("end\n")
            except TypeError:
                lines[-1] = f"%{key} {list(indict[key].keys())[0]}\n"

        # next, write the geometry input lines
        # geometry definition line (e.g. "* xyzfile 0 1 input.xyz" / "* xyz 0 1")
        lines.append(f"* {reduce(lambda x, y: f'{x} {y}', indict['geom']['def'])}\n")

        # write coordinates if "xyz" keyword is used
        # if "xyzfile" is used, nothing more has to be done
        if indict["geom"].get("coord", False):
            for coord in indict["geom"]["coord"]:
                lines.append(f"{reduce(lambda x, y: f'{x} {y}', coord)}\n")
            lines.append("*\n")

        # lastly, write all the keywords and options that should be placed after the geometry input (e.g. NMR settings)
        for key in allkeys[allkeys.index("geom") + 1:]:
            lines.append(f"%{key}\n")
            for option in indict[key].keys():
                lines.append(f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n")
            lines.append("end\n")

        return lines


class OrcaProc(QmProc):
    """
    Perform calculations with ORCA
    - create orca.inp input
    - single-point calculation
    - smd_gsolv calculation
    - geometry optimization with xTB as driver
    - shielding constant calculations
    - coupling constant calculations
    - writing of generic output for shielding and coupling constants
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # expand jobtypes with special orca jobtypes
        """self._jobtypes = {
            **self._jobtypes, **{
                "nmrS": self._nmrS,
                "nmrJ": self._nmrJ,
                "uvvis": self._uvvis,
            }
        }"""

        # contains grid settings for ORCA 5.0+ (True) and older versions (False)
        # can be chosen by simple keyword (low/low+/high/high+) 
        self.__gridsettings = {
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

    def __prep(self, job: ParallelJob, jobtype: str, no_solv: bool = False, xyzfile: str = None) -> OrderedDict:
        """
        prepare an OrderedDict to be fed into the parser in order to write an input file
        for jobtype 'jobtype' (e.g. sp)

        can load up a template file from user assets folder

        xyzfile: name of the xyz-file to be used for the calculation

        NOTE: the xyzfile has to already exist, this function just bluntly writes the name into the input

        TODO - NMR/OR/UVVis etc preparation steps
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
                    os.path.join(USER_ASSETS_PATH, f"{self.instructions['part_name']}.orca.template"))
            except FileNotFoundError:
                raise FileNotFoundError(f"Could not find template file {self.instructions['part_name']}.orca.template.")

        # prepare the main line of the orca input
        indict = self.__prep_main(indict, jobtype, orca5)

        # TODO - add special lines if this is a retry

        # prepare all options that are supposed to be placed before the geometry definition
        indict = self.__prep_pregeom(indict, orca5, no_solv, job.omp)

        # prepare the geometry
        indict = self.__prep_geom(indict, job.conf, xyzfile)

        # indict = self.__prep_postgeom(indict, jobtype, orca5)

        return indict

    def __prep_main(self, indict: OrderedDict, jobtype: str, orca5: bool) -> OrderedDict:

        if "main" not in indict:
            indict["main"] = []

        # grab func, basis
        func = self.instructions["func_name"]
        basis = self.instructions["basis"]
        indict["main"].append(func)

        if "composite" not in self.instructions["func_type"]:
            indict["main"].append(basis)
            # set  RI def2/J,   RIJCOSX def2/J
            # this is only set if no composite DFA is used
            # settings for double hybrids
            if self.instructions["func_type"] == "double":
                indict["main"].extend(["def2/J", "RIJCOSX"])

                if "nmr" in jobtype:
                    indict["main"].append("NOFROZENCORE")
                    indict = od_insert(
                        indict,
                        "mp2",
                        {"RI": ["true"], "density": ["relaxed"]},
                        list(indict.keys()).index("main") + 1
                    )
                else:
                    indict["main"].append("frozencore")
                    indict = od_insert(
                        indict,
                        "mp2",
                        {"RI": ["true"]},
                        list(indict.keys()).index("main") + 1
                    )

                def2cbasis = ("def2-svp", "def2-tzvp", "def2-tzvpp", "def2-qzvpp")
                if self.instructions["basis"].lower() in def2cbasis:
                    indict["main"].append(f"{self.instructions['basis']}/C")
                    if not orca5:
                        indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])
                else:
                    indict["main"].append("def2-TZVPP/C")
                    if not orca5:
                        indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])

            # settings for hybrids
            elif "hybrid" in self.instructions["func_type"]:
                indict["main"].append("RIJCOSX")
                if not orca5:
                    indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])

            # settings for (m)ggas
            elif "gga" in self.instructions["func_type"]:
                indict["main"].append("RI")

        # use 'grid' setting from instructions to quickly choose the grid settings
        indict["main"].extend(self.__gridsettings[orca5][self.instructions["grid"]])

        # add dispersion
        # dispersion correction information
        # FIXME - temporary solution (not very nice)
        mapping = {
            "d3bj": "d3bj",
            "d3(0)": "D3ZERO",
            "d4": "d4",
            "nl": "NL",
        }

        disp = self.instructions["disp"]
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
        if self.instructions["gcp"] and self.instructions["basis"].lower() in gcp_keywords.keys():
            indict["main"].append(f"GCP(DFT/{gcp_keywords[self.instructions['basis'].lower()]})")
        elif self.instructions["gcp"]:
            # TODO - error handling
            global logger
            logger.warning(f"{f'worker{os.getpid()}:':{WARNLEN}}Selected basis not available for GCP.")

        # add job keyword for geometry optimizations
        # with ANCOPT
        if jobtype == "xtb_opt":
            indict["main"].append("ENGRAD")
        # for standard geometry optimization
        elif jobtype == "opt":
            indict["main"].append("OPT")

        return indict

    def __prep_pregeom(self, indict: OrderedDict, orca5: bool, no_solv: bool, nprocs: int) -> OrderedDict:

        """used in nmr calculation
        if self.job["moread"] is not None:
            #! MORead
            #%moinp "jobname2.gbw"
            orcainput["moread"] = self.job["moread"]"""

        if orca5 and nprocs > 1:
            indict = od_insert(
                indict,
                "pal",
                {"nprocs": [nprocs]},
                list(indict.keys()).index("main") + 1
            )

        # TODO - maybe limit TRAH macro steps (or when TRAH activates) to avoid single jobs clogging everything up

        # set keywords for the selected solvent model
        if not self.instructions["gas-phase"] and not no_solv and "sm" in self.instructions.keys():
            if self.instructions["sm"] == "smd":
                indict = od_insert(
                    indict,
                    "cpcm",
                    {"smd": ["true"], "smdsolvent": [f"\"{self.instructions['solvent_key_prog']}\""]},
                    list(indict.keys()).index("main") + 1
                )
            elif self.instructions["sm"] == "cpcm":
                indict["main"].append(f"CPCM({self.instructions['solvent_key_prog']})")

        return indict

    def __prep_geom(self, indict: OrderedDict, conf: GeometryData, xyzfile: str) -> OrderedDict:
        # unpaired, charge, and coordinates
        # by default coordinates are written directly into input file
        if xyzfile is None:
            indict["geom"] = {
                "def": ["xyz", self.instructions["charge"], self.instructions["unpaired"] + 1],
                "coord": conf.toorca(),
            }
        else:
            indict["geom"] = {
                "def": ["xyzfile", self.instructions["charge"], self.instructions["unpaired"] + 1, xyzfile],
            }

        return indict

    def __prep_postgeom(self, indict: OrderedDict, jobtype: str, orca5: bool) -> OrderedDict:
        pass

    @QmProc._create_jobdir
    def _sp(self, job: ParallelJob, silent=False, filename="sp", no_solv: bool = False, jobtype: str = "sp") -> dict[
        str, float | None]:
        """
        ORCA single-point calculation

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

        # set in/out path
        jobdir = os.path.join(self.workdir, job.conf.name, jobtype)
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        indict = self.__prep(job, "sp", no_solv=no_solv)

        # check for flags raised for this jobtype
        if jobtype in job.flags:
            if job.flags[jobtype] == "scf_not_converged":
                indict = self.__apply_flags(indict, "scf_not_converged")

        # write input into file "{filename}.inp" in a subdir created for the conformer
        parser = OrcaParser()
        parser.write_input(inputpath, indict)

        # check, if there is an existing .gbw file and copy it if option 'copy_mo' is true
        global logger
        if self.instructions["copy_mo"]:
            if job.mo_guess is not None and os.path.isfile(job.mo_guess):
                logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {job.mo_guess}.")
                shutil.copy(job.mo_guess, os.path.join(jobdir, f"{filename}.gbw"))

        if not silent:
            logger.info(f"{f'worker{os.getpid()}:':{WARNLEN}}Running ORCA single-point in {inputpath}")
        else:
            logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Running ORCA single-point in {inputpath}")

        # call orca
        call = [self._paths["orcapath"], f"{filename}.inp"]
        self._make_call(call, outputpath, jobdir)

        # do not check returncode, since orca doesn't give meaningful returncodes
        # check the outputfile instead, basically if the file doesn't say "SCF CONVERGED" or "SCF NOT CONVERGED" somewhere, something went wrong
        # easily remedied:
        # "SCF NOT CONVERGED"
        # "Error encountered when trying to calculate the atomic fitting density!"

        # read output
        with open(outputpath, "r", encoding=CODING, newline=None) as out:
            lines = out.readlines()

            # Get final energy
            result["energy"] = next((float(line.split()[4]) for line in lines if "FINAL SINGLE POINT ENERGY" in line),
                                    None)

            # Check if scf is converged
            meta["success"] = next((True for line in lines if "SCF CONVERGED" in line), False)

        # TODO - this is not really correct, might not mean that the scf didn't converge
        if not meta["success"] is True:
            meta["error"] = "SCF not converged"

        if self.instructions["copy_mo"]:
            # store the path to the current .gbw file for this conformer if possible
            if os.path.isfile(os.path.join(jobdir, f"{filename}.gbw")):
                job.meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # TODO - clean up?
        job.meta[jobtype].update(meta)

        return result

    @QmProc._create_jobdir
    def _gsolv(self, job: ParallelJob) -> dict[str, Any | None]:
        """
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
        }

        jobdir = os.path.join(self.workdir, job.conf.name, "gsolv")

        global logger
        logger.info(f"{f'worker{os.getpid()}:':{WARNLEN}}Running ORCA Gsolv calculation in {jobdir}.")

        # calculate gas phase
        # TODO - this is redundant since a single point was probably already calculated before
        # TODO - does this need it's own folder? this is a bit messy
        spres = self._sp(job, silent=True, filename="sp_gas", no_solv=True, jobtype="gsolv")

        if job.meta["gsolv"]["success"]:
            result["energy_gas"] = spres["energy"]
        else:
            meta["success"] = False
            meta["error"] = "SCF not converged"
            job.meta["gsolv"].update(meta)
            return result

        # calculate in solution
        spres = self._sp(job, silent=True, filename="sp_solv", jobtype="gsolv")

        if job.meta["gsolv"]["success"]:
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
    @QmProc._create_jobdir
    def _xtb_opt(self, job: ParallelJob, filename: str = "xtb_opt") -> dict[str, Any]:
        """
        ORCA geometry optimization using ANCOPT

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
        # 'energy' contains the final energy of the optimization (converged or unconverged)
        # 'geom' stores the optimized geometry in GeometryData.xyz format
        result = {
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "geom": None,
        }

        meta = {
            "success": None,
            "error": None,
        }

        jobdir = os.path.join(self.workdir, job.conf.name, "xtb_opt")
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

        # write orca input into file "xtb_opt.inp" in a subdir created for the conformer
        parser.write_input(inputpath, indict)

        # append some additional lines to the coord file for ancopt
        with open(
                os.path.join(jobdir, f"{filename}.coord"), "a", newline=None
        ) as newcoord:
            newcoord.writelines([
                "$external\n",
                f"   orca input file= {filename}.inp\n",
                f"   orca bin= {self._paths['orcapath']}\n",
                "$end\n"
            ])

        # prepare configuration file for ancopt (xcontrol file)
        with open(
                os.path.join(jobdir, xcontrolname), "w", newline=None
        ) as out:
            out.write("$opt \n")
            if self.instructions["opt_spearman"]:
                out.write(f"maxcycle={self.instructions['optcycles']} \n")
                out.write(f"microcycle={self.instructions['optcycles']} \n")

            out.writelines([
                "average conv=true \n",
                f"hlow={self.instructions['hlow']} \n",
                "s6=30.00 \n",
                # remove unnecessary sp/gradient call in xTB
                "engine=lbfgs\n",
                "$external\n",
                f"   orca input file= {filename}.inp\n",
                f"   orca bin= {self._paths['orcapath']} \n",
                "$end \n",
            ])

        # check, if there is an existing .gbw file and copy it if option 'copy_mo' is true
        global logger
        if self.instructions["copy_mo"]:
            if job.mo_guess is not None and os.path.isfile(job.mo_guess):
                logger.debug(f"{f'worker{os.getpid()}:':{WARNLEN}}Copying .gbw file from {job.mo_guess}.")
                shutil.copy(job.mo_guess, os.path.join(jobdir, f"{filename}.gbw"))

        # prepare xtb call
        call = [
            self._paths["xtbpath"],
            f"{filename}.coord",  # name of the coord file generated above
            "--opt",
            self.instructions["optlevel"],
            "--orca",
            "-I",
            xcontrolname
        ]

        # set path to the ancopt output file
        outputpath = os.path.join(jobdir, f"{filename}.out")

        logger.info(f"{f'worker{os.getpid()}:':{WARNLEN}}Running optimization in {jobdir}.")

        # call xtb
        returncode = self._make_call(call, outputpath, jobdir)

        # check if optimization finished without errors
        if returncode != 0:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_opt"
            job.meta["xtb_opt"].update(meta)
            return result

        # read output
        with (open(outputpath, "r", encoding=CODING, newline=None) as file):
            lines = file.readlines()

        result["ecyc"] = []
        result["cycles"] = 0

        # Substrings indicating error in xtb
        error_ind = [
            "external code error",
            "|grad| > 500, something is totally wrong!",
            "abnormal termination of xtb"
        ]

        # Check if xtb terminated normally (if there are any error indicators in the output)
        meta["success"] = False if next((x for x in lines if any(y in x for y in error_ind)),
                                        None) is not None else True
        if not meta["success"]:
            meta["error"] = "what went wrong in xtb_opt"
            job.meta["xtb_opt"].update(meta)
            return result

        # check convergence
        if next((True for x in lines if "GEOMETRY OPTIMIZATION CONVERGED" in x), None) is True:
            result["converged"] = True
        elif next((True for x in lines if "FAILED TO CONVERGE GEOMETRY" in x), None) is True:
            result["converged"] = False

        # Get the number of cycles
        if result["converged"] is not None:
            # tmp is one of the values from the dict defined below
            tmp = {
                True: ("GEOMETRY OPTIMIZATION CONVERGED", 7),
                False: ("FAILED TO CONVERGE GEOMETRY", 5),
            }
            tmp = tmp[result["converged"]]

            result["cycles"] = next(x for x in lines if tmp[0] in x).split()[tmp[1]]

            # Get energies for each cycle
            result["ecyc"].extend(float(line.split("->")[-1]) for line in filter(lambda x: "av. E: " in x, lines))

            # Get the gradient norm
            result["grad_norm"] = float(next((x for x in lines if " :: gradient norm      " in x), None).split()[3])

            # store the final energy of the optimization in 'energy'
            result["energy"] = result["ecyc"][-1]
            meta["success"] = True

        if self.instructions["copy_mo"]:
            # store the path to the current .gbw file for this conformer
            job.meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # read out optimized geometry und update conformer geometry with this
        job.conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
        result["geom"] = job.conf.xyz

        # TODO - this might be a case where it would be reasonable to raise an exception
        try:
            assert result["converged"] is not None
        except AssertionError:
            meta["success"] = False
            meta["error"] = "what went wrong in xtb_opt"

        job.meta["xtb_opt"].update(meta)

        return result

    @staticmethod
    def __apply_flags(indict: OrderedDict[str, Any], *args) -> OrderedDict[str, Any]:
        """
        apply flags to an orca input

        this is very much work in progress
        """
        flag_to_setting = {
            "scf_not_converged": {
                "scf": {"maxiter": ["300"]},
                "main": ["slowconv"],
            },
        }

        for flag in args:
            if flag in flag_to_setting:
                # insert all settings dictated by the flags after the main input line (TODO)
                for key, value in flag_to_setting[flag].items():
                    if key == "main":
                        indict[key].extend(value)
                    else:
                        indict = od_insert(indict, key, value, 1)

        return indict
