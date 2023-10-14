"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""
from collections import OrderedDict
import os
import sys
import time
import subprocess
from typing import Any, Dict, List
from functools import reduce
import shutil

from censo.cfg import (
    CODING,
    ENVIRON,
    WARNLEN,
    PLANCK,
    C,
    SOLV_MODS,
    GSOLV_MODS,
)
from censo.utilities import last_folders, t2x, x2t, print
from censo.qm_processor import QmProc
from censo.errorswarnings import ParsingError
from censo.datastructure import GeometryData


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

        # remove all comments from read lines
        self.__remove_comments(lines)

        for i, line in enumerate(lines):
            # strip leading whitespace
            line = line.lstrip()

            # main input line found
            if line.startswith("!"):
                # get rid of '!'
                line =  line[1:]

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
                    for line2 in lines[i+1:end]:
                        split = line2.split()
                        option = split[0]
                        converted[setting][option] = split[1:]
            # geometry input line found
            elif line.startswith("*") and "geom" not in converted.keys():
                converted["geom"] = {}
                if "xyzfile" in line:
                    converted["geom"]["def"] = ["xyzfile"]
                # the 'xyz' keyword should be one of the first two substrings
                elif "xyz" in line.split()[0] or "xyz" in line.split()[1]:
                    converted["geom"]["def"] = ["xyz"]
                else:
                    raise ParsingError()
                
                # add the remaining line to the dict
                # get rid of '*'
                line = line[1:]
                converted["geom"]["def"].extend(line.split()[1:])

                # consume the remaining geometry information
                if converted["geom"]["def"][0] == "xyz":
                    converted["geom"]["coord"] = []
                    # find end of definition block
                    # start search from next line since geometry definition starts with an '*'
                    end = i + self.__eob(lines[i+1:], endchar="*") + 1

                    for line2 in lines[i+1:end]:
                        converted["geom"]["coord"].append(line2.split())

        return converted


    def __eob(self, lines: List[str], endchar: str = "end") -> int:
        """
        find end of definition block
        """
        end = None
        for index, line in enumerate(lines):
            if endchar in line:
                end = index
                break

        if end is None:
            raise ParsingError()
        else:
            return end


    def __remove_comments(self, inlist: List) -> None:
        """
        remove all comments from an orca input
        """
        for i in range(len(inlist)):
            if "#" in inlist[i]:
                index = inlist[i].index("#")
                inlist[i] = inlist[i][:index]
    
    def __tolines(self, indict: OrderedDict) -> List[str]:
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
            for option in indict[key].keys():
                lines.append(f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n")
            lines.append("end\n")
        
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
        for key in allkeys[allkeys.index("geom")+1:]:
            lines.append(f"%{key}\n")
            for option in indict[key].keys():
                lines.append(f"    {option} {reduce(lambda x, y: f'{x} {y}', indict[key][option])}\n")
            lines.append("end\n")

        return lines


# TODO - keep output if any job fails
class OrcaProc(QmProc):
    """
    Perform calculations with ORCA
    - create orca.inp input
    - single-point calculation
    - smd_gsolv calculation
    - optimization with xTB as driver
    - shielding constant calculations
    - coupling constant calculations
    - writing of generic output for shielding and coupling constants
    """

    # TODO - Lookup dicts for solvent models and functionals/dispersion corretions
    _solvent_models = {
        "sm": SOLV_MODS["orca"],
        "smgsolv": GSOLV_MODS["orca"],
    }

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
        

    # TODO - make this better
    def __prep(self, conf: GeometryData, jobtype: str, no_solv: bool = False, xyzfile: str = None) -> OrderedDict:
        """
        prepare an OrderedDict to be fed into the parser in order to write an input file
        for jobtype 'jobtype' (e.g. sp)

        xyzfile: name of the xyz-file to be used for the calculation

        NOTE: the xyzfile has to already exist, this function just bluntly writes the name into the input!

        TODO - NMR/OR/UVVis etc preparation steps
        """

        # check ORCA version
        orca5 = True if self.instructions["orcaversion"].startswith("5") else False

        indict = OrderedDict()
        indict["main"] = []

        # grab func, basis
        func = self.instructions["func_name"]
        basis = self.instructions["basis"] 
        indict["main"].append(func)

        # TODO - b3lyp-3c removed
        if "composite" not in self.instructions["func_type"]:
            indict["main"].append(basis)
        ####################### SET RI ###########################
            # set  RI def2/J,   RIJCOSX def2/J
            # this is only set if no composite DFA is used
            # settings for double hybrids
            if self.instructions["func_type"] == "double":
                indict["main"].extend(["def2/J", "RIJCOSX"])

                if "nmr" in jobtype:
                    indict["main"].append("NOFROZENCORE")
                    indict["mp2"] = {
                        "RI": ["true"],
                        "density": ["relaxed"],
                    }
                else:
                    indict["main"].append("frozencore")
                    indict["mp2"] = {
                        "RI": ["true"],
                    }

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

        ########################## SET GRID ############################ 

        """
        used in nmr calculation
        if self.job["moread"] is not None:
            #! MORead
            #%moinp "jobname2.gbw"
            orcainput["moread"] = self.job["moread"]"""

        # use 'grid' setting from instructions to quickly choose the grid settings 
        indict["main"].extend(self.__gridsettings[orca5][self.instructions["grid"]])
        
        ########################## DISPERSION ####################### TODO TODO TODO TODO TODO
        # add dispersion
        # dispersion correction information
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

        """
        else:
            print(
                f" {self.dfa_settings.functionals.get(self.instructions['func']).get('disp')} unknown dispersion option!"
            )
        """

        ###################### PARALLEL ############################
        if orca5 and self.instructions["omp"] > 1:
            indict["pal"] = {"nprocs": [self.instructions["omp"]]}

        ###################### SOLVENT/GEOM ########################

        # set keywords for the selected solvent model
        if not self.instructions["gas-phase"] and not no_solv and "sm" in self.instructions.keys():
            if self.instructions["sm"] == "smd":
                indict["cpcm"] = {
                    "smd": ["true"],
                    "smdsolvent": [f"\"{self.instructions['solvent_key_prog']}\""],
                }
            elif self.instuctions["sm"] == "cpcm":
                indict["main"].append(f"CPCM({self.instructions['solvent_key_prog']})")

        # unpaired, charge, and coordinates
        # by default coordinates are written directly into input file
        if xyzfile is None:
            indict["geom"] = {
                "def": ["xyz", self.instructions["charge"], self.instructions["unpaired"]+1],
                "coord": conf.toorca(),
            }
        else:
            indict["geom"] = {
                "def": ["xyzfile", self.instructions["charge"], self.instructions["unpaired"]+1, xyzfile],
            }

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
            print("Selected basis not available for GCP.")

        # add job keyword for geometry optimizations
        # with ANCOPT
        if jobtype == "xtb_opt":
            indict["main"].append("ENGRAD")
        # for standard optimization
        elif jobtype == "opt":
            indict["main"].append("OPT")

        return indict


    @QmProc._create_jobdir
    def _sp(self, conf: GeometryData, silent=False, filename="sp", no_solv: bool = False) -> Dict[str, Any]:
        """
        ORCA single-point calculation

        result = {
            "energy": None,
            "success": None,
            "mo_path": None,
        }
        """
        # set results
        result = {
            "energy": None,
            "success": None,
            "mo_path": None,
        }

        # set in/out path
        jobdir = os.path.join(self.workdir, conf.name, self._sp.__name__[1:])
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # prepare input dict
        indict = self.__prep(conf, "sp", no_solv=no_solv)
        
        # write input into file "{filename}.inp" in a subdir created for the conformer
        parser = OrcaParser()
        parser.write_input(inputpath, indict)

        # check, if there is an existing .gbw file and copy it if option 'copy_mo' is true
        if self.instructions["copy_mo"]:
            if getattr(conf, "mo_path") is not None and os.path.isfile(conf.mo_path):
                print(f"Copying .gbw file from {conf.mo_path}.")
                shutil.copy(conf.mo_path, os.path.join(jobdir, f"{filename}.gbw"))

        if not silent:
            print(f"Running ORCA single-point in {inputpath}")
        
        # call orca
        call = [self.instructions["orcapath"], f"{filename}.inp"]
        returncode = self._make_call(call, outputpath, jobdir)

        # TODO - check returncode?

        # read output
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as out:
                lines = out.readlines()
                for i, line in enumerate(lines):
                    if "FINAL SINGLE POINT ENERGY" in line:
                        result["energy"] = float(line.split()[4])

                    # check if scf is converged:
                    if "ORCA TERMINATED NORMALLY" in line:
                        result["success"] = True
                    
            if not result["success"]:
                result["success"] = False
                print(
                    f"{'ERROR:':{WARNLEN}}ORCA single-point not converged for {conf.name}."
                )
        else:
            # TODO - error handling
            result["success"] = False
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")

        # store the path to the current .gbw file for this conformer
        result["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # TODO - clean up

        return result


    def _gsolv(self, conf: GeometryData):
        """
        result = {
            "success": None,
            "gsolv": None,
            "energy_gas": None,
            "energy_solv": None,
        }
        """
        # what is returned in the end
        result = {
            "success": None,
            "gsolv": None,
            "energy_gas": None,
            "energy_solv": None,
        }

        # set in/out path

        print(
            f"Running SMD_gsolv calculation in "
            f"{last_folders(self.workdir, 2)}."
        )

        # calculate gas phase
        # TODO - this is redundant since a single point was probably already calculated before
        # TODO - does this need it's own folder?
        res = self._sp(conf, silent=True, filename="sp_gas", no_solv=True)

        if res["success"]:
            result["energy_gas"] = res["energy"]
        else:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}in gas phase single-point "
                f"of {last_folders(self.workdir, 2):18}"
            )
            return result

        # calculate in solution
        res = self._sp(conf, silent=True, filename="sp_solv")

        if res["success"]:
            result["energy_solv"] = res["energy"]
        else:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}in gsolv single-point "
                f"of {last_folders(self.workdir, 2):18}"
            )
            return result

        # calculate solvation free enthalpy
        result["gsolv"] = result["energy_solv"] - result["energy_gas"]        
        result["success"] = True

        """
        TODO - same thing as for xtb_gsolv
                self.job["energy"] = energy_gas
                self.job["energy2"] = energy_solv - energy_gas
                if self.job["trange"]:
                    tmp = {}
                    for temperature in self.job["trange"]:
                        tmp[temperature] = energy_solv - energy_gas
                    tmp[self.job["temperature"]] = energy_solv - energy_gas
                    self.job["erange1"] = tmp
                else:
                    self.job["erange1"] = {
                        self.job["temperature"]: energy_solv - energy_gas
                    }
                self.job["success"] = True
                """

        return result


    @QmProc._create_jobdir
    def _xtb_opt(self, conf: GeometryData, filename: str = "xtb_opt"):
        """
        ORCA geometry optimization using ANCOPT

        result = {
            "success": None,
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "mo_path": None,
            "geom": None,
        }
        """
        # NOTE: some "intuitivity problems":
        # the first call of _xtb_opt (probably in spearman opt) generates a coord file, which is then updated externally by xtb
        # the content of this coord file is converted into conf.xyz to be used by orca
        
        # prepare result
        # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
        # 'energy' contains the final energy of the optimization (converged or unconverged)
        # 'success' should only be False if the external program encounters an error
        # 'geom' stores the optimized geometry in GeometryData.xyz format
        result = {
            "success": None,
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
            "grad_norm": None,
            "mo_path": None,
            "geom": None,
        }

        jobdir = os.path.join(self.workdir, conf.name, "xtb_opt")
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
        conf.tocoord(os.path.join(jobdir, f"{filename}.coord"))
        
        # write xyz-file for orca
        conf.toxyz(os.path.join(jobdir, f"{filename}.xyz"))

        # set orca in path
        inputpath = os.path.join(jobdir, f"{filename}.inp")

        # prepare input dict
        parser = OrcaParser()
        indict = self.__prep(conf, filename, xyzfile=f"{filename}.xyz")
        
        # write orca input into file "xtb_opt.inp" in a subdir created for the conformer
        parser.write_input(inputpath, indict)

        # append some additional lines to the coord file for ancopt
        with open(
            os.path.join(jobdir, f"{filename}.coord"), "a", newline=None
        ) as newcoord:
            newcoord.writelines([
                    "$external\n",
                    f"   orca input file= {filename}.inp\n",
                    f"   orca bin= {self.instructions['orcapath']}\n",
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
                f"hlow={self.instructions['gradthr']} \n",
                "s6=30.00 \n",
                # remove unnecessary sp/gradient call in xTB
                "engine=lbfgs\n",
                "$external\n",
                f"   orca input file= {filename}.inp\n",
                f"   orca bin= {self.instructions['orcapath']} \n",
                "$end \n",
            ])

        # prepare xtb call
        call = [
            self.instructions["xtbpath"],
            f"{filename}.coord", # name of the coord file generated above
            "--opt",
            self.instructions["optlevel"],
            "--orca",
            "-I",
            xcontrolname
        ]
        
        # set path to the ancopt output file
        outputpath = os.path.join(jobdir, f"{filename}.out")

        print(f"Running optimization in {last_folders(jobdir, 2):18}")

        # call xtb
        returncode = self._make_call(call, outputpath, jobdir)

        # check if optimization finished without errors
        if returncode != 0:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}optimization "
                f"in {last_folders(self.workdir, 2):18} failed!"
            )
            return result

        # read output
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as file:
                lines = file.readlines()

                result["ecyc"] = []
                result["cycles"] = 0

                for line in lines:
                    # the following is probably self explanatory
                    if (
                        "external code error" in line
                        or "|grad| > 500, something is totally wrong!" in line
                        or "abnormal termination of xtb" in line
                    ):
                        # TODO - error handling
                        print(
                            f"{'ERROR:':{WARNLEN}}optimization in "
                            f"{last_folders(self.workdir, 2):18} failed!"
                        )
                        result["success"] = False
                        return result

                    result["success"] = True
                    
                    if " FAILED TO CONVERGE GEOMETRY " in line:
                        result["cycles"] += int(line.split()[7])
                        result["converged"] = False
                    elif "*** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
                        result["cycles"] += int(line.split()[5])
                        result["converged"] = True
                    elif "av. E: " in line and "->" in line:
                        try:
                            result["ecyc"].append(float(line.split("->")[-1]))
                        except ValueError as e:
                            # TODO - error handling
                            print(
                                f"{'ERROR:':{WARNLEN}}in {conf.name} calculation:\n{e}"
                            )
                            result["success"] = False
                            return result
                    elif " :: gradient norm      " in line:
                        try:
                            result["grad_norm"] = float(line.split()[3])
                        except ValueError as e:
                            # TODO - error handling
                            print(
                                f"{'ERROR:':{WARNLEN}}in {conf.name} calculation:\n{e}"
                            )
                            result["success"] = False
                            return result
        else:
            # TODO - error handling
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")
            result["success"] = False
            return result
        
        # store the final energy of the optimization in 'energy' 
        result["energy"] = result["ecyc"][-1]
        result["success"] = True

        # store the path to the current .gbw file for this conformer
        result["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

        # read out optimized geometry und update conformer geometry with this
        conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
        result["geom"] = conf.xyz

        try:
            assert result["converged"] is not None
        except AssertionError:
            # TODO - error handling
            # this should never happen
            result["success"] = False
            print(f"{'ERROR:':{WARNLEN}}in CONF{self.id} calculation\n{e}")

        return result