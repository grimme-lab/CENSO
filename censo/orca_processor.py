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

from censo.cfg import (
    CODING,
    ENVIRON,
    WARNLEN,
    editable_ORCA_input,
    PLANCK,
    C
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

        # TODO - externalize line appending?

        return lines


# TODO - keep output if any job fails
# TODO - __prep middle
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
    def __prep(self, conf: GeometryData, jobtype: str, grid: str, no_solv: bool = False) -> OrderedDict:
        """
        prepare an OrderedDict to be fed into the parser in order to write an input file
        for jobtype 'jobtype' (e.g. sp)

        grid: mandatory argument determining the grid settings chosen from 'gridsettings' dict

        TODO - NMR/OR/UVVis etc preparation steps
        TODO - externalize comp parameter setups
        """

        # check ORCA version
        orca5 = True if self.paths["orcaversion"].startswith("5") else False

        indict = OrderedDict()
        indict["main"] = []

        # grab func, basis
        func = self.dfa_settings.functionals.get(self.instructions["func"]).get("orca")
        basis = self.instructions.get("basis") 
        indict["main"].append(func)

        # TODO - b3lyp-3c removed
        if self.instructions["func"] not in self.dfa_settings.composites:
            indict["main"].append(basis)
        ####################### SET RI ###########################
            # set  RI def2/J,   RIJCOSX def2/J
            # this is only set if no composite DFA is used
            # settings for double hybrids
            if self.instructions["func"] in self.dfa_settings.doublehs:
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
            elif self.instructions["func"] in self.dfa_settings.hybrids:
                indict["main"].extend(["def2/J", "RIJCOSX"])
                if not orca5:
                    indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])
            
            # settings for (m)ggas
            elif self.instructions["func"] in self.dfa_settings.ggas:
                indict["main"].extend(["RI", "def2/J"])    

        ########################## SET GRID ############################ 

        """
        used in nmr calculation
        if self.job["moread"] is not None:
            #! MORead
            #%moinp "jobname2.gbw"
            orcainput["moread"] = self.job["moread"]"""

        # use 'grid' setting from function signature to quickly choose the grid settings 
        indict["main"].extend(self.__gridsettings[orca5][grid])
        
        ########################## DISPERSION ####################### TODO TODO TODO TODO TODO
        # add dispersion
        # dispersion correction information
        mapping = {
            "d3bj": "d3bj",
            "d3(0)": "D3ZERO",
            "d4": "d4",
            "nl": "NL",
        }

        disp = self.dfa_settings.functionals.get(self.instructions["func"]).get("disp")
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

        ###################### SOLVENT/GEOM ########################

        # set keywords for the selected solvent model
        if not self.instructions["gas-phase"] and not no_solv:
            if self.instructions["sm"] in ("smd"):
                orcainput["cpcm"] = {
                    "smd": ["true"],
                    "smdsolvent": [f"{self.solvents_dict['smd'][1]}"],
                }
            elif self.instuctions["sm"] == "cpcm":
                indict["main"].append(f"CPCM({self.solvents_dict['cpcm'][1]})")

        # unpaired, charge, and coordinates
        # by default xyzfile format is unused, coordinates written directly into input file
        indict["geom"] = {
            "def": ["xyz", self.instructions["charge"], self.instructions["unpaired"]+1],
            "coord": conf.toorca(),
        }
        """
        don't use this for now
        else:
            # xyz geometry
            geom, _ = t2x(self.workdir)
            orcainput["geom"] = {
                "def": ["xyz", self.instructions["charge"], self.instructions["unpaired"]+1]
                "coord": []
            }
        """
        return indict


    #@self.__prep TODO
    def _sp(self, conf: GeometryData, silent=False, filename="sp", no_solv: bool = False) -> Dict[str, Any]:
        """
        ORCA single-point calculation

        result = {
            "energy": None,
            "success": None,
        }
        """
        # set results
        result = {
            "energy": None,
            "success": None,
        }

        # create dir for conf
        confdir = os.path.join(self.workdir, conf.name)
        if os.path.isdir(confdir):
            # TODO - error handling warning? stderr?
            print(f"Folder {confdir} already exists. Potentially overwriting files.")
        elif os.system(f"mkdir {confdir}") != 0 and not os.path.isdir(confdir):
            print(f"Workdir for conf {conf.name} could not be created.")
            result["success"] = False
            return result

        # set in/out path
        inputpath = os.path.join(self.workdir, conf.name, f"{filename}.inp")
        outputpath = os.path.join(self.workdir, conf.name, f"{filename}.out")

        # prepare input dict
        parser = OrcaParser()
        indict = self.__prep(conf, "sp", "low+", no_solv=no_solv) # TODO - IMPORTANT not every sp should use low+ gridsize
        
        # write input into file "{filename}.inp" in a subdir created for the conformer
        parser.write_input(inputpath, indict)

        if not silent:
            print(f"Running single-point in {inputpath}")
        
        # start SP calculation
        with open(outputpath, "w", newline=None) as outputfile:
            # make external call to orca with "{filename}.inp" as argument
            call = [self.paths["orcapath"], f"{filename}.inp"]
            subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=os.path.join(self.workdir, conf.name),
                stdout=outputfile,
                env=ENVIRON,
            )

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
                print(
                    f"{'ERROR:':{WARNLEN}}scf in {inputpath} "
                    "not converged!"
                )
        else:
            # TODO - error handling
            result["success"] = False
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")
        
        return result


    def _gsolv(self, conf: GeometryData):
        """
        Calculate SMD_gsolv, needs ORCA
        if optimization is not performed with ORCA, only the density 
        functional for optimization is employed, 
        from my understanding smd is parametrized at 298 K, therefore it should only
        be used at this temperature.

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

        print(
            f"Running SMD_gsolv calculation in "
            f"{last_folders(self.workdir, 2)}."
        )

        # calculate gas phase
        res = self._sp(silent=True, filename="sp_gas", no_solv=True)

        if self.result["success"]:
            result["energy_gas"] = res["energy"]
        else:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}in gas phase single-point "
                f"of {last_folders(self.workdir, 2):18}"
            )
            return result

        # calculate in solution
        res = self._sp(silent=True, filename="sp_solv")

        if self.result["success"]:
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


    def _xtbopt(self, conf: GeometryData):
        """
        ORCA geometry optimization using ANCOPT
        implemented within xtb, generates inp.xyz, inp (orca-input) 
        and adds information to coord (xtb can then tell which file 
        orca has to use).

        result = {
            "success": None,
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
        }

        uses:
        fullopt --> outputname decision
        workdir --> folder of calculation

        TODO - condense this function
        """
        # prepare result
        # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
        # 'energy' contains the final energy of the optimization (converged or unconverged)
        result = {
            "success": None,
            "energy": None,
            "cycles": None,
            "converged": None,
            "ecyc": None,
        }
        
        # FIXME
        if self.job["fullopt"]:
            output = "opt-part2.out"
        else:
            output = "opt-part1.out"
        
        outputpath = os.path.join(self.workdir, output)
        
        print(f"Running optimization in {last_folders(self.workdir, 2):18}")
        files = [
            "xtbrestart",
            "xtbtopo.mol",
            "xcontrol-inp",
            "wbo",
            "charges",
            "gfnff_topo",
        ]
        
        # remove potentially preexisting files to avoid confusion
        for file in files:
            if os.path.isfile(os.path.join(self.workdir, file)):
                os.remove(os.path.join(self.workdir, file))
        
        # add inputfile information to coord (xtb as a driver)
        with open(
            os.path.join(self.workdir, "coord"), "r", newline=None
        ) as coord:
            tmp = coord.readlines()
        with open(
            os.path.join(self.workdir, "coord"), "w", newline=None
        ) as newcoord:
            for line in tmp[:-1]:
                newcoord.write(line)
            newcoord.write("$external\n")
            newcoord.write("   orca input file= inp\n")
            newcoord.write(
                f"   orca bin= {os.path.join(external_paths['orcapath'], 'orca')} \n"
            )
            newcoord.write("$end\n")

        # TODO - generate orca input file
        self.__prep(conf, "xtbopt")

        # generate coord file to be used by xtb
        conf.tocoord(os.path.join(self.workdir, "coord"))

        # prepare configuration file for ancopt (xcontrol file)
        with open(
            os.path.join(self.workdir, "opt.inp"), "w", newline=None
        ) as out:
            out.write("$opt \n")
            if (
                self.job["optcycles"] is not None
                and float(self.job["optcycles"]) > 0
            ):
                out.write(f"maxcycle={str(self.job['optcycles'])} \n")
                out.write(f"microcycle={str(self.job['optcycles'])} \n")

            out.writelines([
                "average conv=true \n",
                f"hlow={self.job.get('hlow', 0.01)} \n",
                "s6=30.00 \n",
                # remove unnecessary sp/gradient call in xTB
                "engine=lbfgs\n",
                "$external\n",
                "   orca input file= inp\n",
                    f"   orca bin= {self.paths['orcapath']} \n",
                "$end \n",
            ])

        # TODO - ???
        time.sleep(0.02)

        call = [
            self.paths["xtbpath"],
            "coord", # name of the coord file generated above
            "--opt",
            self.instructions["optlevel"],
            "--orca",
            "-I",
            "opt.inp", # name of the xcontrol file generated above
        ]

        # make xtb call, write into 'outputfile'
        with open(outputpath, "w", newline=None) as outputfile:
            returncode = subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
                env=ENVIRON,
            )

        # check if optimization finished correctly
        if returncode != 0:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}optimization "
                f"in {last_folders(self.workdir, 2):18} not converged"
            )
            return result

        # TODO - ???
        time.sleep(0.02)

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
                            f"{last_folders(self.workdir, 2):18} not converged"
                        )
                        result["success"] = False
                        return result
                    
                    if " FAILED TO CONVERGE GEOMETRY " in line:
                        self.result["cycles"] += int(line.split()[7])
                        self.result["converged"] = False
                    
                    if "*** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
                        self.result["cycles"] += int(line.split()[5])
                        self.result["converged"] = True
                    
                    if "av. E: " in line and "->" in line:
                        try:
                            self.result["ecyc"].append(float(line.split("->")[-1]))
                        except ValueError as e:
                            # TODO - error handling
                            print(
                                f"{'ERROR:':{WARNLEN}}in CONF{self.id} calculation:\n{e}"
                            )
                            result["success"] = False
                            return result
                    
                    if " :: gradient norm      " in line:
                        try:
                            self.result["grad_norm"] = float(line.split()[3])
                        except ValueError as e:
                            # TODO - error handling
                            print(
                                f"{'ERROR:':{WARNLEN}}in CONF{self.id} calculation:\n{e}"
                            )
                            result["success"] = False
                            return result
        else:
            # TODO - error handling
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")
            result["success"] = False
            return result
        
        # store the final energy of the optimization in 'energy' 
        # TODO - can an unconverged optimization be called a 'success'?
        self.result["energy"] = self.result["ecyc"][-1]
        self.result["success"] = True

        """if not self.job["onlyread"]:
            # convert optimized xyz to cinpoord file
            x2t(self.workdir, infile="inp.xyz")
        return """

        return result