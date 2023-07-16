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
        "some setting after geom": {"some option": [...], ...}

        comments get removed from the read input lines and are therefore lost

        TODO - cannot deal with %coord or internal coordinates
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
                lines.append(f"    {option} {reduce(lambda x, y: x + ' ' + y, indict[key][option])}\n")
            lines.append("end\n")
        
        # next, write the geometry input lines
        # geometry definition line (e.g. "* xyzfile 0 1 input.xyz" / "* xyz 0 1")
        lines.append(f"* {reduce(lambda x, y: x + ' ' + y, indict['geom']['def'])}\n")
        
        # write coordinates if "xyz" keyword is used
        if indict["geom"].get("coord", False):
            for coord in indict["geom"]["coord"]:
                lines.append(f"{reduce(lambda x, y: x + ' ' + y, coord)}\n")
            lines.append("*\n")

        # lastly, write all the keywords and options that should be placed after the geometry input (e.g. NMR settings)
        for key in allkeys[allkeys.index("geom")+1:]:
            lines.append(f"%{key}\n")
            for option in indict[key].keys():
                lines.append(f"    {option} {reduce(lambda x, y: x + ' ' + y, indict[key][option])}\n")
            lines.append("end\n")

        # TODO - externalize line appending?

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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # expand jobtypes with special orca jobtypes
        self._jobtypes = {
            **self.__jobtypes, **{
                "nmrS": self.__nmrS,
                "nmrJ": self.__nmrJ,
                "uvvis": self.__uvvis,
            }
        }


    def __prep(self, jobtype: str, xyzfile: str = None, no_solv: bool = False) -> OrderedDict:
        """
        prepare an OrderedDict to be fed into the parser in order to write an input file
        for jobtype 'jobtype' (e.g. sp)

        TODO - NMR/OR/UVVis etc preparation steps
        TODO - support for xyz preparation?
        TODO - externalize comp parameter setups
        """
        indict = OrderedDict()
        indict["main"] = []

        # grab func, basis, solvation, computational parameters
        indict["main"].append(instructions["func"])
        indict["main"].append(instructions["basis"])

        # FIXME
        if instructions["func"] in self.dfa_settings.composites:
            if self.job["func"] == "b3lyp-3c":
                orcainput["functional"] = [f"! {'b3lyp'}"]
                orcainput["basis"] = [f"! {self.job['basis']}"]
                orcainput["gcp"] = [f"! GCP(DFT/SV(P))"]
            else:
                orcainput["functional"] = [
                    f"! {dfa_settings.composites.get(self.job['func']).get('orca')}"
                ]
        else:
        ####################### SET RI/GRID ###########################
        # if not composite method
            # set  RI def2/J,   RIJCOSX def2/J gridx6 NOFINALGRIDX,  RIJK def2/JK
            # settings for double hybrids
            if self.job["func"] in dfa_settings().dh_dfa():
                indict["main"].extend(["def2/J", "RIJCOSX"])

                if "nmr" in jobtype:
                    indict["main"].append("NOFROZENCORE")
                    indict["mp2"] = {
                        "RI": ["true"],
                        "density": ["relaxed"],
                    }
                else:
                    indict["main"].append("frozencore")
                    orcainput["mp2"] = ["%mp2", "    RI true", "end"]

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
            # TODO - what about b3lyp-3c??
            elif self.instructions["func"] in self.dfa_settings.hybrids:
                indict["main"].extend(["def2/J", "RIJCOSX"])
                if not orca5:
                    indict["main"].extend(["GRIDX6", "NOFINALGRIDX"])
            
            # settings for (m)ggas
            elif self.instructions["func"] in self.dfa_settings.ggas:
                indict["main"].extend(["RI", "def2/J"])    

        ########################## SET GRID ############################

        # set grid
        if (
            self.job["func"] in dfa_settings().dh_dfa()
            or self.job["func"] in dfa_settings().hybrid_dfa()
        ):
            if orca5:
                orcainput["grid"] = ["! DEFGRID2"]
            else:
                orcainput["grid"] = ["! grid5 nofinalgrid"]
        else:
            if orca5:
                orcainput["grid"] = ["! DEFGRID1"]
            else:
                orcainput["grid"] = ["! grid4 nofinalgrid"]
        if self.job["moread"] is not None:
            #! MORead
            #%moinp "jobname2.gbw"
            orcainput["moread"] = self.job["moread"]
        orcainput["scfconv"] = ["! scfconv6"]

        extension = {
            "low": {"grid": ["! grid4 nofinalgrid"], "scfconv": ["! loosescf"]},
            "low+": {"grid": ["! grid4 nofinalgrid"], "scfconv": ["! scfconv6"]},
            "high": {"grid": ["! grid4 nofinalgrid"], "scfconv": ["! scfconv7"]},
            "high+": {"grid": ["! grid5 nofinalgrid"], "scfconv": ["! scfconv7"]},
        }
        extension5 = {
            "low": {"grid": ["! DEFGRID2"], "scfconv": ["! loosescf"]},
            "low+": {"grid": ["! DEFGRID2"], "scfconv": ["! scfconv6"]},
            "high": {"grid": ["! DEFGRID2"], "scfconv": ["! scfconv7"]},
            "high+": {"grid": ["! DEFGRID2"], "scfconv": ["! scfconv7"]},
        }
        if self.job["prepinfo"]:
            if isinstance(self.job["prepinfo"], list):
                if self.job["prepinfo"][0] in extension.keys():
                    if orca5:
                        orcainput["grid"] = extension5[self.job["prepinfo"][0]]["grid"]
                        orcainput["scfconv"] = extension5[self.job["prepinfo"][0]]["scfconv"]
                    else:
                        orcainput["grid"] = extension[self.job["prepinfo"][0]]["grid"]
                        orcainput["scfconv"] = extension[self.job["prepinfo"][0]]["scfconv"]

        
        ########################## DISPERSION #######################
        # add dispersion
        # dispersion correction information
        if dfa_settings.functionals.get(self.job["func"]).get("disp") == "composite":
            if self.job["func"] == "b3lyp-3c":
                orcainput["disp"] = ["! d3bj"]
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "d3bj":
            if self.job["func"] not in ("b97-d3(0)", "b97-d3"):
                orcainput["disp"] = ["! d3bj"]
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "d3(0)":
            if self.job["func"] not in ("b97-d3(0)", "b97-d3"):
                orcainput["disp"] = ["! D3ZERO"]
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "d4":
            orcainput["disp"] = ["! D4"]
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "nl":
            # b3lyp NL
            if orca5:
                orcainput["disp"] = ["! NL "]
            else:
                orcainput["disp"] = ["! NL vdwgrid3"]
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "novdw":
            pass
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "included":
            pass
        else:
            print(
                f" {dfa_settings.functionals.get(self.job['func']).get('disp')} unknown dispersion option!"
            )


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
        if xyzfile is not None:
            indict["geom"] = {
                "def": ["xyzfile", self.instructions["charge"], self.instructions["unpaired"]+1, xyzfile]
            }
        else:
            raise Exception("no xyzfile") # TODO - error handling
        """
        don't use this for now
        else:
            # xyz geometry
            geom, _ = t2x(self.job["workdir"])
            orcainput["geom"] = {
                "def": ["xyz", self.instructions["charge"], self.instructions["unpaired"]+1]
                "coord": []
            }
        """
        return indict


    #@self.__prep TODO
    def _sp(self, silent=False, filename="sp.out", no_solv: bool = False) -> Dict[str, Any]:
        """
        ORCA single-point calculation

        result = {
            "energy": None,
            "success": None,
        }
        """
        # set output path
        outputpath = os.path.join(self.workdir, filename)

        # prepare input dict
        parser = OrcaParser()
        indict = self.__prep("sp", no_solv=no_solv)
        
        # write input into file "inp"
        parser.write_input(os.path.join(self.workdir, "inp"), indict)

        if not silent:
            print(f"Running single-point in {last_folders(self.workdir, 2)}")
        
        # start SP calculation
        with open(outputpath, "w", newline=None) as outputfile:
            # make external call to orca with "inp" as argument
            call = [self.paths["orcapath"], "inp"]
            subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
                env=ENVIRON,
            )

        # set results
        result = {
            "energy": None,
            "success": None,
        }

        # read output
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as inp:
                lines = inp.readlines()
                for i, line in enumerate(lines):
                    if "FINAL SINGLE POINT ENERGY" in line:
                        result["energy"] = float(line.split()[4])

                    # check if scf is converged:
                    if "ORCA TERMINATED NORMALLY" in line:
                        result["success"] = True
                    
            if not self.job["success"]:
                result["success"] = False
                print(
                    f"{'ERROR:':{WARNLEN}}scf in {last_folders(self.workdir, 2)} "
                    "not converged!"
                )
        else:
            result["success"] = False
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")
        
        return result


    def _gsolv(self):
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
        print(
            f"Running SMD_gsolv calculation in "
            f"{last_folders(self.workdir, 2)}."
        )

        # what is returned in the end
        result = {
            "success": None,
            "gsolv": None,
            "energy_gas": None,
            "energy_solv": None,
        }

        # calculate gas phase
        res = self._sp(silent=True, filename="sp_gas.out", no_solv=True)

        if self.job["success"]:
            result["energy_gas"] = res["energy"]
        else:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}in gas phase single-point "
                f"of {last_folders(self.workdir, 2):18}"
            )
            return result

        # calculate in solution
        res = self._sp(silent=True, filename="sp_solv.out")

        if self.job["success"]:
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


    def _xtbopt(self):
        """
        ORCA input generation and geometry optimization using ANCOPT
        implemented within xtb, generates inp.xyz, inp (orca-input) 
        and adds information to coord (xtb can then tell which file 
        orca has to use).

        uses:
        fullopt --> outputname decision
        workdir --> folder of calculation

        return:
        cycles --> number of optimization cycles
        ecyc --> energy at cycle
        energy --> energy at last step
        success --> calulation without crash
        converged --> geometry optimization converged
        """
        error_logical = False
        if self.job["fullopt"]:
            output = "opt-part2.out"
        else:
            output = "opt-part1.out"
        outputpath = os.path.join(self.job["workdir"], output)
        if not self.job["onlyread"]:
            print(f"Running optimization in {last_folders(self.job['workdir'], 2):18}")
            files = [
                "xtbrestart",
                "xtbtopo.mol",
                "xcontrol-inp",
                "wbo",
                "charges",
                "gfnff_topo",
            ]
            for file in files:
                if os.path.isfile(os.path.join(self.job["workdir"], file)):
                    os.remove(os.path.join(self.job["workdir"], file))
            # convert coord to xyz, write inp.xyz
            t2x(self.job["workdir"], writexyz=True, outfile="inp.xyz")
            # add inputfile information to coord (xtb as a driver)
            with open(
                os.path.join(self.job["workdir"], "coord"), "r", newline=None
            ) as coord:
                tmp = coord.readlines()
            with open(
                os.path.join(self.job["workdir"], "coord"), "w", newline=None
            ) as newcoord:
                for line in tmp[:-1]:
                    newcoord.write(line)
                newcoord.write("$external\n")
                newcoord.write("   orca input file= inp\n")
                newcoord.write(
                    f"   orca bin= {os.path.join(external_paths['orcapath'], 'orca')} \n"
                )
                newcoord.write("$end\n")

            with open(
                os.path.join(self.job["workdir"], "inp"), "w", newline=None
            ) as inp:
                for line in self._prep_input(xyzfile="inp.xyz"):
                    inp.write(line + "\n")
            # Done writing input!
            callargs = [
                external_paths["xtbpath"],
                "coord",
                "--opt",
                self.job["optlevel"],
                "--orca",
                "-I",
                "opt.inp",
            ]
            with open(
                os.path.join(self.job["workdir"], "opt.inp"), "w", newline=None
            ) as out:
                out.write("$opt \n")
                if (
                    self.job["optcycles"] is not None
                    and float(self.job["optcycles"]) > 0
                ):
                    out.write(f"maxcycle={str(self.job['optcycles'])} \n")
                    out.write(f"microcycle={str(self.job['optcycles'])} \n")
                out.write("average conv=true \n")
                out.write(f"hlow={self.job.get('hlow', 0.01)} \n")
                out.write("s6=30.00 \n")
                # remove unnecessary sp/gradient call in xTB
                out.write("engine=lbfgs\n")
                out.write("$external\n")
                out.write("   orca input file= inp\n")
                out.write(
                    f"   orca bin= {os.path.join(external_paths['orcapath'], 'orca')} \n"
                )
                out.write("$end \n")
            time.sleep(0.02)
            with open(outputpath, "w", newline=None) as outputfile:
                returncode = subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
            if returncode != 0:
                error_logical = True
                print(
                    f"{'ERROR:':{WARNLEN}}optimization "
                    f"in {last_folders(self.job['workdir'], 2):18} not converged"
                )
            time.sleep(0.02)
        # read output
        # check if optimization finished correctly:
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as inp:
                stor = inp.readlines()
                for line in stor:
                    if (
                        "external code error" in line
                        or "|grad| > 500, something is totally wrong!" in line
                        or "abnormal termination of xtb" in line
                    ):
                        print(
                            f"{'ERROR:':{WARNLEN}}optimization in "
                            f"{last_folders(self.job['workdir'], 2):18} not converged"
                        )
                        error_logical = True
                        break
                    elif " FAILED TO CONVERGE GEOMETRY " in line:
                        self.job["cycles"] += int(line.split()[7])
                        self.job["converged"] = False
                    elif "*** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
                        self.job["cycles"] += int(line.split()[5])
                        self.job["converged"] = True
                ################
                for line in stor:
                    if "av. E: " in line and "->" in line:
                        try:
                            self.job["ecyc"].append(float(line.split("->")[-1]))
                        except ValueError as e:
                            error_logical = True
                            print(
                                f"{'ERROR:':{WARNLEN}}in CONF{self.id} calculation:\n{e}"
                            )
                            break
                    if " :: gradient norm      " in line:
                        try:
                            self.job["grad_norm"] = float(line.split()[3])
                        except ValueError as e:
                            error_logical = True
                            print(
                                f"{'ERROR:':{WARNLEN}}in CONF{self.id} calculation:\n{e}"
                            )
                            break
        else:
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")
            error_logical = True
        if not error_logical:
            try:
                self.job["energy"] = self.job["ecyc"][-1]
                self.job["success"] = True
            except Exception:
                error_logical = True
        if error_logical:
            self.job["energy"] = 0.0
            self.job["success"] = False
            self.job["converged"] = False
            self.job["ecyc"] = []
            self.job["grad_norm"] = 10.0

        if not self.job["onlyread"]:
            # convert optimized xyz to coord file
            x2t(self.job["workdir"], infile="inp.xyz")
        return

    def _nmrS(self):
        """
        ORCA NMR shielding constant calculation
        """
        if not self.job["onlyread"]:
            with open(
                os.path.join(self.job["workdir"], "inpS"), "w", newline=None
            ) as inp:
                for line in self._prep_input():
                    inp.write(line + "\n")
            # Done input!
            # shielding calculation
            print(
                "Running shielding calculation in {:18}".format(
                    last_folders(self.job["workdir"], 2)
                )
            )
            with open(
                os.path.join(self.job["workdir"], "orcaS.out"), "w", newline=None
            ) as outputfile:
                call = [os.path.join(external_paths["orcapath"], "orca"), "inpS"]
                subprocess.call(
                    call,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
            time.sleep(0.1)
        # check if calculation was successfull:
        with open(
            os.path.join(self.job["workdir"], "orcaS.out"),
            "r",
            encoding=CODING,
            newline=None,
        ) as inp:
            store = inp.readlines()
            self.job["success"] = False
            for line in store:
                if "ORCA TERMINATED NORMALLY" in line:
                    self.job["success"] = True
        if not self.job["success"]:
            print(
                f"{'ERROR:':{WARNLEN}}shielding calculation in "
                f"{last_folders(self.job['workdir'], 1):18} failed!"
            )
        return

    def _nmrJ(self):
        """
        ORCA NMR coupling constant calculation

        uses:
        prepinfo nmrJ
        workdir
        progpath
        success
        """
        if not self.job["onlyread"]:
            # generate input   # double hybrids not implemented
            with open(
                os.path.join(self.job["workdir"], "inpJ"), "w", newline=None
            ) as inp:
                for line in self._prep_input():
                    inp.write(line + "\n")
            # Done input!
            # start coupling calculation
            print(
                "Running coupling calculation in {}".format(
                    last_folders(self.job["workdir"], 2)
                )
            )
            with open(
                os.path.join(self.job["workdir"], "orcaJ.out"), "w", newline=None
            ) as outputfile:
                call = [os.path.join(external_paths["orcapath"], "orca"), "inpJ"]
                subprocess.call(
                    call,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
            time.sleep(0.1)
        # check if calculation was successfull:
        with open(
            os.path.join(self.job["workdir"], "orcaJ.out"),
            "r",
            encoding=CODING,
            newline=None,
        ) as inp:
            store = inp.readlines()
            self.job["success"] = False
            for line in store:
                if "ORCA TERMINATED NORMALLY" in line:
                    self.job["success"] = True
        if not self.job["success"]:
            print(
                f"{'ERROR:':{WARNLEN}}coupling calculation "
                f"in {last_folders(self.job['workdir'], 1):18} failed!"
            )
        return

    def _genericoutput(self):
        """
        ORCA read shielding and coupling constants and write them to plain output
        """
        fnameshield = "orcaS.out"
        atom = []
        sigma = []
        try:
            with open(
                os.path.join(self.job["workdir"], fnameshield),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                data = inp.readlines()
            for line in data:
                if "CHEMICAL SHIELDING SUMMARY (ppm)" in line:
                    start = data.index(line)
            for line in data[(start + 6) :]:
                splitted = line.split()
                if len(splitted) == 4:
                    atom.append(int(splitted[0]) + 1)
                    sigma.append(float(splitted[2]))
                else:
                    break
        except FileNotFoundError:
            print(
                f"{'INFORMATION:':{WARNLEN}}Missing file "
                f"{fnameshield} in {last_folders(self.job['workdir'], 2)}"
            )
            self.job["success"] = False
        self.job["success"] = True
        fnamecoupl = "orcaJ.out"
        atom1 = []
        atom2 = []
        jab = []
        try:
            with open(
                os.path.join(self.job["workdir"], fnamecoupl),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                data = inp.readlines()
            for line in data:
                if "NMR SPIN-SPIN COUPLING CONSTANTS" in line:
                    start = int(data.index(line)) + 6
                if " ****ORCA TERMINATED NORMALLY****" in line:
                    end = int(data.index(line))

            for line in data[start:end]:
                if "NUCLEUS" in line:
                    tmpsplitted = line.split()
                    atom1.append(int(tmpsplitted[4]) + 1)
                    atom2.append(int(tmpsplitted[9]) + 1)
                elif "Total" in line and "iso= " in line:
                    splitted = line.split()
                    jab.append(float(splitted[5]))
                else:
                    pass
        except FileNotFoundError:
            print(
                f"{'INFORMATION:':{WARNLEN}}Missing file "
                f"{fnamecoupl} in {last_folders(self.job['workdir'], 2)}"
            )
            self.job["success"] = False
        self.job["success"] = True
        with open(
            os.path.join(self.job["workdir"], "nmrprop.dat"), "w", newline=None
        ) as out:
            s = sorted(zip(atom, sigma))
            atom, sigma = map(list, zip(*s))
            self.shieldings = dict(zip(atom, sigma))
            for i in range(len(atom)):
                out.write("{:{digits}} {}\n".format(atom[i], sigma[i], digits=4))
            for i in range(self.job["nat"] - len(atom)):
                out.write("\n")
            for i in range(len(atom1)):
                out.write(
                    "{:{digits}} {:{digits}}   {}\n".format(
                        atom1[i], atom2[i], jab[i], digits=4
                    )
                )
        time.sleep(0.02)
        return


    def _uvvis(self):
        """
        calculation of uvvis spectra
        """
        # read uvvis excitation energies and oscillator strengths
        if self.job["calc_uvvis"] and "ABSORPTION SPECTRUM" in line:
            # offset of 5 because of text in orca output
            # TODO - catch error in reading uvvis data
            uvvis_table = stor[i+5:i+5+self.job["nroots"]]
            for row in uvvis_table:
                tmp = {"wavelength": 0.0, "energy": 0.0, "osc_str": 0.0}
                tmp["energy"] = float(row.split()[1]) * PLANCK * C
                tmp["wavelength"] = float(row.split()[2])
                tmp["osc_str"] = float(row.split()[3])
                self.job["excitations"].append(tmp)