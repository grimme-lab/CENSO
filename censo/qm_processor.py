"""
Contains QmProc base class, 
Additionally contains functions which should be present irrespective of the QM
code. (xTB always available)
"""
import os
from typing import Any, Callable, Dict, List, Union
from censo.datastructure import GeometryData

from censo.orca_processor import OrcaProc
from censo.tm_processor import TmProc
from math import isclose
import time
import subprocess
import json
from censo.cfg import (
    ENVIRON,
    CODING,
    WARNLEN,
    rot_sym_num,
)
from censo.utilities import last_folders, print, frange
        

class QmProc:
    """
    QmProc base class with xtb as driver (see _xtb methods)
    """

    def __init__(self, paths: Dict[str, str], solvents_dict: Dict = None):
        # jobtype is basically an ordered (!!!) (important for example if sp is required before the next job)
        # list containing the instructions of which computations to do        
        self.jobtype: List[str]
        
        # stores instructions, meaning the settings that should be applied for all jobs
        # e.g. 'gfnv' (GFN version for xtb_sp/xtb_rrho/xtb_gsolv)
        self.instructions: Dict[str, Any]
        
        # absolute path to the folder where jobs should be executed in
        self.workdir: str

        # dict to map the jobtypes to their respective methods
        self._jobtypes: Dict[str, Callable] = {
            "sp": self._sp,
            "gsolv": self._gsolv,
            "opt": self._opt,
            "genericoutput": self._genericoutput,
            "xtb_sp": self._xtb_sp,
            "xtb_gsolv": self._xtb_gsolv,
            "xto_opt": self._xtb_opt,
            "xtb_rrho": self._xtbrrho,
        }

        # stores the solvent lookup dict
        self.solvents_dict = solvents_dict

        # stores lookup dict for external paths
        self.paths: Dict[str, str] = paths


    def run(self, conformer: GeometryData) -> Dict[int, Dict[str, Any]]:
        """
        run methods depending on jobtype
        DO NOT OVERRIDE OR OVERLOAD! this will possibly break e.g. ProcessHandler.execute
        """
        
        res = {conformer.id: {}}
        
        for job in self.jobtype:
            res[conformer.id][job] = self._jobtypes[job](conformer)
            
        # returns dict e.g.: {140465474831616: {"sp": ..., "gsolv": ..., etc.}}
        return res


    def print(self):
        """
        print method for each part, should be implemented if needed
        """
        pass


    def _prep(self):
        """
        input preparation
        """
        pass


    def _sp(self, silent=False):
        """
        single-point calculation
        """
        pass

    def _opt(self):
        """
        geometry optimization
        """
        pass
    
    
    def _gsolv(self):
        """
        gsolv calculation using the respective solvent model
        """
        pass
    

    def _genericoutput(self):
        """
        Read shielding and coupling constants and write them to plain output
        The first natom lines contain the shielding constants, and from
        line natom +1 the coupling constants are written.
        """
        pass

    def _get_sym_num(self, sym=None, linear=None):
        """Get rotational symmetry number from Schoenfließ symbol"""
        if sym is None:
            sym = "c1"
        if linear is None:
            linear = False
        symnum = 1
        if linear and "c" in sym.lower()[0]:
            symnum = 1
            return symnum
        elif linear and "d" in sym.lower()[0]:
            symnum = 2
            return symnum
        for key in rot_sym_num.keys():
            if key in sym.lower():
                symnum = rot_sym_num.get(key, 1)
                break
        return symnum


    def _xtb_sp(self, filename: str = "xtb_sp.out", no_solv: bool = False, silent: bool = False) -> Dict[str, Any]:
        """
        Get single-point energy from xtb
        result = {
            "energy": None,
            "success": None,
        }
        """
        
        if not silent:
            print(
                f"Running xtb_sp calculation in "
                f"{last_folders(self.workdir, 3)}"
            )

        outputpath = os.path.join(self.workdir, filename)
        # if not self.job["onlyread"]: ??? TODO
        files = [
            "xtbrestart",
            "xtbtopo.mol",
            "xcontrol-inp",
            "wbo",
            "charges",
            "gfnff_topo",
            filename,
        ]

        # remove potentially preexisting files to avoid confusion
        for file in files:
            if os.path.isfile(os.path.join(self.workdir, file)):
                os.remove(os.path.join(self.workdir, file))

        # setup call for xtb single-point
        call = [
            self.paths["xtbpath"],
            "coord",
            "--" + self.instructions["gfnv"],
            "--sp",
            "--chrg",
            self.instructions["charge"],
            "--norestart",
            "--parallel",
            self.instructions["omp"],
        ]

        # add solvent to xtb call if not a gas-phase sp 
        # (set either through run settings or by call kwarg e.g. for _xtb_gsolv)
        # note on solvents_dict (or rather censo_solvents.json): 
        # [0] is the normal name of the solvent, if it is available, [1] is the replacement
        if not self.instructions.get("gas-phase", False) or no_solv:
            call.extend(
                [
                    "--" + self.instructions["sm_rrho"],
                    self.solvents_dict["xtb"][1], # auto-choose replacement solvent by default
                    "reference",
                ]
            )

            # set gbsa grid if gbsa is used, do nothing for alpb
            if self.instructions["sm_rrho"] == "gbsa":
                call.extend(
                    [
                        "-I",
                        "xcontrol-inp",
                    ]
                )
                with open(
                    os.path.join(self.workdir, "xcontrol-inp"), "w", newline=None
                ) as xcout:
                    xcout.write("$gbsa\n")
                    xcout.write("  gbsagrid=tight\n")
                    xcout.write("$end\n")

        # call xtb and write output into outputfile
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

        # set results
        result = {
            "energy": None,
            "success": None,
        }

        # if returncode != 0 then some error happened in xtb
        if returncode != 0:
            result["energy"] = None
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}{self.instructions['gfnv'].upper()} error in "
                f"{last_folders(self.workdir, 2)}"
            )
            return result

        # FIXME - ???
        time.sleep(0.02)

        # read energy from outputfile
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as outputfile:
                for line in outputfile.readlines():
                    if "| TOTAL ENERGY" in line:
                        try:
                            result["energy"] = float(line.split()[3])
                            result["success"] = True
                        except Exception:
                            if not silent:
                                print(
                                    f"{'ERROR:':{WARNLEN}}while converting "
                                    f"single-point in: {last_folders(self.workdir, 2)}"
                                )
                            result["energy"] = None
                            result["success"] = False
                            return result
        else:
            result["energy"] = None
            result["success"] = False
            if not silent:
                print(
                    f"{'ERROR:':{WARNLEN}}{self.instructions['gfnv'].upper()} error in "
                    f"{last_folders(self.workdir, 2)}"
                )
            return result

        # everything went fine, return result
        if not silent:
            print(
                f"{self.instructions['gfnv'].upper()} energy for {last_folders(self.workdir, 2)}"
                f" = {result['energy']:>.7f}"
            )

        return result


    def _xtb_gsolv(self) -> Dict[str, Any]:
        """
        Calculate additive GBSA or ALPB solvation contribution by
        Gsolv = Esolv - Egas, using GFNn-xTB or GFN-FF
        result = {
            "success": None,
            "gsolv": None,
            "energy_xtb_gas": None,
            "energy_xtb_solv": None,
        }
        """
        print(
            f"Running xtb_gsolv calculation in "
            f"{last_folders(self.workdir, 3)}"
        )
        tmp_gas = None
        tmp_solv = None

        # what is returned in the end
        result = {
            "success": None,
            "gsolv": None,
            "energy_xtb_gas": None,
            "energy_xtb_solv": None,
        }

        # run gas-phase GFN single-point
        res = self._xtb_sp(filename="gas.out", silent=True, no_solv=True)
        if res["success"]:
            tmp_gas = res["energy"]
        else:
            print(
                f"{'ERROR:':{WARNLEN}}Gas phase {self.instructions['gfnv'].upper()} error in "
                f"{last_folders(self.workdir, 3)}"
            )
            result["success"] = False
            return result

        # run single-point in solution:
        # ''reference'' corresponds to 1\;bar of ideal gas and 1\;mol/L of liquid
        #   solution at infinite dilution,
        res = self._xtb_sp(filename="solv.out", silent=True)
        if res["success"]:
            tmp_solv = res["energy"]
        else:
            print(
                f"{'ERROR:':{WARNLEN}}Solution phase {self.instructions['gfnv'].upper()} error in "
                f"{last_folders(self.workdir, 3)}"
            )
            result["success"] = False
            return result

        # only reached if both gas-phase and solvated sp succeeded   
        result["gsolv"] = tmp_solv - tmp_gas

        # TODO - what is the sense behind this?
        # leave this out for now, only 'calculate' this when needed
        """if self.instructions.get("trange", None):
            result["erange1"] = {
                temperature: tmp_solv - tmp_gas 
                for temperature in frange(self.instructions["trange"][0], self.instructions["trange"][1], self.instructions["trange"][2])
            }
        else:
            self.job["erange1"] = None"""

        result["success"] = True
        result["energy_xtb_gas"] = tmp_gas
        result["energy_xtb_solv"] = tmp_solv

        return result


    def _xtb_opt(self):
        """
        geometry optimization using xtb as driver, has to be implemented for each qm code
        """
        pass


    def _xtbrrho(self, filename="hess.out"):
        """
        mRRHO contribution with GFNn-xTB/GFN-FF
        result = {
            "energy": None,
            "success": None,
            "rmsd": None,
            "gibbs": None,
            "enthalpy": None,
            "entropy": None,
            "symmetry": None,
            "symnum": None,
        }
        """
        outputpath = os.path.join(self.workdir, filename)
        xcontrolpath = os.path.join(self.workdir, "xcontrol-inp")
        # if not self.job["onlyread"]: # TODO ???
        print(
            f"Running {str(self.instructions['gfnv']).upper()} mRRHO in "
            f"{last_folders(self.workdir, 2)}"
        )

        # TODO - is this list complete?
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

        # if a temperature in the trange is close (+-0.6) to the fixed temperature then replace this value in trange
        # don't do this for now, to make the program more predictable
        """ for t in trange:
            if isclose(self.instructions["temperature"], t, abs_tol=0.6):
                trange[trange.index(t)] = self.instructions["temperature"]
                break """

        # setup xcontrol
        with open(xcontrolpath, "w", newline=None) as xcout:
            xcout.write("$thermo\n")
            if self.instructions.get("trange", False):
                trange = frange(self.instructions['trange'][0], self.instructions['trange'][1], self.instructions['trange'][2])
                xcout.write(
                    f"    temp="
                    f"{','.join([str(i) for i in trange])}\n")
            else:
                xcout.write(f"    temp={self.instructions['temperature']}\n")


            xcout.write(f"    sthr={self.instructions['sthr']}\n")

            xcout.write(f"    imagthr={self.instructions['imagthr']}\n")

            # TODO - when do you actually want to set the scale manually?
            # don't set scale manually for now
            """ if self.instructions.get("scale", "automatic") != "automatic":
                # for automatic --> is method dependant leave it to xTB e.g. GFNFF has a
                # different scaling factor than GFN2
                xcout.write(f"    scale={self.instructions['scale']}\n") """
            
            xcout.write("$symmetry\n")
            xcout.write("     maxat=1000\n")
            # always consider symmetry
            # xcout.write("    desy=0.1\n") # taken from xtb defaults
            # xcout.write("    desy=0.0\n")

            # set gbsa grid
            if not self.instructions["gas-phase"] and self.instructions["sm_rrho"] == "gbsa":
                xcout.write("$gbsa\n")
                xcout.write("  gbsagrid=tight\n")
            
            xcout.write("$end\n")

        if self.instructions["bhess"]:
            # set ohess or bhess
            dohess = "--bhess"
            olevel = "vtight"
        else:
            dohess = "--ohess"
            olevel = "vtight"

        # FIXME ???
        time.sleep(0.02)
        
        with open(outputpath, "w", newline=None) as outputfile:
            call = [
                self.paths["xtbpath"],
                "coord",
                "--" + self.instructions["gfnv"],
                dohess,
                olevel,
                "--chrg",
                self.instructions["charge"],
                "--enso",
                "--norestart",
                "-I",
                "xcontrol-inp",
                "--parallel",
                self.instructions["omp"],
            ]
            if not self.instructions["gas-phase"]:
                call.extend(
                    [
                        "--" + self.instructions["sm_rrho"],
                        self.solvents_dict["xtb"][1], # auto-choose replacement solvent by default
                    ]
                )
            if self.instructions["rmsdbias"]:
                # move one dir up to get to cwd (FIXME)
                cwd = os.path.join(os.path.split(self.workdir)[::-1][1:][::-1])

                call.extend(
                    [
                        "--bias-input",
                        os.path.join(cwd, "rmsdpot.xyz"),
                    ]
                )

            # run xtb
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

        # what is returned in the end
        result = {
            "energy": None,
            "success": None,
            "rmsd": None,
            "gibbs": None,
            "enthalpy": None,
            "entropy": None,
            "symmetry": None,
            "symnum": None,
        }

        # check if converged:
        if returncode != 0:
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}{self.instructions['gfnv'].upper()} ohess/bhess error in "
                f"{last_folders(self.workdir, 2):18}"
            )
            return result

        # start reading output
        # TODO - error handling
        if not os.path.isfile(outputpath):
            result["success"] = False
            print(
                f"{'ERROR:':{WARNLEN}}file {self.workdir}/hess.out could not be found!"
            )
            return result

        # read output and store lines
        with open(outputpath, "r", encoding=CODING, newline=None) as outputfile:
            lines = outputfile.readlines()

        # get gibbs energy, enthalpy and entropy for given temperature range
        # TODO - there should be no case where the trange setting is unset
        trange = frange(self.instructions["trange"][0], self.instructions["trange"][1], self.instructions["trange"][2])
        # gibbs energy
        gt = {}

        # enthalpy
        ht = {}
        
        # rotational entropy
        rotS = {}

        for line in lines:
            # get rotational entropy
            if "VIB" in line:
                try:
                    index = lines.index(line) + 1
                    T = float(line.split()[0])
                    rotS[T] = float(lines[index].split()[4])
                except (KeyError, ValueError):
                    # TODO - error handling !!!
                    pass

            # get gibbs energy and enthalpy
            if "T/K" in line:
                for line2 in lines[lines.index(line)+2:]:
                    if "----------------------------------" in line2:
                        break
                    else:
                        try:
                            T = float(line2.split()[0])
                            gt[T] = float(line2.split()[4])
                            ht[T] = float(line2.split()[2])
                        except (ValueError, KeyError):
                            # TODO - error handling
                            print(f"{'ERROR:':{WARNLEN}}can not convert G(T)")

            # extract rmsd
            if "final rmsd / " in line and self.instructions["bhess"]:
                try:
                    result["rmsd"] = float(line.split()[3])
                except (ValueError, IndexError):
                    # TODO - error handling ?
                    result["rmsd"] = None

            # extract symmetry
            if ":  linear? " in line:
                # linear needed for symmetry and S_rot (only if considersym is turned off)
                try:
                    if line.split()[2] == "false":
                        result["linear"] = False
                    elif line.split()[2] == "true":
                        result["linear"] = True
                except (IndexError, Exception) as e:
                    # TODO - error handling
                    print(e)

        # check if xtb calculated the temperature range correctly
        if len(trange) == len(gt) and len(trange) == len(ht) and len(trange) == len(rotS):
            result["gibbs"] = gt
            result["enthalpy"] = ht
            result["entropy"] = rotS
        else:
            result["success"] = False
            return result

        # xtb_enso.json is generated by xtb by using the '--enso' argument *only* when using --bhess or --ohess (when a hessian is calculated)
        # contains output from xtb in json format to be more easily digestible by CENSO
        if os.path.isfile(os.path.join(self.workdir, "xtb_enso.json")):
            with open(
                os.path.join(self.workdir, "xtb_enso.json"),
                "r",
                encoding=CODING,
                newline=None,
            ) as f:
                data = json.load(f)

            # read number of imaginary frequencies and print warning
            if "number of imags" in data:
                if data["number of imags"] > 0:
                    print(
                        f"{'WARNING:':{WARNLEN}}found {data['number of imags']} significant"
                        f" imaginary frequencies in "
                        f"{last_folders(self.workdir, 2)}"
                    )

            # get gibbs energy
            if "G(T)" in data:
                # TODO - error handling
                if float(self.instructions["temperature"]) == 0:
                    result["success"] = True
                    result["energy"] = data.get("ZPVE", 0.0)
                    result["gibbs"][self.instructions["temperature"]] = data.get("ZPVE", 0.0)
                    result["enthalpy"][self.instructions["temperature"]] = data.get("ZPVE", 0.0)
                    result["entropy"][self.instructions["temperature"]] = None # set this to None for predictability
                else:
                    result["success"] = True
                    result["energy"] = data.get("G(T)", 0.0)
                    result["gibbs"][self.instructions["temperature"]] = data.get("G(T)", 0.0)
                    result["enthalpy"][self.instructions["temperature"]] = None # set this to None for predictability
                    result["entropy"][self.instructions["temperature"]] = None # set this to None for predictability

                # only determine symmetry if all the needed information is there
                if "point group" and "linear" in data.keys():
                    result["symmetry"] = data["point group"]
                    result["linear"] = data.get("linear", result["linear"])
                    # calculate symnum
                    result["symnum"] = self._get_sym_num(
                        sym=result["symmetry"], linear=result["linear"]
                    )
                else:
                    """could not determine symnum"""
                    # TODO
            else:
                # TODO - error handling
                print(
                    f"{'ERROR:':{WARNLEN}}while reading xtb_enso.json in: {last_folders(self.workdir, 2)}"
                )
                result["success"] = False
        else:
            # TODO - error handling
            print(
                f"{'WARNING:':{WARNLEN}}File "
                f"{os.path.join(self.workdir, 'xtb_enso.json')} doesn't exist!"
            )
            result["success"] = False


class ProcessorFactory:
    
    # for now these are the only available processor types
    __proctypes: Dict[str, type] = {
        "orca": OrcaProc,
        "tm": TmProc,
    }
    
    @classmethod
    def create_processor(cls, *args, **kwargs) -> QmProc:
        """
        returns an instance of the requested processor type (mapping with the 'prog' setting)
        for now the QmProc class uses xtb as driver (this method should be changed if additional drivers are implemented)
        available processor types are mapped in ProcessorFactory.__proctypes

        example: create_processor(orca, external_paths, solvents_dict=...) 
        (lookup the constructors for the processor types for further documentation)
        """
        type_t: type = cls.__proctypes.get(prog, None)
        
        if not type_t is None:
            return type_t(*args, **kwargs)
        else:
            raise TypeError(f"No processor type was found for {prog} in {cls.__proctypes}.")