"""
Contains OrcaProc class for calculating ORCA related properties of conformers.
"""
from collections import OrderedDict
import os
import sys
import time
import subprocess
from typing import Any, Dict

from censo.cfg import (
    CODING,
    ENVIRON,
    WARNLEN,
    external_paths,
    dfa_settings,
    editable_ORCA_input,
    PLANCK,
    C
)
from censo.utilities import last_folders, t2x, x2t, print
from censo.qm_processor import QmProc

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

    def __init__(self):
        super().__init__()
        
        # expand jobtypes with special orca jobtypes
        self.jobtypes = {
            **self.jobtypes, **{
                "nmrS": self._nmrS,
                "nmrJ": self._nmrJ,
                "uvvis": self._uvvis,
            }
        }
        
        
    def _prep(self, xyzfile=False, returndict=False):
        """
        cefine preparation step analogue

        use:
        xyzfile --> if set : * xyzfile ...  xyzfile.xyz 
        returns
        call --> list with settings
        """

        # understood commands in prepinfo:
        # DOGCP = uses gcp if basis known for gcp

        orcainput_start = OrderedDict(
            [
                ("moread", None),
                ("functional", None),
                ("disp", None),
                ("basis", None),
                ("gcp", None),
                ("RI-approx", None),
                ("grid", None),
                ("scfconv", None),
                ("frozencore", None),
                ("mp2", None),
                (
                    "default",
                    editable_ORCA_input.get('default',
                    [
                        "! smallprint printgap noloewdin",
                        "! NOSOSCF",
                        "%MaxCore 8000",
                        "%output",
                        "       print[P_BondOrder_M] 1",
                        "       print[P_Mayer] 1",
                        "       print[P_basis] 2",
                        "end",
                    ]),
                ),
                ("job", None),
                ("optthreshold", None),
                ("parallel", None),
                ("solvation", None),
                ("geom", None),
                ("couplings", None),
                ("shieldings", None),
                ("uvvis", None),
                ("uvvis_nroots", None)
            ]
        )
        if "nmrJ" or "nmrS" in self.job["prepinfo"]:
            nmrprop = True
        else:
            nmrprop = False


        # check ORCA5 or older versions:
        try:
            if int(external_paths['orcaversion'].split('.')[0]) >= 5:
                orca5 = True
            else:
                orca5 = False
        except:
            print(f"{'ERROR:':{WARNLEN}}Can not convert the orcaversion, needed for input generation!")

        # build up call:
        orcainput = orcainput_start.copy()
        # set functional
        if self.job["func"] in dfa_settings.composite_method_basis.keys() and self.job["basis"] == dfa_settings.composite_method_basis.get(self.job["func"], "NONE"):
            if self.job["func"] == "b3lyp-3c":
                orcainput["functional"] = [f"! {'b3lyp'}"]
                orcainput["basis"] = [f"! {self.job['basis']}"]
                orcainput["gcp"] = [f"! GCP(DFT/SV(P))"]
            else:
                orcainput["functional"] = [
                    f"! {dfa_settings.functionals.get(self.job['func']).get('orca')}"
                ]
        else:
            if self.job["func"] in ("kt2", "kt2-novdw"):
                orcainput["functional"] = [
                    "%method",
                    "  method dft",
                    "  functional gga_xc_kt2",
                    "end",
                ]
            else:
                orcainput["functional"] = [
                    f"! {dfa_settings.functionals.get(self.job['func']).get('orca')}"
                ]
            # set basis set
            orcainput["basis"] = [f"! {self.job['basis']}"]
            # set gcp:
            if "DOGCP" in self.job["prepinfo"] or self.job["func"] == "b3lyp-3c":
                gcp_keywords = {
                    "minis": "MINIS",
                    "sv": "SV",
                    "6-31g(d)": "631GD",
                    "def2-sv(p)": "SV(P)",
                    "def2-svp": "SVP",
                    "def2-tzvp": "TZ",
                }
                if self.job["basis"].lower() in gcp_keywords.keys():
                    if self.job[
                        "func"
                    ] in dfa_settings.composite_method_basis.keys() and self.job[
                        "basis"
                    ] == dfa_settings.composite_method_basis.get(
                        self.job["func"], "NONE"
                    ):
                        pass
                    else:
                        orcainput["gcp"] = [
                            f"! GCP(DFT/{gcp_keywords[self.job['basis'].lower()]})"
                        ]
            # set  RI def2/J,   RIJCOSX def2/J gridx6 NOFINALGRIDX,  RIJK def2/JK
            if self.job["func"] in dfa_settings().dh_dfa():
                if nmrprop:
                    orcainput["frozencore"] = ["!NOFROZENCORE"]
                else:
                    orcainput["frozencore"] = ["! Frozencore"]
                def2cbasis = ("Def2-SVP", "Def2-TZVP", "Def2-TZVPP", "Def2-QZVPP")
                if str(self.job["basis"]).upper() in def2cbasis:
                    # --> decide cosx or RIJK
                    if orca5:
                        orcainput["RI-approx"] = [
                            f"! def2/J {str(self.job['basis'])}/C RIJCOSX"
                        ]
                    else:
                        orcainput["RI-approx"] = [
                            f"! def2/J {str(self.job['basis'])}/C RIJCOSX GRIDX7 NOFINALGRIDX"
                        ]
                    # call.append(f"! RIJK def2/JK {str(self.job['basis'])}/C")
                else:
                    if orca5:
                        orcainput["RI-approx"] = [
                            f"! def2/J def2-TZVPP/C RIJCOSX"
                        ]
                        # call.append(f"! RIJK def2/JK def2-TZVPP/C ")                    
                    else:
                        orcainput["RI-approx"] = [
                            f"! def2/J def2-TZVPP/C RIJCOSX GRIDX7 NOFINALGRIDX"
                        ]
                        # call.append(f"! RIJK def2/JK def2-TZVPP/C ")
                if nmrprop:
                    orcainput["mp2"] = [
                        "%mp2",
                        "    RI true",
                        "    density relaxed",
                        "end",
                    ]
                else:
                    orcainput["mp2"] = ["%mp2", "    RI true", "end"]
            elif (
                self.job["func"] in dfa_settings().hybrid_dfa()
                and self.job["func"] != "pbeh-3c"
            ):
                if orca5:
                    orcainput["RI-approx"] = [f"! def2/J RIJCOSX"]
                else:
                    orcainput["RI-approx"] = [f"! def2/J RIJCOSX GRIDX6 NOFINALGRIDX"]
            elif self.job["func"] in dfa_settings.composite_method_basis.keys():
                pass
            else:  # essentially gga
                orcainput["RI-approx"] = ["! RI def2/J"]
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
        # set scfconv or convergence threshold e.g. loosescf or scfconv6

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

        # optimization ancopt or pure orca
        if self.job["jobtype"] == "xtbopt":
            orcainput["job"] = ["! ENGRAD"]
        elif self.job["jobtype"] == "opt":
            orcainput["job"] = ["! OPT"]
            # add thresholds
            orcainput["optthreshold"] = []
        # nprocs
        if int(self.job["omp"]) >= 1:
            orcainput["parallel"] = [
                "%pal",
                "    nprocs {}".format(self.job["omp"]),
                "end",
            ]
        # solvent model
        # upd_solvent = {
        #     "chcl3": "chloroform",
        #     "h2o": "water",
        #     "ch2cl2":"dichloromethane",
        #     "octanol": "1-octanol",
        #     "hexadecane": "N-HEXADECANE",
        # }
        # solventexch = {
        #     "acetone": "Acetone",
        #     "chcl3": "Chloroform",
        #     "acetonitrile": "Acetonitrile",
        #     "ch2cl2": "CH2Cl2",
        #     "dmso": "DMSO",
        #     "h2o": "Water",
        #     "methanol": "Methanol",
        #     "thf": "THF",
        #     "toluene": "Toluene",
        #     "octanol": "Octanol",
        # }
        # if self.job['solvent'] != 'gas':
        #     if self.job['sm'] in ('smd', 'smd_gsolv'):
        #         self.job['solvent'] = upd_solvent.get(self.job['solvent'], self.job['solvent'])
        #         orcainput['solvation'] = [
        #             '%cpcm',
        #             '    smd    true',
        #             (f'    smdsolvent '
        #             f'"{solventexch.get(self.job["solvent"],self.job["solvent"])}"'),
        #             'end',
        #         ]
        #     elif self.job['sm'] == 'cpcm':
        #         orcainput['solvation'] = [(
        #             f"! CPCM("
        #             f"{solventexch.get(self.job['solvent'],self.job['solvent'])})"
        #             ),
        #         ]
        if self.job["solvent"] != "gas":
            if self.job["sm"] in ("smd", "smd_gsolv"):
                orcainput["solvation"] = [
                    "%cpcm",
                    "    smd    true",
                    (
                        f"    smdsolvent "
                        f'"{censo_solvent_db[self.job["solvent"]]["smd"][1]}"'
                    ),
                    "end",
                ]
            elif self.job["sm"] == "cpcm":
                orcainput["solvation"] = [
                    (f"! CPCM(" f"{censo_solvent_db[self.job['solvent']]['cpcm'][1]})")
                ]
        # unpaired, charge, and coordinates
        if xyzfile:
            orcainput["geom"] = [
                (
                    f"* xyzfile {self.job['charge']} "
                    f"{self.job['unpaired']+1} {str(xyzfile)}"
                )
            ]
        else:
            # xyz geometry
            geom, _, _ = t2x(self.job["workdir"])
            orcainput["geom"] = [f"*xyz {self.job['charge']} {self.job['unpaired']+1}"]
            orcainput["geom"].extend(geom)
            orcainput["geom"].append("*")
        # nmr kt2 disp
        if self.job["func"] == "kt2" and ("nmrJ" or "nmrS" in self.job["prepinfo"]):
            orcainput["disp"] = []
        # couplings
        if nmrprop and "nmrJ" in self.job["prepinfo"]:
            if not any(
                [
                    self.job["h_active"],
                    self.job["c_active"],
                    self.job["f_active"],
                    self.job["si_active"],
                    self.job["p_active"],
                ]
            ):
                self.job["h_active"] = True
                self.job["c_active"] = True
                self.job["f_active"] = True
                self.job["si_active"] = True
                self.job["p_active"] = True
            tmp = []
            tmp.append("%eprnmr")
            if self.job["h_active"]:
                tmp.append(" Nuclei = all H { ssfc }")
            if self.job["c_active"]:
                tmp.append(" Nuclei = all C { ssfc }")
            if self.job["f_active"]:
                tmp.append(" Nuclei = all F { ssfc }")
            if self.job["si_active"]:
                tmp.append(" Nuclei = all Si { ssfc }")
            if self.job["p_active"]:
                tmp.append(" Nuclei = all P { ssfc }")
            tmp.append(" SpinSpinRThresh 8.0")
            tmp.append("end")
            orcainput["couplings"] = tmp
        # shielding
        if nmrprop and "nmrS" in self.job["prepinfo"]:
            if not any(
                [
                    self.job["h_active"],
                    self.job["c_active"],
                    self.job["f_active"],
                    self.job["si_active"],
                    self.job["p_active"],
                ]
            ):
                self.job["h_active"] = True
                self.job["c_active"] = True
                self.job["f_active"] = True
                self.job["si_active"] = True
                self.job["p_active"] = True
            tmp = []
            tmp.append("%eprnmr")
            if self.job["h_active"]:
                tmp.append(" Nuclei = all H { shift }")
            if self.job["c_active"]:
                tmp.append(" Nuclei = all C { shift }")
            if self.job["f_active"]:
                tmp.append(" Nuclei = all F { shift }")
            if self.job["si_active"]:
                tmp.append(" Nuclei = all Si { shift }")
            if self.job["p_active"]:
                tmp.append(" Nuclei = all P { shift }")
            tmp.append(" origin giao")
            tmp.append(" giao_2el giao_2el_same_as_scf")
            tmp.append(" giao_1el giao_1el_analytic")
            tmp.append("end")
            orcainput["shieldings"] = tmp

        # uvvis
        if self.job["calc_uvvis"]:
            orcainput["uvvis"] = [
                "%tddft",
                f"  nroots {self.job['nroots']}",
                "end",
            ]

        error_logical = False
        if not orcainput["functional"]:
            error_logical = True
        elif (
            not orcainput["basis"]
            and self.job["func"] not in dfa_settings.composite_method_basis.keys()
        ):
            error_logical = True
        elif not orcainput["geom"]:
            error_logical = True
        if error_logical:
            print("unusable input!")

        tmp = []
        for value in orcainput.values():
            if value:
                tmp.extend(value)
        if returndict:
            return tmp, orcainput
        else:
            return tmp

    def _sp(self, silent=False, filename="sp.out"):
        """
        ORCA input generation and single-point calculation
        """
        outputpath = os.path.join(self.job["workdir"], filename)
        if not self.job["onlyread"]:
            # write input into file "inp"
            with open(
                os.path.join(self.job["workdir"], "inp"), "w", newline=None
            ) as inp:
                for line in self._prep():
                    inp.write(line + "\n")

            if not silent:
                print(f"Running single-point in {last_folders(self.job['workdir'], 2)}")
            # start SP calculation
            with open(outputpath, "w", newline=None) as outputfile:
                # make external call to orca with "inp" as argument
                call = [os.path.join(external_paths["orcapath"], "orca"), "inp"]
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
        # read output
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as inp:
                stor = inp.readlines()
                for i, line in enumerate(stor):
                    if "FINAL SINGLE POINT ENERGY" in line:
                        self.job["energy"] = float(line.split()[4])
                    # check if scf is converged:
                    if "ORCA TERMINATED NORMALLY" in line:
                        self.job["success"] = True
                    
            if not self.job["success"]:
                self.job["energy"] = 0.0
                self.job["success"] = False # FIXME - redundant?
                print(
                    f"{'ERROR:':{WARNLEN}}scf in {last_folders(self.job['workdir'], 2)} "
                    "not converged!"
                )
        else:
            self.job["energy"] = 0.0
            self.job["success"] = False
            print(f"{'WARNING:':{WARNLEN}}{outputpath} doesn't exist!")
        return

    def _gsolv(self):
        """
        Calculate SMD_gsolv, needs ORCA
        if optimization is not performed with ORCA, only the density 
        functional for optimization is employed, 
        from my understanding smd is parametrized at 298 K, therefore it should only
        be used at this temperature.
        energy --> gas phase
        energy2 --> smd_gsolv gsolv contribution
        """
        energy_gas = None
        energy_solv = None
        print(
            f"Running SMD_gsolv calculation in "
            f"{last_folders(self.job['workdir'], 2)}."
        )
        # calculate gas phase
        keepsolv = self.job["solvent"]
        keepsm = self.job["sm"]
        self.job["solvent"] = "gas"
        self.job["sm"] = "gas-phase"
        self._sp(silent=True, filename="sp_gas.out")

        if self.job["success"] == False:
            self.job["energy"] = 0.0
            self.job["energy2"] = 0.0
            print(
                f"{'ERROR:':{WARNLEN}}in gas phase single-point "
                f"of {last_folders(self.job['workdir'], 2):18}"
            )
            return
        else:
            energy_gas = self.job["energy"]
            self.job["energy"] = 0.0
        # calculate in solution
        self.job["solvent"] = keepsolv
        self.job["sm"] = keepsm
        self._sp(silent=True, filename="sp_solv.out")
        if self.job["success"] == False:
            self.job["energy"] = 0.0
            self.job["energy2"] = 0.0
            print(
                f"{'ERROR:':{WARNLEN}}in gas solution phase single-point "
                f"of {last_folders(self.job['workdir'], 2):18}"
            )
            return
        else:
            energy_solv = self.job["energy"]
            self.job["energy"] = 0.0
        if self.job["success"]:
            if energy_solv is None or energy_gas is None:
                self.job["energy"] = 0.0
                self.job["energy"] = 0.0
                self.job["success"] = False
                print(
                    f"{'ERROR:':{WARNLEN}}in SMD_Gsolv calculation "
                    f"{last_folders(self.job['workdir'], 2):18}"
                )
            else:
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
        return

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