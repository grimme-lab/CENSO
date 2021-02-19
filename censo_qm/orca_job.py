"""
Contains OrcaJob class for calculating ORCA related properties of conformers.
"""
from collections import OrderedDict
import os
import time
import subprocess
import shutil
from .cfg import (
    CODING,
    ENVIRON,
    WARNLEN,
    censo_solvent_db,
    external_paths,
    composite_method_basis,
    composite_dfa,
    gga_dfa,
    hybrid_dfa,
    dh_dfa,
    disp_already_included_in_func,
)
from .utilities import last_folders, t2x, x2t, print
from .qm_job import QmJob


class OrcaJob(QmJob):
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

    def __init__(self, rank, *args, **kwargs):
        QmJob.__init__(self, rank, *args, **kwargs)

    def _prep_input(self, xyzfile=False, returndict=False):
        """
        cefine preparation step analogue

        use:
        xyzfile --> if set : * xyzfile ...  xyzfile.xyz 
        returns
        call --> list with settings
        """

        # understood commands in prepinfo:
        # DOGCP = uses gcp if basis known for gcp
        # 

        orcainput_start = OrderedDict(
            [
                ("functional", None),
                ("disp", None),
                ("basis", None),
                ("gcp", None),
                ("RI-approx", None),
                ("grid", None),
                ("scfconv", None),
                ("frozencore", None),
                ("mp2", None),
                ("default", [
                            "! smallprint printgap noloewdin",
                            "! NOSOSCF",
                            "%MaxCore 8000",
                            "%output",
                            "       print[P_BondOrder_M] 1",
                            "       print[P_Mayer] 1",
                            "       print[P_basis] 2",
                            "end",
                        ]
                ),
                ("job", None),
                ("optthreshold", None),
                ("parallel", None),
                ("solvation", None),
                ("geom", None),
                ("couplings", None),
                ("shieldings", None),
            ]
        )
        if "nmrJ" or "nmrS" in self.job["prepinfo"]:
            nmrprop = True
        else:
            nmrprop = False

        # definitions taken from .cfg:
        # composite_dfa
        # gga_dfa
        # hybrid_dfa
        # dh_dfa
        # disp_already_included_in_func


        # build up call:
        orcainput = orcainput_start.copy()
        # set functional
        if (self.job["func"] in composite_dfa and 
            self.job['basis'] == composite_method_basis.get(self.job["func"],'NONE')):
            orcainput["functional"] = [f"! {self.job['func']}"]
        else:
            if self.job["func"] == "kt2":
                orcainput["functional"] = [
                    "%method",
                    "  method dft",
                    "  functional gga_xc_kt2",
                    "end",
                ]
            elif self.job["func"] == "dsd-blyp":
                orcainput["functional"] = [f"! ri-{self.job['func']}"]
            else:
                orcainput["functional"] = [f"! {self.job['func']}"]
            # set basis set
            orcainput["basis"] = [f"! {self.job['basis']}"]
            # set gcp:
            if "DOGCP" in self.job['prepinfo']:
                gcp_keywords  = {
                        'minis': "MINIS",
                        "sv": "SV",
                        "6-31g(d)": "631GD",
                        'def2-sv(p)': "SV(P)",
                        'def2-svp': "SVP",
                        'def2-tzvp': "TZ",
                }
                if self.job['basis'].lower() in gcp_keywords.keys():
                    if self.job['func'] in composite_dfa and self.job['basis'] == composite_method_basis.get(self.job["func"],'NONE'):
                        pass
                    else:
                        orcainput["gcp"] = [f"! GCP(DFT/{gcp_keywords[self.job['basis'].lower()]})"]
            # set  RI def2/J,   RIJCOSX def2/J gridx6 NOFINALGRIDX,  RIJK def2/JK
            if self.job["func"] in dh_dfa:
                if nmrprop:
                    orcainput["frozencore"] = ["!NOFROZENCORE"]
                else:
                    orcainput["frozencore"] = ["! Frozencore"]
                def2cbasis = ("Def2-SVP", "Def2-TZVP", "Def2-TZVPP", "Def2-QZVPP")
                if str(self.job["basis"]).upper() in def2cbasis:
                    # --> decide cosx or RIJK
                    orcainput["RI-approx"] = [
                        f"! def2/J {str(self.job['basis'])}/C RIJCOSX GRIDX7 NOFINALGRIDX"
                    ]
                    # call.append(f"! RIJK def2/JK {str(self.job['basis'])}/C")
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
            elif self.job["func"] in hybrid_dfa:
                orcainput["RI-approx"] = [f"! def2/J RIJCOSX GRIDX6 NOFINALGRIDX"]
            elif self.job["func"] in composite_dfa:
                pass
            else:  # essentially gga
                orcainput["RI-approx"] = ["! RI def2/J"]
        # set grid
        if self.job["func"] in dh_dfa or self.job["func"] in hybrid_dfa:
            orcainput["grid"] = ["! grid5 nofinalgrid"]
        else:
            orcainput["grid"] = ["! grid4 nofinalgrid"]

        orcainput["scfconv"] = ["! scfconv6"]
        # set scfconv or convergence threshold e.g. loosescf or scfconv6

        extension = {
            "low": {"grid": ["! grid4 nofinalgrid"], "scfconv": ["! loosescf"]},
            "low+": {"grid": ["! grid4 nofinalgrid"], "scfconv": ["! scfconv6"]},
            "high": {"grid": ["! grid4 nofinalgrid"], "scfconv": ["! scfconv7"]},
            "high+": {"grid": ["! grid5 nofinalgrid"], "scfconv": ["! scfconv7"]},
        }
        if self.job["prepinfo"]:
            if isinstance(self.job["prepinfo"], list):
                if self.job["prepinfo"][0] in extension.keys():
                    orcainput["grid"] = extension[self.job["prepinfo"][0]]["grid"]
                    orcainput["scfconv"] = extension[self.job["prepinfo"][0]]["scfconv"]
            else:
                pass

        # add dispersion
        if self.job["func"] not in disp_already_included_in_func:
            orcainput["disp"] = ["! d3bj"]
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
            geom, _ = t2x(self.job["workdir"])
            orcainput["geom"] = [f"*xyz {self.job['charge']} {self.job['unpaired']+1}"]
            orcainput["geom"].extend(geom)
            orcainput["geom"].append("*")
        #nmr kt2 disp
        if self.job["func"] == "kt2" and ("nmrJ" or "nmrS" in self.job["prepinfo"]):
            orcainput["disp"] = []
        # couplings
        if nmrprop and "nmrJ" in self.job["prepinfo"]:
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

        error_logical = False
        if not orcainput["functional"]:
            error_logical = True
        elif not orcainput["basis"] and self.job["func"] not in composite_dfa:
            error_logical = True
        elif not orcainput["geom"]:
            error_logical = True
        if error_logical:
            print("unusable input!")

        tmp = []
        for key, value in orcainput.items():
            if value:
                tmp.extend(value)
        if returndict:
            return tmp, orcainput
        else:
            return tmp

    def _sp(self, silent=False):
        """
        ORCA input generation and single-point calculation
        """
        if not self.job["onlyread"]:
            with open(
                os.path.join(self.job["workdir"], "inp"), "w", newline=None
            ) as inp:
                for line in self._prep_input():
                    inp.write(line + "\n")

            # Done writing input!
            time.sleep(0.02)
            if not silent:
                print(f"Running single-point in {last_folders(self.job['workdir'], 2)}")
            # start SP calculation
            with open(
                os.path.join(self.job["workdir"], "sp.out"), "w", newline=None
            ) as outputfile:
                call = [os.path.join(external_paths["orcapath"], "orca"), "inp"]
                subprocess.call(
                    call,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                )
            time.sleep(0.05)
        # check if scf is converged:
        if os.path.isfile(os.path.join(self.job["workdir"], "sp.out")):
            with open(
                os.path.join(self.job["workdir"], "sp.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                for line in stor:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        self.job["energy"] = float(line.split()[4])
                    if "ORCA TERMINATED NORMALLY" in line:
                        self.job["success"] = True
            if not self.job["success"]:
                self.job["energy"] = 0.0
                self.job["success"] = False
                print(
                    f"{'ERROR:':{WARNLEN}}scf in {last_folders(self.job['workdir'], 2)} "
                    "not converged!"
                )
        else:
            self.job["energy"] = 0.0
            self.job["success"] = False
            print(
                f"{'WARNING:':{WARNLEN}}{os.path.join(self.job['workdir'], 'sp.out')} "
                "doesn't exist!"
            )
        return

    def _smd_gsolv(self):
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
        self._sp(silent=True)

        if self.job["success"] == False:
            self.job["success"] = False
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
        # mv inp inp_solv sp.out sp_solv.out
        try:
            shutil.move(
                os.path.join(self.job["workdir"], "inp"),
                os.path.join(self.job["workdir"], "inp_gas"),
            )
            shutil.move(
                os.path.join(self.job["workdir"], "sp.out"),
                os.path.join(self.job["workdir"], "sp_gas.out"),
            )
        except FileNotFoundError:
            pass
        # calculate in solution
        self.job["solvent"] = keepsolv
        self.job["sm"] = keepsm
        self._sp(silent=True)
        if self.job["success"] == False:
            self.job["success"] = False
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
        # mv inp inp_solv sp.out sp_solv.out
        try:
            shutil.move(
                os.path.join(self.job["workdir"], "inp"),
                os.path.join(self.job["workdir"], "inp_solv"),
            )
            shutil.move(
                os.path.join(self.job["workdir"], "sp.out"),
                os.path.join(self.job["workdir"], "sp_solv.out"),
            )
        except FileNotFoundError:
            pass
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
                self.job["erange1"] = {self.job["temperature"] :energy_solv - energy_gas}
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
        if not self.job["onlyread"]:
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
                newcoord.write(f"   orca bin= {os.path.join(self.job['progpath'], 'orca')}")
                newcoord.write("$end")

            with open(
                os.path.join(self.job["workdir"], "inp"), "w", newline=None
            ) as inp:
                for line in self._prep_input(xyzfile="inp.xyz"):
                    inp.write(line + "\n")
            # Done writing input!
            callargs = [
                self.job["xtb_driver_path"],
                "coord",
                "--opt",
                self.job["optlevel"],
                "--orca",
                "-I",
                "opt.inp"
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
                out.write(f"   orca bin= {os.path.join(self.job['progpath'], 'orca')}")
                out.write("$end \n")
            time.sleep(0.02)
            with open(
                os.path.join(self.job["workdir"], output), "w", newline=None
            ) as outputfile:
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
            # check if optimization finished correctly:
        if os.path.isfile(os.path.join(self.job["workdir"], output)):
            with open(
                os.path.join(self.job["workdir"], output),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
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
            with open(
                os.path.join(self.job["workdir"], output),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                for line in inp:
                    if "av. E: " in line:
                        # self.job["ecyc"].append(float(line.split("Eh")[0].split()[-1]))
                        self.job["ecyc"].append(float(line.split("->")[-1]))
                    if " :: gradient norm      " in line:
                        self.job["grad_norm"] = float(line.split()[3])
        else:
            print(
                f"{'WARNING:':{WARNLEN}}"
                f"{os.path.join(self.job['workdir'], output)} doesn't exist!"
            )
            error_logical = True
        if not error_logical:
            try:
                self.job["energy"] = self.job["ecyc"][-1]
                self.job["success"] = True
            except:
                error_logical = True
        if error_logical:
            self.job["energy"] = 0.0
            self.job["success"] = False
            self.job["converged"] = False
            self.job["ecyc"] = []
            self.job["grad_norm"] = 10.0

        # convert optimized xyz to coord file
        x2t(self.job["workdir"], infile="inp.xyz")
        return

    def _nmrS(self):
        """
        ORCA NMR shielding constant calculation
        """
        with open(os.path.join(self.job["workdir"], "inpS"), "w", newline=None) as inp:
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
            call = [os.path.join(self.job["progpath"], "orca"), "inpS"]
            subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.job["workdir"],
                stdout=outputfile,
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
        # generate input   # double hybrids not implemented
        with open(os.path.join(self.job["workdir"], "inpJ"), "w", newline=None) as inp:
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
            call = [os.path.join(self.job["progpath"], "orca"), "inpJ"]
            subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.job["workdir"],
                stdout=outputfile,
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

    def execute(self):
        """
        Choose what to execute for the jobtype
        use:
        prep --> ignore 
        sp --> _sp
        cosmors --> not with orca
        opt --> pure opt with ORCA
        xtbopt --> opt with xtb as driver
        rrhoxtb --> _rrho()
        """
        if self.job["jobtype"] == "prep":
            self.job["success"] = True
            pass
        elif self.job["jobtype"] == "xtb_sp":
            self._xtb_sp()
        elif self.job["jobtype"] in ("sp", "sp_implicit"):
            self._sp()
        elif self.job["jobtype"] == "opt":
            print("RUNNING xtbopt!!!")
            # self._opt()
            self._xtbopt()
        elif self.job["jobtype"] == "xtbopt":
            self._xtbopt()
        elif self.job["jobtype"] == "rrhoxtb":
            self._xtbrrho()
        # elif self.job['jobtype'] == "rrhoorca":
        #     self._rrho()
        elif self.job["jobtype"] == "smd_gsolv":
            self._smd_gsolv()
        elif self.job["jobtype"] == "couplings_sp":
            self._nmrJ()
        elif self.job["jobtype"] == "shieldings_sp":
            self._nmrS()
        elif self.job["jobtype"] == "genericout":
            self._genericoutput()
        elif self.job["jobtype"] in ("gbsa_gsolv", "alpb_gsolv"):
            if self.job["prepinfo"]:
                tmp_solvent = self.job["solvent"]
                self.job["solvent"] = "gas"
                self._sp()
                if not self.job["success"]:
                    return
                self.job["solvent"] = tmp_solvent
            self._xtb_gsolv()
        else:
            print(f"JOBTYPE {self.job['jobtype']} UNKNOWN!")
