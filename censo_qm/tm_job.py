"""
Contains TmJob class for calculating TM related properties of conformers.
"""
import os
import math
import sys

try:
    from math import isclose
except ImportError:
    from .utilities import isclose
import time
import subprocess
import shutil
from .cfg import (
    CODING,
    ENVIRON,
    AU2KCAL,
    WARNLEN,
    censo_solvent_db,
    external_paths,
    cosmors_param,
    dfa_settings,
)
from .utilities import last_folders, print
from .qm_job import QmJob


class TmJob(QmJob):
    """
    Perform calculations with TM
    - create input with cefine
    - single-point calculation
    - COSMO-RS calculation
    - optimization with xTB as driver
    - shielding constant calculations
    - coupling constant calculations
    - writing of generic output for shielding and coupling constants
    """

    def __init__(self, rank, *args, **kwargs):
        QmJob.__init__(self, rank, *args, **kwargs)

    def _prep_cefine(self):
        """
        Run define for Turbomole calculation using comandline program cefine.
        """
        controlappend = []
        if self.job["basis"] == "def2-QZVP(-gf)":
            self.job["basis"] = "def2-QZVP"
            removegf = True
        else:
            removegf = False
        if self.job["basis"] == "def2-TZVPD(-f)":
            self.job["basis"] = "def2-TZVPD"
            remove_f = True
        else:
            remove_f = False

        # build cefine call:
        minimal_call = [
            external_paths["cefinepath"],
            "-func",
            dfa_settings.functionals.get(self.job["func"]).get("tm"),
            "-bas",
            str(self.job["basis"]),
            "-sym",
            "c1",
            "-noopt",
        ]
        if remove_f:
            minimal_call.extend(["-fpol"])
        if removegf:
            # remove -fg functions from def2-QZVP basis set
            minimal_call.extend(["-gf"])

        extension = {
            "clear": [],
            "low": ["-grid", "m3", "-scfconv", "6"],
            "low+": ["-grid", "m4", "-scfconv", "6"],
            "high": ["-grid", "m4", "-scfconv", "7"],
            "high+": ["-grid", "m5", "-scfconv", "7"],
        }

        dogcp = False  # uses gcp with basis and is added to controlappend
        call = minimal_call
        if self.job["prepinfo"]:
            if isinstance(self.job["prepinfo"], list):
                if self.job["prepinfo"][0] in extension.keys():
                    call.extend(extension[self.job["prepinfo"][0]])
                if "DOGCP" in self.job["prepinfo"]:
                    _ = self.job["prepinfo"].pop(self.job["prepinfo"].index("DOGCP"))
                    dogcp = True
                if len(self.job["prepinfo"]) > 1:
                    call.extend(self.job["prepinfo"][1:])
            else:
                call.extend(extension["low"])
        else:
            call.extend(extension["low"])

        # r2scan-3c use at least grid m4
        if self.job["func"] == "r2scan-3c":
            if "m3" in call:
                call[call.index("m3")] = "m4"

        # settings which request no dispersion:
        # kt1 and kt2 ar only used in part4 (NMR) and fail when using unparameterized disp
        # wb97x-v already includes VV10

        # dispersion correction information
        if dfa_settings.functionals.get(self.job["func"]).get("disp") == "composite":
            pass
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "d3bj":
            call.extend(["-d3"])
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "d3(0)":
            call.extend(["-zero"])
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "d4":
            call.extend(["-d4"])
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "nl":
            call.extend(["-novdw"])
            controlappend.append("$donl")
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "novdw":
            call.extend(["-novdw"])
        elif dfa_settings.functionals.get(self.job["func"]).get("disp") == "included":
            call.extend(["-novdw"])
        else:
            print(
                f" {dfa_settings.functionals.get(self.job['func']).get('disp')} unknown dispersion option!"
            )

        if "-novdw" in call:
            requestnovdw = True
        else:
            requestnovdw = False

        # update unpaired electrons
        if int(self.job["unpaired"]) > 0:
            call = call + ["-uhf", str(self.job["unpaired"])]
        # update charge:
        if int(self.job["charge"]) != 0:
            call = call + ["-chrg", str(self.job["charge"])]

        broken = False
        # call cefine:
        for _ in range(2):
            #print(call)
            tmp = subprocess.check_output(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.job["workdir"],
                env=ENVIRON,
            )
            time.sleep(0.08)
            output = tmp.decode("utf-8").splitlines()
            for line in output:
                if "define ended abnormally" in line:
                    self.job["success"] = False
                    broken = True
                elif "define_huge" in line:
                    print(f"{'ERROR:':{WARNLEN}}define_huge: not found!")
                    self.job["success"] = False
                    broken = True
                elif "Could not find the beginning of the MO-eigenvalue data" in line:
                    self.job["success"] = False
                    broken = True
            # check if wrong functional was written by cefine
            if not os.path.isfile(os.path.join(self.job["workdir"], "control")):
                line = (
                    f"{'ERROR:':{WARNLEN}}CONF{self.id}: The control file was not generated"
                    " by cefine! Can not proceed with calculation!"
                )
                self.job["internal_error"].append(line)
                self.job["success"] = False
                broken = True
            if broken:
                path_to_define = shutil.which("define")
                with open(os.path.join(self.job["workdir"], "cefine.out"),
                        "w",
                        encoding=CODING,
                        newline=None,
                    ) as out:
                    for line in output:
                        out.write(line+'\n')
                    out.write("define path: "+ path_to_define+'\n')
                    out.write("cefine path: "+ external_paths["cefinepath"]+'\n')
                return
            # check if wrong functional was written by cefine
            with open(
                os.path.join(self.job["workdir"], "control"),
                "r",
                encoding=CODING,
                newline=None,
            ) as control:
                checkup = control.readlines()
            for line in checkup:
                if "functional" in line:
                    testfunc = [
                        dfa_settings.functionals.get(self.job["func"]).get("tm")
                    ]
                    if dfa_settings.functionals.get(self.job["func"]).get("tm") in (
                        "b973c",
                        "b97-3c",
                    ):
                        testfunc.extend(["b973c", "b97-3c"])
                    elif self.job["func"] == "b3lyp-3c":
                        testfunc.extend(["b3-lyp"])
                    if not any(func in line for func in testfunc):
                        print(
                            "Wrong functional in control file"
                            " in {}".format(last_folders(self.job["workdir"], 2))
                        )
                        self.job["success"] = False
                        self.job["internal_error"].append("prep-failed")
                    else:
                        self.job["success"] = True
                        break
        if not self.job["success"]:
            return

        if dfa_settings.functionals.get(self.job["func"]).get("tm") in ("kt2", "kt1"):
            # update functional to kt2
            with open(
                os.path.join(self.job["workdir"], "control"), "r", newline=None
            ) as inp:
                tmp = inp.readlines()
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as out:
                for line in tmp:
                    if "functional" in line:
                        out.write("   functional xcfun set-gga  \n")
                        if (
                            dfa_settings.functionals.get(self.job["func"]).get("tm")
                            == "kt2"
                        ):
                            out.write("   functional xcfun kt2 1.0  \n")
                        elif (
                            dfa_settings.functionals.get(self.job["func"]).get("tm")
                            == "kt1"
                        ):
                            out.write("   functional xcfun kt1 1.0  \n")
                    else:
                        out.write(line + "\n")
            time.sleep(0.02)
        # modify the control file
        # solvent_dcosmors = {
        #     "acetone": [" epsilon= 20.7", "$dcosmo_rs file=propanone_25.pot"],
        #     "chcl3": [" epsilon= 4.8", "$dcosmo_rs file=chcl3_25.pot"],
        #     "acetonitrile": [" epsilon= 36.6", "$dcosmo_rs file=acetonitrile_25.pot"],
        #     "ch2cl2": [" epsilon= 9.1", "$dcosmo_rs file=chcl3_25.pot"],
        #     "dmso": [" epsilon= 47.2", "$dcosmo_rs file=dimethylsulfoxide_25.pot"],
        #     "h2o": [" epsilon= 80.1", "$dcosmo_rs file=h2o_25.pot"],
        #     "methanol": [" epsilon= 32.7", "$dcosmo_rs file=methanol_25.pot"],
        #     "thf": [" epsilon= 7.6", "$dcosmo_rs file=thf_25.pot"],
        #     "toluene": [" epsilon= 2.4", "$dcosmo_rs file=toluene_25.pot"],
        #     "octanol": [" epsilon= 9.86", "$dcosmo_rs file=octanol_25.pot"],
        #     "woctanol": [" epsilon= 8.1", "$dcosmo_rs file=wet-octanol_25.pot"],
        #     "hexadecane": [" epsilon= 2.08", "$dcosmo_rs file=hexadecane_25.pot"],
        #     "dmf": [" epsilon= 38.3", "$dcosmo_rs file="],  # not in standard TM parameter folder!
        # }
        if self.job["solvent"] not in ("gas", "gas-phase", None):
            solvent_dcosmors = {
                self.job["solvent"]: [
                    f" epsilon= {censo_solvent_db[self.job['solvent']]['DC']}",
                    f"$dcosmo_rs file={censo_solvent_db[self.job['solvent']]['dcosmors'][1]}_25.pot",
                ]
            }
        # handle solvents:

        if self.job["sm"] == "dcosmors" and self.job["solvent"] != "gas":
            try:
                filename = solvent_dcosmors.get(self.job["solvent"])[1].split("=")[1]
            except IndexError:
                filename = ""
            if not os.path.isfile(
                os.path.join(
                    os.path.dirname(
                        os.path.dirname(os.path.dirname(shutil.which("ridft")))
                    ),
                    "parameter/" + filename,
                )
            ):
                if os.path.isfile(
                    os.path.expanduser(os.path.join("~/.censo_assets/", filename))
                ):
                    tmp = solvent_dcosmors.get(self.job["solvent"], "")
                    tmp[4] = "$dcosmo_rs file=" + str(
                        os.path.expanduser(os.path.join("~/.censo_assets/", filename))
                    )
                    solvent_dcosmors[self.job["solvent"]] = tmp
                else:
                    line = (
                        f"{'WARNING:':{WARNLEN}}DCOSMO-RS potential file not found!"
                        " Trying file without verification!"
                    )
                    # print(line)
                    self.job["internal_error"].append(line)
        if self.job["solvent"] != "gas" and self.job["sm"] in ("cosmo", "dcosmors"):
            if solvent_dcosmors.get(self.job["solvent"], "not found!") == "not found!":
                print(
                    f"{'ERROR:':{WARNLEN}} Solvent {self.job['solvent']} is not known for cefine!"
                )
                self.job["success"] = False
                self.job["internal_error"].append("prep-failed")
                return
            else:
                controlappend.append("$cosmo")
                # write epsilon (dielectric constant)
                controlappend.append(solvent_dcosmors.get(self.job["solvent"], "")[0])
                if self.job["jobtype"] not in ("opt-rot", "opt-rot_sp"):
                    controlappend.append(" cavity closed")
                    controlappend.append(" use_contcav")
                    controlappend.append(" nspa=272")
                    controlappend.append(" nsph=162")
                    controlappend.append("$cosmo_isorad")
                # write parameterfile for dcosmors
                if self.job["sm"] == "dcosmors":
                    controlappend.append(
                        solvent_dcosmors.get(self.job["solvent"], "")[1]
                    )
        if self.job["jobtype"] in ("opt-rot", "opt-rot_sp"):
            controlappend.append("$scfinstab dynpol nm")
            for i in self.job["freq_or"]:
                controlappend.append(f" {i}")  # e.g. 589
            controlappend.append("$velocity gauge")
            controlappend.append("$rpaconv 4")

        if dogcp:
            if dfa_settings.functionals.get(self.job["func"]).get(
                "tm"
            ) in dfa_settings.composite_method_basis.keys() and self.job[
                "basis"
            ] == dfa_settings.composite_method_basis.get(
                dfa_settings.functionals.get(self.job["func"]).get("tm"), "NONE"
            ):
                pass
            else:
                if self.job["basis"] == "def2-SV(P)":
                    controlappend.append("$gcp dft/sv(p)")
                else:
                    controlappend.append(
                        f"$gcp dft/{self.job['basis'].lower().replace('-', '')}"
                    )

        # write to control file:
        with open(
            os.path.join(self.job["workdir"], "control"),
            "r",
            encoding=CODING,
            newline=None,
        ) as control:
            tmp = control.readlines()
        # check if dispersion is found
        nodisp = True
        replacewatm = ""
        needatm = ("b97-3c", "b973c")
        for line in tmp:
            if "$disp" in line:
                nodisp = False
                if dfa_settings.functionals.get(self.job["func"])["tm"] in needatm:
                    if "$disp3 -bj -abc" not in line:
                        replacewatm = "$disp3 -bj -abc"
        if nodisp and (self.job["func"] not in needatm) and not requestnovdw:
            controlappend.append("$disp3 -bj")
        elif nodisp and (self.job["func"] in needatm) and not requestnovdw:
            controlappend.append("$disp3 -bj -abc")
        if nodisp and requestnovdw:
            replacewatm = " "
        if controlappend:
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as newcontrol:
                for line in tmp[:-1]:
                    if "$end" in line:
                        pass
                    if "$disp" in line and replacewatm:
                        newcontrol.write(replacewatm + "\n")
                    else:
                        newcontrol.write(line)
                for line in controlappend:
                    newcontrol.write(line + "\n")
                newcontrol.write("$end\n")
        if self.job["copymos"]:
            if self.job["unpaired"] > 0:
                molist = ["alpha", "beta"]
            else:
                molist = ["mos"]
            try:
                for item in molist:
                    tmp_from = os.path.join(
                        "CONF" + str(self.id), self.job["copymos"], item
                    )
                    tmp_to = os.path.join(self.job["workdir"], item)
                    shutil.copy(tmp_from, tmp_to)
            except FileNotFoundError:
                pass
            for item in molist:
                if (
                    not os.path.isfile(os.path.join(self.job["workdir"], item))
                    or os.stat(os.path.join(self.job["workdir"], item)).st_size == 0
                ):
                    print(f"{'ERROR:':{WARNLEN}} {item} is missing!")
        # NMR part
        if self.job["jobtype"] in (
            "couplings",
            "couplings_sp",
            "shieldings",
            "shieldings_sp",
        ):
            with open(
                os.path.join(self.job["workdir"], "control"),
                "r",
                encoding=CODING,
                newline=None,
            ) as control:
                tmp = control.readlines()
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as newcontrol:
                for line in tmp:
                    if "rpacor" in line:
                        rpacor = 10000
                        try:
                            tmpval = float(line.strip().split()[-1])
                            if tmpval > rpacor:
                                rpacor = tmpval
                        except (ValueError, IndexError):
                            pass
                        tmp[tmp.index(line)] = f"$rpacor {str(rpacor)} \n"
                    elif "$cosmo_isorad" in line:
                        tmp.pop(tmp.index(line))
                for line in tmp[:-1]:
                    if "$end" in line:
                        pass
                    else:
                        newcontrol.write(line)
                newcontrol.write("$ncoupling\n")
                newcontrol.write(" simple\n")
                # fc sd pso dso nofcsdcross
                newcontrol.write(" thr=0.0\n")
                nucsel1 = "$nucsel "
                nucsel2 = "$nucsel2 "
                if self.job["h_active"]:
                    nucsel1 = nucsel1 + '"h" '
                    nucsel2 = nucsel2 + '"h" '
                if self.job["c_active"]:
                    nucsel1 = nucsel1 + '"c" '
                    nucsel2 = nucsel2 + '"c" '
                if self.job["f_active"]:
                    nucsel1 = nucsel1 + '"f" '
                    nucsel2 = nucsel2 + '"f" '
                if self.job["si_active"]:
                    nucsel1 = nucsel1 + '"si" '
                    nucsel2 = nucsel2 + '"si" '
                if self.job["p_active"]:
                    nucsel1 = nucsel1 + '"p" '
                    nucsel2 = nucsel2 + '"p" '
                if any(
                    [
                        self.job["h_active"],
                        self.job["c_active"],
                        self.job["f_active"],
                        self.job["si_active"],
                        self.job["p_active"],
                    ]
                ):
                    newcontrol.write(nucsel1 + "\n")
                    newcontrol.write(nucsel2 + "\n")
                elif all(
                    [
                        self.job["h_active"],
                        self.job["c_active"],
                        self.job["f_active"],
                        self.job["si_active"],
                        self.job["p_active"],
                    ]
                ):
                    # don't write nucsel, every shielding, coupling will be calculated
                    pass
                newcontrol.write("$rpaconv 8\n")
                newcontrol.write("$end")
        time.sleep(0.15)

    # ****************************end cefine************************************

    def _sp(self, silent=False, outfile="ridft.out"):
        """
        Turbomole single-point calculation, needs previous cefine run
        """
        fullpath = os.path.join(self.job["workdir"], outfile)
        if not self.job["onlyread"]:
            # perform calculation:
            if not silent:
                print(
                    f"Running single-point in {last_folders(self.job['workdir'], 2):18}"
                )
            with open(fullpath, "w", newline=None) as outputfile:
                subprocess.call(
                    ["ridft"],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
        # process output:
        time.sleep(0.02)
        # check if scf is converged:
        if os.path.isfile(fullpath):
            with open(fullpath, "r", encoding=CODING, newline=None) as inp:
                stor = inp.readlines()
                if " ENERGY CONVERGED !\n" not in stor:
                    print(
                        f"{'ERROR:':{WARNLEN}}scf in "
                        f"{last_folders(self.job['workdir'], 2):18} not converged!"
                    )
                    self.job["success"] = False
                    self.job["energy"] = 0.0
                    self.job["energy2"] = 0.0
                    return
                # get energy:
                for line in stor:
                    if "|  total energy      = " in line:
                        try:
                            self.job["energy"] = float(line.strip().split()[4])
                            self.job["success"] = True
                            if self.job["jobtype"] == "sp_implicit":
                                self.job["energy2"] = 0.0
                        except (ValueError, IndexError):
                            print(
                                f"{'ERROR:':{WARNLEN}}while converting energy "
                                f"in: {last_folders(self.job['workdir'], 2):18}"
                            )
                            self.job["success"] = False
                            self.job["energy"] = 0.0
                            self.job["energy2"] = 0.0
                            return
                if not self.job["success"]:
                    self.job["energy"] = 0.0
                    self.job["energy2"] = 0.0
                    return
        else:
            print(f"{'WARNING:':{WARNLEN}}{fullpath} doesn't exist!")
            self.job["success"] = False
            self.job["energy"] = 0.0
            self.job["energy2"] = 0.0
            return

    # ****************************end _sp***************************************

    def _cosmors(self):
        """
        Run COSMO-RS from within censo.
        calculates directly in the workdir. folder COSMO has to be created 
        beforehand.
        energy --> gas phase scf energy
        energy2 --> gsolv contribution

        cosmorsparam = fine or normal
        ctd-param  = which parametrization version e.g. BP_TZVP_19.ctd
        """
        if not self.job["onlyread"]:
            print(
                f"Running COSMO-RS calculation in "
                f"{last_folders(self.job['workdir'], 3):18}"
            )
            # parametrisation
            # ctd and database have to match fine -> fine  or non-fine -> non-fine
            if self.job["ctd-param"] == "automatic":
                if self.job["cosmorsparam"] == "fine":
                    if (
                        "FINE" not in self.job["cosmorssetup"]
                        and "1201" not in self.job["cosmorssetup"]
                        and "1301" not in self.job["cosmorssetup"]
                    ):
                        self.job["cosmorssetup"] = self.job["cosmorssetup"].replace(
                            "BP_TZVP", "BP_TZVPD_FINE"
                        )
                    elif (
                        "FINE" not in self.job["cosmorssetup"]
                        and "1201" in self.job["cosmorssetup"]
                    ):
                        self.job["cosmorssetup"] = self.job["cosmorssetup"].replace(
                            "BP_TZVP", "BP_TZVPD_FINE_HB2012"
                        )
                    elif (
                        "FINE" not in self.job["cosmorssetup"]
                        and "1301" in self.job["cosmorssetup"]
                    ):
                        self.job["cosmorssetup"] = self.job["cosmorssetup"].replace(
                            "BP_TZVP", "BP_TZVPD_FINE_HB2012"
                        )
                elif self.job["cosmorsparam"] == "normal":
                    if "FINE" in self.job["cosmorssetup"]:
                        ## normal cosmors
                        tmp = self.job["cosmorssetup"]
                        tmp = tmp.replace("_FINE", "")
                        self.job["cosmorssetup"] = tmp.replace("BP_TZVPD", "BP_TZVP")
            else:
                if self.job["cosmorsparam"] == "fine":
                    find = "normal"
                    replacewith = "fine"
                else:
                    find = "fine"
                    replacewith = "normal"
                self.job["ctd-param"] = self.job["ctd-param"].replace(find, replacewith)
                tmp_dat = self.job["cosmorssetup"].split()
                tmp_new = f"ctd = {cosmors_param[self.job['ctd-param']]} "
                for item in tmp_dat:
                    if ".ctd" in item:
                        tmp_new = tmp_new + " ".join(tmp_dat[tmp_dat.index(item) + 1 :])
                self.job["cosmorssetup"] = tmp_new
            # run two single-points:
            # copy mos first
            if self.job["copymos"]:
                if self.job["unpaired"] > 0:
                    molist = ["alpha", "beta"]
                else:
                    molist = ["mos"]
                try:
                    for item in molist:
                        tmp_from = os.path.join(
                            "CONF" + str(self.id), self.job["copymos"], item
                        )
                        tmp_to = os.path.join(self.job["workdir"], item)
                        shutil.copy(tmp_from, tmp_to)
                except FileNotFoundError:
                    pass
                for item in molist:
                    if (
                        not os.path.isfile(os.path.join(self.job["workdir"], item))
                        or os.stat(os.path.join(self.job["workdir"], item)).st_size == 0
                    ):
                        print(f"{'ERROR:':{WARNLEN}}{item} is missing!")
            # first single-point in gas phase!
            tmp_solvent = self.job["solvent"]
            tmp_sm = self.job["solvent"]
            self.job["solvent"] = "gas"
            self.job["sm"] = "gas-phase"
            self._prep_cefine()
            if not self.job["success"]:
                return
            # running single-point in gas phase
            self._sp(silent=True)
            self.job["solvent"] = tmp_solvent
            self.job["sm"] = tmp_sm
            if not self.job["success"]:
                print(
                    f"{'ERROR:':{WARNLEN}}gas-phase single-point calculation failed in: "
                    f"{last_folders(self.job['workdir'], 3):18}"
                )
                return
            # try:
            #     tmp_from = os.path.join(self.job["workdir"], "ridft.out")
            #     tmp_to = os.path.join(self.job["workdir"], "gas_"+ self.job["partx"] +".out")
            #     shutil.copy(tmp_from, tmp_to)
            # except (FileNotFoundError, KeyError):
            #     pass
            with open(
                os.path.join(self.job["workdir"], "out.energy"), "w", newline=None
            ) as out:
                out.write(str(self.job["energy"]) + "\n")
            gas_phase_energy = self.job["energy"]
            self.job["energy"] = 0.0
            # running single-point in ideal conductor!
            with open(
                os.path.join(self.job["workdir"], "control"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                tmp = inp.readlines()
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as out:
                for line in tmp[:-1]:
                    out.write(line + "\n")
                if self.job["cosmorsparam"] == "normal":
                    # normal
                    out.write("$cosmo \n")
                    out.write(" epsilon=infinity \n")
                    out.write(" use_contcav \n")
                    out.write(" cavity closed \n")
                    out.write(" nspa=272 \n")
                    out.write(" nsph=162 \n")
                    out.write("$cosmo_out file=out.cosmo \n")
                    out.write("$end \n")
                else:
                    # fine
                    out.write("$cosmo \n")
                    out.write(" epsilon=infinity \n")
                    out.write(" use_contcav \n")
                    out.write(" cavity closed \n")
                    out.write(" nspa=272 \n")
                    out.write(" nsph=162 \n")
                    out.write("$cosmo_out file=out.cosmo \n")
                    # out.write("$cosmo_isorad \n")
                    # out.write("$cosmo_isodens \n")
                    out.write("$end \n")
            self.job["success"] = False
            self._sp(silent=True)
            if not self.job["success"]:
                print(
                    f"{'ERROR:':{WARNLEN}}single-point in ideal conductor calculation failed in: "
                    f"{last_folders(self.job['workdir'], 3):18}"
                )
                return
            # info from .ensorc # replacement for cosmothermrc
            # fdir=/software/cluster/COSMOthermX16/COSMOtherm/DATABASE-COSMO/BP-TZVP-COSMO autoc
            # cosmors_solv = {
            #     "acetone": ["f = propanone.cosmo "],
            #     "h2o": ["f = h2o.cosmo "],
            #     "chcl3": ["f = chcl3.cosmo "],
            #     "ch2cl2": ["f = ch2cl2.cosmo "],
            #     "acetonitrile": ["f = acetonitrile_c.cosmo "],
            #     "dmso": ["f = dimethylsulfoxide.cosmo "],
            #     "methanol": ["f = methanol.cosmo "],
            #     "thf": ["f = thf.cosmo "],
            #     "toluene": ["f = toluene_c0.cosmo "],
            #     "octanol": ["f = 1-octanol "],
            #     "hexadecane": ["f = n-hexadecane "],
            #     "woctanol": ["f = h2o.cosmo ", "f = 1-octanol "],
            # }

            if self.job["solvent"] not in ("gas", "gas-phase", None):
                if self.job["solvent"] == "woctanol":
                    cosmors_solv = {
                        "woctanol": ["f = h2o.cosmo ", "f = 1-octanol.cosmo "]
                    }
                elif self.job['vapor_pressure']:
                    cosmors_solv = {f"{self.job['solvent']}": [f"f = {'out.cosmo'} "]}
                else:
                    tmp_1 = os.path.splitext(
                        censo_solvent_db[self.job["solvent"]]["cosmors"][1]
                    )[0]
                    filename = f"{tmp_1}.cosmo"
                    cosmors_solv = {f"{self.job['solvent']}": [f"f = {filename} "]}

            mixture = {"woctanol": ["0.27 0.73"]}

            if external_paths.get("dbpath", None) is not None:
                if self.job["cosmorsparam"] == "fine":
                    solv_data = external_paths.get("dbpath_fine")

                else:
                    solv_data = external_paths.get("dbpath_normal")
            else:
                # old behaviour
                if self.job["cosmorsparam"] == "fine":
                    solv_data = os.path.join(
                        os.path.split(self.job["cosmorssetup"].split()[5].strip('"'))[
                            0
                        ],
                        "DATABASE-COSMO/BP-TZVPD-FINE",
                    )
                else:
                    solv_data = os.path.join(
                        os.path.split(self.job["cosmorssetup"].split()[5].strip('"'))[
                            0
                        ],
                        "DATABASE-COSMO/BP-TZVP-COSMO",
                    )
            # test = ['ctd = BP_TZVP_C30_1601.ctd cdir = "/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES"']
            with open(
                os.path.join(self.job["workdir"], "cosmotherm.inp"), "w", newline=None
            ) as out:
                out.write(self.job["cosmorssetup"] + "\n")
                # write from ensorc
                out.write("EFILE VPFILE \n")
                out.write("!!\n")  # for jobname in cosmors
                if len(cosmors_solv[self.job["solvent"]]) > 1:
                    mix = mixture[self.job["solvent"]][0]
                    for line in cosmors_solv[self.job["solvent"]]:
                        out.write(line + "fdir=" + solv_data + " autoc \n")
                elif len(cosmors_solv[self.job["solvent"]]) == 1:
                    mix = "1.0 0.0"
                    if self.job['vapor_pressure']:
                        out.write(
                            cosmors_solv[self.job["solvent"]][0]
                            + " \n")
                    else:
                        out.write(
                            cosmors_solv[self.job["solvent"]][0]
                            + " fdir="
                            + solv_data
                            + " autoc \n"
                        )
                out.write("f = out.cosmo \n")

                if self.job["trange"]:
                    tmp1 = self.job["trange"]
                    tinside = False
                    for temp in tmp1:
                        if isclose(self.job["temperature"], temp, abs_tol=0.6):
                            tinside = True
                    if not tinside:
                        tmp1.append(self.job["temperature"])
                else:
                    tmp1 = [self.job["temperature"]]
                tlist = [str("{:.2f}".format(i - 273.15)) for i in tmp1]
                # henry = "henry  xh={ "+mix+" }  tc=25.0 Gsolv"
                for i in tlist:
                    henry = "henry  xh={ " + mix + " }  tc=" + i + " Gsolv"
                    out.write(henry + "\n")
            time.sleep(0.01)
            # running COSMOtherm
            with open(
                os.path.join(self.job["workdir"], "cosmotherm.out"), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    ["cosmotherm", "cosmotherm.inp"],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
            time.sleep(0.1)
        # calculations done, get output
        if self.job["onlyread"]:
            # if self.job.get("partx", None) is None:
            # get gas phase energy:
            gaspath = os.path.join(self.job["workdir"], "out.energy")
            if os.path.exists(gaspath) and os.path.getsize(gaspath) > 0:
                with open(gaspath, "r", newline=None) as inp:
                    try:
                        gas_phase_energy = float(inp.readline().strip())
                    except Exception:
                        self.job["success"] = False
                        return
            else:
                self.job["success"] = False
                return

        # get T and Gsolv for version > cosmothermX16
        ## volumework:
        R = 1.987203585e-03  # kcal/(mol*K)
        videal = (
            24.789561955 / 298.15
        )  # molar volume for ideal gas at 298.15 K 100.0 kPa
        gsolvt = {}

        cosmothermtab = os.path.join(self.job["workdir"], "cosmotherm.tab")
        try:
            with open(cosmothermtab, "r", encoding=CODING, newline=None) as inp:
                stor = inp.readlines()
            for line in stor:
                if "T=" in line:
                    T = float(line.split()[5])
                    vwork = R * T * math.log(videal * T)
                elif " out " in line:
                    gsolvt[T] = float(line.split()[-1]) / AU2KCAL + vwork / AU2KCAL
            self.job["erange1"] = gsolvt
        except (FileNotFoundError, ValueError):
            print(
                f"{'ERROR:':{WARNLEN}}cosmotherm.tab was not written, this error can be "
                "due to a missing licensefile information, or wrong path "
                "to the COSMO-RS Database."
            )
            self.job["energy"] = 0.0
            self.job["energy2"] = 0.0
            self.job["erange1"] = {}
            self.job["success"] = False
            return
        except IndexError:
            print(f"{'ERROR:':{WARNLEN}}IndexERROR in cosmotherm.tab!")
            self.job["energy"] = 0.0
            self.job["energy2"] = 0.0
            self.job["erange1"] = {}
            self.job["success"] = False
            return
        # cosmothermrd
        gsolv_out = 0.0
        for temp in gsolvt.keys():
            if isclose(self.job["temperature"], temp, abs_tol=0.6):
                gsolv_out = gsolvt[temp]
        temp = float(self.job["temperature"])
        ## volumework:
        R = 1.987203585e-03  # kcal/(mol*K)
        videal = (
            24.789561955 / 298.15
        )  # molar volume for ideal gas at 298.15 K 100.0 kPa
        volwork = R * temp * math.log(videal * temp)
        if not self.job["onlyread"]:
            with open(
                os.path.join(os.path.dirname(self.job["workdir"]), "cosmors.out"),
                "w",
                newline=None,
            ) as out:
                out.write(
                    "This is cosmothermrd (python version in ENSO) (SG,FB,SAW, 06/18)\n"
                )
                out.write("final thermochemical solvation properties in kcal/mol\n")
                out.write(
                    "----------------------------------------------------------\n"
                )
                out.write(
                    " Gsolv({} K)= {:10.3f}\n".format(
                        temp, gsolv_out * AU2KCAL - volwork
                    )
                )
                out.write(" VWork({} K)= {:10.3f}\n".format(temp, volwork))
                out.write(
                    " Gsolv+VWork({} K)= {:10.3f}\n".format(
                        temp, (gsolv_out * AU2KCAL)  # volwork already included!
                    )
                )
            time.sleep(0.01)
        self.job["energy"] = gas_phase_energy
        self.job["energy2"] = gsolv_out  # VOLWORK INCLUDED
        self.job["erange1"][self.job["temperature"]] = gsolv_out  # VOLWORK INCLUDED
        self.job["success"] = True

    # ********************************end _cosmors***********************************
    def _xtbopt(self):
        """
        Turbomole optimization using the ANCOPT optimizer implemented in xTB
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

            callargs = [
                external_paths["xtbpath"],
                "coord",
                "--opt",
                self.job["optlevel"],
                "--tm",
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
                out.write("$end \n")
                callargs.append("-I")
                callargs.append("opt.inp")
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
                    f"{'ERROR:':{WARNLEN}}optimization in "
                    f"{last_folders(self.job['workdir'], 2):18} not converged"
                )
            time.sleep(0.02)
        # read output:
        # check if converged:
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
                            "ERROR: optimization in {:18} not converged".format(
                                last_folders(self.job["workdir"], 2)
                            )
                        )
                        error_logical = True
                        break
                    elif " FAILED TO CONVERGE GEOMETRY " in line:
                        self.job["cycles"] += int(line.split()[7])
                        self.job["converged"] = False
                    elif "*** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
                        self.job["cycles"] += int(line.split()[5])
                        self.job["converged"] = True
                #####
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

    ##### VERSION BEFORE AVERAGING KEEP FOR NOW
    # def _xtbopt(self):
    #     """
    #     Turbomole optimization using the ANCOPT optimizer implemented in xTB
    #     """
    #     error_logical = False
    #     if self.job["fullopt"]:
    #         output = "opt-part2.out"
    #     else:
    #         output = "opt-part1.out"
    #     if not self.job["onlyread"]:
    #         print(f"Running optimization in {last_folders(self.job['workdir'], 2):18}")
    #         files = [
    #             "xtbrestart",
    #             "xtbtopo.mol",
    #             "xcontrol-inp",
    #             "wbo",
    #             "charges",
    #             "gfnff_topo",
    #         ]
    #         for file in files:
    #             if os.path.isfile(os.path.join(self.job["workdir"], file)):
    #                 os.remove(os.path.join(self.job["workdir"], file))

    #         callargs = [
    #             self.job["xtb_driver_path"],
    #             "coord",
    #             "--opt",
    #             self.job["optlevel"],
    #             "--tm",
    #         ]
    #         with open(
    #                 os.path.join(self.job["workdir"], "opt.inp"), "w", newline=None
    #             ) as out:
    #             out.write("$opt \n")
    #             if self.job["optcycles"] is not None and float(self.job["optcycles"]) > 0:
    #                 out.write(f"maxcycle={str(self.job['optcycles'])} \n")
    #                 out.write(f"microcycle={str(self.job['optcycles'])} \n")
    #             out.write("average conv=true \n")
    #             out.write(f"hlow={self.job.get('hlow', 0.01)} \n")
    #             out.write("s6=30.00 \n")
    #             # remove unnecessary sp/gradient call in xTB
    #             out.write("engine=lbfgs\n")
    #             out.write("$end \n")
    #             callargs.append("-I")
    #             callargs.append("opt.inp")
    #         time.sleep(0.02)
    #         with open(
    #             os.path.join(self.job["workdir"], output), "w", newline=None
    #         ) as outputfile:
    #             returncode = subprocess.call(
    #                 callargs,
    #                 shell=False,
    #                 stdin=None,
    #                 stderr=subprocess.STDOUT,
    #                 universal_newlines=False,
    #                 cwd=self.job["workdir"],
    #                 stdout=outputfile,
    #                 env=ENVIRON,
    #             )
    #         if returncode != 0:
    #             error_logical = True
    #             print(
    #                 "ERROR: optimization in {:18} not converged".format(
    #                     last_folders(self.job["workdir"], 2)
    #                 )
    #             )
    #         time.sleep(0.02)
    #     # check if converged:
    #     if os.path.isfile(os.path.join(self.job["workdir"], output)):
    #         with open(
    #             os.path.join(self.job["workdir"], output),
    #             "r",
    #             encoding=CODING,
    #             newline=None,
    #         ) as inp:
    #             stor = inp.readlines()
    #             for line in stor:
    #                 if (
    #                     "external code error" in line
    #                     or "|grad| > 500, something is totally wrong!" in line
    #                     or "abnormal termination of xtb" in line
    #                 ):
    #                     print(
    #                         "ERROR: optimization in {:18} not converged".format(
    #                             last_folders(self.job["workdir"], 2)
    #                         )
    #                     )
    #                     error_logical = True
    #                     break
    #                 elif " FAILED TO CONVERGE GEOMETRY " in line:
    #                     self.job["cycles"] += int(line.split()[7])
    #                     # self.job['cycles'] =  int(line.split()[7]) + self.job['cycles']
    #                     self.job["converged"] = False
    #                 elif "*** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
    #                     self.job["cycles"] += int(line.split()[5])
    #                     # self.job['cycles'] = int(line.split()[5]) + self.job['cycles']
    #                     self.job["converged"] = True
    #         with open(
    #             os.path.join(self.job["workdir"], output),
    #             "r",
    #             encoding=CODING,
    #             newline=None,
    #         ) as inp:
    #             for line in inp:
    #                 if "total energy  :" in line and not "gain" in line:
    #                     self.job["ecyc"].append(float(line.split("Eh")[0].split()[-1]))
    #     else:
    #         print(
    #             "WARNING: {} doesn't exist!".format(
    #                 os.path.join(self.job["workdir"], output)
    #             )
    #         )
    #         error_logical = True
    #     if os.path.isfile(os.path.join(self.job["workdir"], "energy")):
    #         with open(
    #             os.path.join(self.job["workdir"], "energy"),
    #             "r",
    #             encoding=CODING,
    #             newline=None,
    #         ) as energy:
    #             storage = energy.readlines()
    #         try:
    #             self.job["energy"] = float(storage[-2].split()[1])
    #             self.job["success"] = True
    #         except ValueError:
    #             print(
    #                 "ERROR while converting energy in {:18}".format(
    #                     last_folders(self.job["workdir"], 2)
    #                 )
    #             )
    #     else:
    #         error_logical = True
    #     if error_logical:
    #         self.job["energy"] = 0.0
    #         self.job["success"] = False
    #         self.job["converged"] = False
    #         self.job["ecyc"] = []

    # ********************************end _xtbopt***********************************

    def _opt(self):
        """
        Turbomole optimization using JOBEX, with adapted thresholds!
        """
        pass

    def _nmr_coupling(self):
        """
        Turbomole coupling constant calculation.
        """
        print(
            f"Running couplings calculation in {last_folders(self.job['workdir'], 2)}"
        )
        outputpath = os.path.join(self.job["workdir"], "escf.out")
        if not self.job["onlyread"]:
            # escf doesnot allow for mgrids!
            with open(
                os.path.join(self.job["workdir"], "control"), "r", newline=None
            ) as inp:
                tmp = inp.readlines()
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as out:
                for line in tmp:
                    if "gridsize" in line:
                        out.write(f"   gridsize {5}   \n")
                    else:
                        out.write(line + "\n")

            with open(outputpath, "w", newline=None) as outputfile:
                subprocess.call(
                    [external_paths["escfpath"], "-smpcpus", str(self.job["omp"])],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
            time.sleep(0.02)
        # check for convergence
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as inp:
                stor = inp.readlines()
                if "   ****  escf : all done  ****\n" in stor:
                    self.job["success"] = True
                else:
                    print(
                        f"{'ERROR:':{WARNLEN}}coupling calculation failed in "
                        f"{last_folders(self.job['workdir'], 1):18}"
                    )
                    self.job["success"] = False
        else:
            print(f"The file escf.out could not be found!")
            self.job["success"] = False

    def _nmr_shielding(self):
        """
        Turbomole shielding constant calculation.
        """
        print(
            "Running shielding calculation in {}".format(
                last_folders(self.job["workdir"], 2)
            )
        )
        outputpath = os.path.join(self.job["workdir"], "mpshift.out")
        if not self.job["onlyread"]:
            # update grid to m5!
            with open(
                os.path.join(self.job["workdir"], "control"), "r", newline=None
            ) as inp:
                tmp = inp.readlines()
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as out:
                for line in tmp:
                    if "gridsize" in line:
                        out.write(f"   gridsize {'m5'}   \n")
                    if "$disp" in line and self.job["func"] in ("kt2", "kt1"):
                        pass
                    else:
                        out.write(line + "\n")
                time.sleep(0.02)

            with open(outputpath, "w", newline=None) as outputfile:
                subprocess.call(
                    [external_paths["mpshiftpath"], "-smpcpus", str(self.job["omp"])],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
                time.sleep(0.02)
        # check if shift calculation is converged:
        if os.path.isfile(outputpath):
            with open(outputpath, "r", encoding=CODING, newline=None) as inp:
                stor = inp.readlines()
                found = False
                for line in stor:
                    if " ****  mpshift : all done  ****" in line:
                        self.job["success"] = True
                        found = True
                if not found:
                    print(
                        f"{'ERROR:':{WARNLEN}}shielding calculation failed in "
                        f"{last_folders(self.job['workdir'], 2):18}"
                    )
                    self.job["success"] = False
        else:
            print(f"The file mpshift.out could not be found!")
            self.job["success"] = False

    def _genericoutput(self):
        """
        Read shielding and coupling constants and write them to plain output.
        """
        fnameshield = "mpshift.out"
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
                if ">>>>> DFT MAGNETIC SHIELDINGS <<<<<" in line:
                    start = data.index(line)
            for line in data[start:]:
                if "ATOM" in line:
                    splitted = line.split()
                    atom.append(int(splitted[2]))
                    sigma.append(float(splitted[4]))
        except FileNotFoundError:
            print(
                "Missing file: {} in {}, Shielding constants are not written.".format(
                    fnameshield, last_folders(self.job["workdir"], 2)
                )
            )
            self.job["success"] = False
        except ValueError:
            print(
                f"{'ERROR:':{WARNLEN}}ValueError in generic_output, nmrprop.dat can be flawed !"
            )
            self.job["success"] = False
        self.job["success"] = True
        fnamecoupl = "escf.out"
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
                if "Nuclear coupling constants" in line:
                    start = int(data.index(line)) + 3
                if "-----------------------------------" in line:
                    end = int(data.index(line))
            for line in data[start:end]:
                if len(line.split()) in (6, 7):
                    splitted = line.split()
                    atom1.append(int(splitted[1]))
                    atom2.append(int(splitted[4].split(":")[0]))
                    jab.append(float(splitted[5]))
        except FileNotFoundError:
            print(
                "Missing file: {} in {}, Coupling constants are not written.".format(
                    fnamecoupl, last_folders(self.job["workdir"], 2)
                )
            )
            self.job["success"] = False
        except ValueError:
            print(
                f"{'ERROR:':{WARNLEN}}ValueError in generic_output, nmrprop.dat can be flawed"
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

    def _optrot(self, silent=False):
        """
        calculate optical rotation
        """
        if not self.job["onlyread"]:
            # update control file
            with open(
                os.path.join(self.job["workdir"], "control"), "r", newline=None
            ) as inp:
                tmp = inp.readlines()
            with open(
                os.path.join(self.job["workdir"], "control"), "w", newline=None
            ) as out:
                for line in tmp:
                    if "functional" in line:
                        if (
                            dfa_settings.functionals[self.job["func2"]]["tm"]
                            == "b3lyp-3c"
                        ):
                            out.write(f"   functional b3-lyp   \n")
                        else:
                            out.write(
                                f"   functional {dfa_settings.functionals[self.job['func2']]['tm']}   \n"
                            )
                    elif "$disp" in line:
                        pass
                    else:
                        out.write(line + "\n")

            if not silent:
                print(
                    f"Running optical-rotation calculation in {last_folders(self.job['workdir'], 2):18}"
                )
            files = ["dipl_a", "dipole_a", "rhs_a"]
            for file in files:
                if os.path.isfile(os.path.join(self.job["workdir"], file)):
                    os.remove(os.path.join(self.job["workdir"], file))
            with open(
                os.path.join(self.job["workdir"], "escf.out"), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    [external_paths["escfpath"]],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.job["workdir"],
                    stdout=outputfile,
                    env=ENVIRON,
                )
            time.sleep(0.02)
        # Read output
        # check if scf is converged:
        if os.path.isfile(os.path.join(self.job["workdir"], "escf.out")):
            with open(
                os.path.join(self.job["workdir"], "escf.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                # -------------testnew
                escf_ok = False
                length_rep = False  # (length representation)
                velocity_rep = False  # (velocity representation)
                do_read_length = False
                do_read_velocity = True
                if self.job["func2"] in dfa_settings().hybrid_dfa():
                    do_read_length = True
                dum = 0
                frequencies = self.job["freq_or"]
                frequencies.sort(reverse=True)
                for line in stor:
                    if "escf ended normally" in line:
                        escf_ok = True
                    if " Frequency / nm:  " in line:
                        freq = float(line.strip().split()[-1])
                        if isclose(frequencies[dum], freq, abs_tol=0.6):
                            freq = float(frequencies[dum])
                        else:
                            print("Can't find freq in nm!")
                        dum += 1
                    if " specific rotation [alpha] in deg*[dm(g/cc)]^(-1)" in line:
                        if not length_rep:
                            length_rep = True
                            velocity_rep = False
                        elif length_rep:
                            velocity_rep = True
                            length_rep = False
                        if velocity_rep and do_read_velocity:
                            self.job["energy"] = 0.0
                            self.job["success"] = True
                            self.job["energy2"] = 0.0
                            self.job["erange1"][freq] = float(line.split("(-1)")[-1])
                        elif length_rep and do_read_length:
                            print("Using length representation.")
                            self.job["energy"] = 0.0
                            self.job["success"] = True
                            self.job["energy2"] = 0.0
                            self.job["erange1"][freq] = float(line.split("(-1)")[-1])
                for freq in self.job["freq_or"]:
                    if freq not in self.job["erange1"].keys():
                        escf_ok = False
                if not escf_ok:
                    print(
                        f"{'ERROR:':{WARNLEN}}in escf.out "
                        f"{last_folders(self.job['workdir'], 2):18} not converged!"
                    )
                    self.job["success"] = False
                    self.job["energy"] = 0.0
                    self.job["energy2"] = 0.0
                    self.job["erange1"] = {}
        else:
            print(
                f"{'WARNING:':{WARNLEN}}"
                f"{os.path.join(self.job['workdir'], 'escf.out')} doesn't exist!"
            )
            self.job["success"] = False
            self.job["energy"] = 0.0
            self.job["energy2"] = 0.0
            self.job["erange1"] = {}

    def execute(self):
        """
        Choose what to execute for the jobtype
        """
        try:
            if self.job["jobtype"] == "rrhoxtb":
                self._xtbrrho()
            elif self.job["jobtype"] == "xtb_sp":
                self._xtb_sp()
            elif self.job["jobtype"] == "prep":
                if not self.job["onlyread"]:
                    self._prep_cefine()
                else:
                    self.job["success"] = True
            elif self.job["jobtype"] in ("sp", "sp_implicit"):
                if self.job["prepinfo"]:
                    # do cefine first
                    if not self.job["onlyread"]:
                        self._prep_cefine()
                    else:
                        self.job["success"] = True
                    if not self.job["success"]:
                        return
                self._sp()
            elif self.job["jobtype"] == "cosmors":
                self._cosmors()
            elif self.job["jobtype"] == "xtbopt":
                self._xtbopt()
            elif self.job["jobtype"] == "genericout":
                self._genericoutput()
            elif self.job["jobtype"] in ("couplings", "couplings_sp"):
                if self.job["prepinfo"]:
                    if not self.job["onlyread"]:
                        self._prep_cefine()
                    else:
                        self.job["success"] = True
                    if not self.job["success"]:
                        return
                    if self.job["jobtype"] == "couplings_sp":
                        self._sp(silent=False)
                        if not self.job["success"]:
                            return
                        else:
                            try:
                                tmp_from = os.path.join(self.job["workdir"], "mos")
                                tmp_to = os.path.join(self.job["workdir"], "mos_J")
                                shutil.copy(tmp_from, tmp_to)
                            except FileNotFoundError:
                                pass
                self._nmr_coupling()
            elif self.job["jobtype"] in ("shieldings", "shieldings_sp"):
                if self.job["prepinfo"]:
                    if not self.job["onlyread"]:
                        self._prep_cefine()
                    else:
                        self.job["success"] = True
                    if not self.job["success"]:
                        return
                    if self.job["copymos"]:
                        # use mos as starting mos if basisJ == basisS but
                        # funcJ != funcS
                        try:
                            tmp_from = os.path.join(self.job["workdir"], "mos_J")
                            tmp_to = os.path.join(self.job["workdir"], "mos")
                            shutil.copy(tmp_from, tmp_to)
                        except FileNotFoundError:
                            pass
                        # print("copied mos")
                    if self.job["jobtype"] == "shieldings_sp":
                        self._sp(silent=False)
                        # print("ran sp")
                        if not self.job["success"]:
                            return
                # print("running shieldings")
                self._nmr_shielding()
            elif self.job["jobtype"] in ("gbsa_gsolv", "alpb_gsolv"):
                if self.job["prepinfo"]:
                    tmp_solvent = self.job["solvent"]
                    self.job["solvent"] = "gas"
                    if not self.job["onlyread"]:
                        self._prep_cefine()
                    else:
                        self.job["success"] = True
                    if not self.job["success"]:
                        return
                    self._sp()
                    if not self.job["success"]:
                        return
                    self.job["solvent"] = tmp_solvent
                self._xtb_gsolv()
            elif self.job["jobtype"] in ("opt-rot", "opt-rot_sp"):
                if self.job["prepinfo"]:
                    if not self.job["onlyread"]:
                        self._prep_cefine()
                    else:
                        self.job["success"] = True
                    if not self.job["success"]:
                        return
                if self.job["jobtype"] == "opt-rot_sp":
                    self._sp()
                    if not self.job["success"]:
                        return
                self._optrot()
            else:
                print(f"JOBTYPE {self.job['jobtype']} UNKNOWN!")
        except OSError as error:
            self.job["success"] = False
            print(f"{'IO-ERROR:':{WARNLEN}}in CONF{self.id} calculation")
            print(error, file=sys.stderr)