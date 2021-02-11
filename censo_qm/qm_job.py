"""
Contains QmJob base class for calculating QM related properties of conformers.
Additionally contains functions which should be present irrespective of the QM
code. (xTB always available)
"""
import os
import math

try:
    from math import isclose
except ImportError:
    from .utilities import isclose
import time
import subprocess
import json
from .cfg import ENVIRON, CODING, censo_solvent_db
from .utilities import last_folders, print
from .datastructure import MoleculeData


class QmJob(MoleculeData):
    """
    QmJob base class for calculating QM related properties of conformers.
    """

    def __init__(self, rank, *args, **kwargs):
        MoleculeData.__init__(self, rank, *args, **kwargs)
        self.reset_job_info()

    def reset_job_info(self):
        """
        Clear information/instructions from the previous job
        """
        self.job = {
            "jobtype": "",
            "prepinfo": [],  # additional info for cefine
            "method": "",  # description of the method
            "method2": "",  # description of the method
            "workdir": "",
            "copymos": "",
            "omp": 1,
            "charge": 0,
            "unpaired": 0,
            "gfn_version": None,
            "bhess": None,
            "consider_sym": False,
            "symmetry": "C1",
            "rmsdbias": False,
            "sm_rrho": None,
            "func": None,
            "func2": None,  # functional used in subsequent property calculation
            "basis": None,
            "solvent": "gas",
            "sm": "gas-phase",
            "rmsd": 0.0,  # rmsd in case of bhess
            "nat": None,  # number of atoms
            "onlyread": False,  # don't calculate, just perform readout
            "progpath": "",
            "xtb_driver_path": "",  # program path to xtb if xtb as driver
            # optimization related:
            "optcycles": None,  # number of cycles that are allowed in the
            "hlow": None,  # setting for ancopt
            "optlevel": None,
            # geometry optimization
            "cycles": 0,  # number of cycles it needed for optimization convergence
            "ecyc": [],
            "decyc": [],
            "grad_norm": 10.0,
            "converged": False,
            # temperature related:
            "trange": [],  # list with temperatures to evaluate G,H,S
            "temperature": 298.15,
            # nmrprop related:
            "h_active": False,
            "c_active": False,
            "f_active": False,
            "p_active": False,
            "si_active": False,
            # optical rotation related:
            "freq_or": [],
            # return values which can be updated:
            "success": False,
            "energy": 0.0,
            "energy2": 0.0,
            "erange1": {},
            "erange2": {},
            "erange3": {},
            "errormessage": [],
            "internal_error": [],
            #
            "cosmorsparam": "",  # normal/fine
        }

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

    def _genericoutput(self):
        """
        Read shielding and coupling constants and write them to plain output
        The first natom lines contain the shielding constants, and from
        line natom +1 the coupling constants are written.
        """
        pass

    def _xtb_sp(self):
        """
        Get single-point energy from xTB
        """
        files = [
            "xtbrestart",
            "xtbtopo.mol",
            "xcontrol-inp",
            "wbo",
            "charges",
            "gfnff_topo",
            "sp.out",
        ]
        for file in files:
            if os.path.isfile(os.path.join(self.job["workdir"], file)):
                os.remove(os.path.join(self.job["workdir"], file))
        # run single-point:

        call = [
            self.job["progpath"],
            "coord",
            "--" + self.job["gfn_version"],
            "--sp",
            "--chrg",
            str(self.job["charge"]),
            "--norestart",
        ]
        if self.job["solvent"] != "gas":
            call.extend(
                [
                    "--" + str(self.job["sm"]),
                    censo_solvent_db[self.job["solvent"]]["xtb"][1],
                ]
            )

        with open(
            os.path.join(self.job["workdir"], "sp.out"), "w", newline=None
        ) as outputfile:
            returncode = subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.job["workdir"],
                stdout=outputfile,
                env=ENVIRON,
            )
        if returncode != 0:
            self.job["energy"] = 0.0
            self.job["success"] = False
            print(
                f"ERROR: {self.job['gfn_version']}-xTB error in "
                f"{last_folders(self.job['workdir'], 2)}"
            )
            return
        # read gas phase energy:
        if os.path.isfile(os.path.join(self.job["workdir"], "sp.out")):
            with open(
                os.path.join(self.job["workdir"], "sp.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                store = inp.readlines()
                for line in store:
                    if "| TOTAL ENERGY" in line:
                        try:
                            self.job["energy"] = float(line.split()[3])
                            self.job["success"] = True
                        except:
                            print(
                                "Error while converting "
                                "single-point in: {}".format(
                                    last_folders(self.job["workdir"], 2)
                                )
                            )
                            self.job["energy"] = 0.0
                            self.job["success"] = False
                            return
        else:
            self.job["energy"] = 0.0
            self.job["success"] = False
            print(
                f"ERROR: {self.job['gfn_version']}-xTB error in "
                f"{last_folders(self.job['workdir'], 2)}"
            )
            return
        print(
            f"{self.job['gfn_version']}-xTB energy for {last_folders(self.job['workdir'], 2)}"
            f" = {self.job['energy']: >.7f}"
        )

    def _xtb_gsolv(self):
        """
        Calculate additive GBSA or ALPB solvation contribution by
        Gsolv = Esolv - Egas,
        using xTB and the GFNn or GFN-FF hamiltonian.
        --> return gsolv at energy2
        """
        tmp_gas = 0
        tmp_solv = 0
        if self.job["jobtype"] == "gbsa_gsolv":
            xtbsm = "gbsa"
        elif self.job["jobtype"] == "alpb_gsolv":
            xtbsm = "alpb"
        print(
            f"Running {self.job['jobtype'].upper()} calculation in "
            f"{last_folders(self.job['workdir'], 3)}"
        )
        files = [
            "xtbrestart",
            "xtbtopo.mol",
            "xcontrol-inp",
            "wbo",
            "charges",
            "gfnff_topo",
            "gas.out",
            "solv.out",
        ]
        for file in files:
            if os.path.isfile(os.path.join(self.job["workdir"], file)):
                os.remove(os.path.join(self.job["workdir"], file))
        # run gas phase single-point:
        with open(
            os.path.join(self.job["workdir"], "gas.out"), "w", newline=None
        ) as outputfile:
            returncode = subprocess.call(
                [
                    self.job["xtb_driver_path"],
                    "coord",
                    "--" + self.job["gfn_version"],
                    "--sp",
                    "--chrg",
                    str(self.job["charge"]),
                    "--norestart",
                ],
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.job["workdir"],
                stdout=outputfile,
                env=ENVIRON,
            )
        if returncode != 0:
            self.job["energy2"] = 0.0
            self.job["success"] = False
            print(
                f"ERROR: Gas phase {self.job['gfn_version']}-xTB error in "
                f"{last_folders(self.job['workdir'], 3)}"
            )
            return
        # read gas phase energy:
        if os.path.isfile(os.path.join(self.job["workdir"], "gas.out")):
            with open(
                os.path.join(self.job["workdir"], "gas.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                store = inp.readlines()
                for line in store:
                    if "| TOTAL ENERGY" in line:
                        try:
                            tmp_gas = float(line.split()[3])
                            self.job["success"] = True
                        except:
                            print(
                                "Error while converting gas phase "
                                "single-point in: {}".format(
                                    last_folders(self.job["workdir"], 3)
                                )
                            )
                            tmp_gas = None
                            self.job["energy2"] = 0.0
                            self.job["success"] = False
                            return
        else:
            self.job["energy2"] = 0.0
            self.job["success"] = False
            print(
                f"ERROR: Gas phase {self.job['gfn_version']}-xTB error in "
                f"{last_folders(self.job['workdir'], 3)}"
            )
            return
        # run single-point in solution:
        # ``reference'' corresponds to 1\;bar of ideal gas and 1\;mol/L of liquid
        #   solution at infinite dilution,
        with open(
            os.path.join(self.job["workdir"], "solv.out"), "w", newline=None
        ) as outputfile:
            returncode = subprocess.call(
                [
                    self.job["xtb_driver_path"],
                    "coord",
                    "--" + self.job["gfn_version"],
                    "--sp",
                    "--" + xtbsm,
                    censo_solvent_db[self.job["solvent"]]["xtb"][1],
                    "reference",
                    "--chrg",
                    str(self.job["charge"]),
                    "--norestart",
                ],
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.job["workdir"],
                stdout=outputfile,
                env=ENVIRON,
            )
            if returncode != 0:
                self.job["energy2"] = 0.0
                self.job["success"] = False
                print(
                    f"ERROR: Solution phase {self.job['gfn_version']}-xTB error in "
                    f"{last_folders(self.job['workdir'], 3)}"
                )
                return
        time.sleep(0.05)
        # #read solv.out
        if os.path.isfile(os.path.join(self.job["workdir"], "solv.out")):
            with open(
                os.path.join(self.job["workdir"], "solv.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                store = inp.readlines()
                for line in store:
                    if "| TOTAL ENERGY" in line:
                        try:
                            tmp_solv = float(line.split()[3])
                            self.job["success"] = True
                        except:
                            print(
                                "Error while converting solution phase "
                                "single-point in: {}".format(
                                    last_folders(self.job["workdir"], 3)
                                )
                            )
                            tmp_solv = None
                            self.job["energy2"] = 0.0
                            self.job["success"] = False
                            return
        else:
            self.job["energy2"] = 0.0
            self.job["success"] = False
            print(
                f"ERROR: Solution phase {self.job['gfn_version']}-xTB error in "
                f"{last_folders(self.job['workdir'], 3)}"
            )
            return
        if self.job["success"]:
            if tmp_solv is None or tmp_gas is None:
                self.job["energy2"] = 0.0
                self.job["success"] = False
            else:
                self.job["energy2"] = tmp_solv - tmp_gas
                self.job["success"] = True
                self.job["energy_xtb_gas"] = tmp_gas
                self.job["energy_xtb_solv"] = tmp_solv

    def _xtbrrho(self):
        """
        mRRHO contribution with GFNn/GFN-FF-XTB
        """
        if not self.job["onlyread"]:
            print(
                f"Running {str(self.job['gfn_version']).upper()}-xTB mRRHO in "
                f"{last_folders(self.job['workdir'], 2)}"
            )
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
            if self.job["trange"]:
                for t in list(self.job["trange"]):
                    if isclose(self.job["temperature"], t, abs_tol=0.6):
                        self.job["trange"].pop(self.job["trange"].index(t))
                self.job["trange"].append(self.job["temperature"])
            with open(
                os.path.join(self.job["workdir"], "xcontrol-inp"), "w", newline=None
            ) as xcout:
                xcout.write("$thermo\n")
                if self.job["trange"]:
                    xcout.write(
                        f"    temp="
                        f"{','.join([str(i) for i in self.job['trange']])}\n"
                    )
                else:
                    xcout.write("    temp={}\n".format(self.job["temperature"]))
                if self.job.get("sthr", "automatic") == 'automatic':
                    xcout.write("    sthr=50.0\n")
                else:
                    xcout.write("    sthr={}\n".format(self.job['sthr']))
                if self.job.get("imagthr", "automatic") == 'automatic':
                    if self.job["bhess"]:
                        xcout.write("    imagthr={}\n".format("-100"))
                    else:
                        xcout.write("    imagthr={}\n".format("-50"))
                else:
                    xcout.write("    imagthr={}\n".format(self.job["imagthr"]))
                if self.job.get("scale", "automatic") == 'automatic':
                    #xcout.write("    scale={}\n".format("1.0"))
                    pass # is method dependant leave it to xTB e.g. GFNFF has a 
                    # different scaling factor than GFN2
                else:
                    xcout.write("    scale={}\n".format(self.job['scale']))
                xcout.write("$symmetry\n")
                if self.job["consider_sym"]:
                    # xcout.write("    desy=0.1\n") # taken from xtb defaults
                    xcout.write("     maxat=1000\n")
                else:
                    xcout.write("    desy=0.0\n")
                xcout.write("$end\n")
            if self.job["bhess"]:
                # set ohess or bhess
                dohess = "--bhess"
                olevel = "vtight"
            else:
                dohess = "--ohess"
                olevel = "vtight"
            time.sleep(0.05)
            with open(
                os.path.join(self.job["workdir"], "ohess.out"), "w", newline=None
            ) as outputfile:
                if self.job["solvent"] != "gas":
                    callargs = [
                        self.job["progpath"],
                        "coord",
                        "--" + str(self.job["gfn_version"]),
                        dohess,
                        olevel,
                        "--" + str(self.job["sm_rrho"]),
                        censo_solvent_db[self.job["solvent"]]["xtb"][1],
                        "--chrg",
                        str(self.job["charge"]),
                        "--enso",
                        "--norestart",
                        "-I",
                        "xcontrol-inp",
                    ]
                else:
                    callargs = [
                        self.job["progpath"],
                        "coord",
                        "--" + str(self.job["gfn_version"]),
                        dohess,
                        olevel,
                        "--chrg",
                        str(self.job["charge"]),
                        "--enso",
                        "--norestart",
                        "-I",
                        "xcontrol-inp",
                    ]
                if self.job["rmsdbias"]:
                    callargs.extend(
                        [
                            "--bias-input",
                            str(os.path.join(self.job["cwd"], "rmsdpot.xyz")),
                        ]
                    )
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
            time.sleep(0.05)
            # check if converged:
            if returncode != 0:
                self.job["energy"] = 0.0
                self.job["success"] = False
                self.job["errormessage"].append(
                    f"ERROR: {str(self.job['gfn_version']).upper()}-xTB ohess error in "
                    f"{last_folders(self.job['workdir'], 2):18}"
                )
                print(self.job["errormessage"][-1])
                return
        # start reading output!
        if self.job["trange"]:
            if not os.path.isfile(os.path.join(self.job["workdir"], "ohess.out")):
                self.job["energy"] = 0.0
                self.job["success"] = False
                self.job["errormessage"].append(
                    f"ERROR: file {self.job['workdir']}/ohess.out could not be found!"
                )
                print(self.job["errormessage"][-1])
                return
            gt = {}
            ht = {}
            rotS = {}
            with open(
                os.path.join(self.job["workdir"], "ohess.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                store = inp.readlines()
                # rotational entropy:
                for line in store:
                    if "VIB" in line:
                        try:
                            realline = store.index(line) + 1
                            T = float(line.split()[0])
                            rotS[T] = float(store[realline].split()[4])
                        except (KeyError, ValueError):
                            pass
                for line in store:
                    if "T/K" in line:
                        start = store.index(line)
                for line in store[start + 2 :]:
                    if "----------------------------------" in line:
                        break
                    else:
                        try:
                            T = float(line.split()[0])
                            gt[T] = float(line.split()[4])
                            ht[T] = float(line.split()[2])
                        except (ValueError, KeyError):
                            print("ERROR: can not convert G(T)")
            if len(self.job["trange"]) == len(gt):
                self.job["success"] = True
                self.job["erange1"] = gt
                self.job["erange2"] = ht
                self.job["erange3"] = rotS
            else:
                self.job["success"] = False
                return
        # end self.trange
        if self.job["bhess"]:
            # read rmsd_info
            with open(
                os.path.join(self.job["workdir"], "ohess.out"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                store = inp.readlines()
                for line in store:
                    if "final rmsd / " in line:
                        try:
                            self.job["rmsd"] = float(line.split()[3])
                        except (ValueError):
                            self.job["rmsd"] = 0.0
        if os.path.isfile(os.path.join(self.job["workdir"], "xtb_enso.json")):
            with open(
                os.path.join(self.job["workdir"], "xtb_enso.json"),
                "r",
                encoding=CODING,
                newline=None,
            ) as f:
                data = json.load(f)
            if "number of imags" in data:
                if data["number of imags"] > 0:
                    print(
                        f"WARNING: found {data['number of imags']} significant"
                        f" imaginary frequencies in "
                        f"{last_folders(self.job['workdir'], 2)}"
                    )
            if "G(T)" in data:
                if float(self.job["temperature"]) == 0:
                    self.job["energy"] = data.get("ZPVE", 0.0)
                    self.job["success"] = True
                    self.job["erange1"][self.job["temperature"]] = data.get("ZPVE", 0.0)
                    self.job["erange2"][self.job["temperature"]] = data.get("ZPVE", 0.0)
                else:
                    self.job["energy"] = data.get("G(T)", 0.0)
                    self.job["erange1"][self.job["temperature"]] = data.get("G(T)", 0.0)
                    # self.job['erange2'][self.job['temperature']] =  data.get("H(T)", 0.0)
                    self.job["success"] = True
                if "point group" in data:
                    self.job["symmetry"] = data["point group"]
            else:
                print(
                    "Error while converting mRRHO in: {}".format(
                        last_folders(self.job["workdir"], 2)
                    )
                )
                self.job["energy"] = 0.0
                self.job["success"] = False
        else:
            print(
                "WARNING: File {} doesn't exist!".format(
                    os.path.join(self.job["workdir"], "xtb_enso.json")
                )
            )
            self.job["energy"] = 0.0
            self.job["success"] = False

    def execute(self):
        """
        Chooses which function to execute based on jobtype.
        """
        pass
