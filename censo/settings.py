import ast
import json
import os
import shutil
import sys
import configparser
from argparse import Namespace
from dataclasses import dataclass
from functools import reduce
from multiprocessing import Lock
from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Type, Union

from censo.cfg import (ASSETS_PATH, DIGILEN, GSOLV_MODS, PARTS, PLENGTH, PROGS, SOLV_MODS,
                            WARNLEN, GRIDOPTIONS, GFNOPTIONS, CENSORCNAME, __version__)
from censo.errorswarnings import LogicError, LogicWarning


class DfaSettings:
    def __init__(self, obj: Dict[str, Dict[str, Dict]]):
        self.__dfa_dict = obj

    
    def find_func(self, part: str, prog=None):
        """
        return all functionals available for a certain part and (optionally) program
        """
        # TODO - turn into filter using filterfunction defined within find_func
        tmp = []
        for k, v in self.__dfa_dict["functionals"].items():
            if part in v["part"]:
                if prog is None:
                    tmp.append(k)
                else:
                    if v[prog] != "":
                        tmp.append(k)
        
        return tmp


    def get_name(self, func: str, prog: str):
        """
        return the name of a certain functional in the given qm program
        """
        return self.__dfa_dict["functionals"][func][prog]


    def get_disp(self, func: str):
        """
        return the dispersion correction of a certain functional
        """
        return self.__dfa_dict["functionals"][func]['disp']


    def get_type(self, func: str):
        """
        return the type of a certain functional
        """
        return self.__dfa_dict["functionals"][func]["type"]


    @property
    def functionals(self) -> Dict[str, Dict]:
        return self.__dfa_dict["functionals"]


class CensoSettings:
    """
    TODO
    """
    
    cosmors_param = {
        "12-normal": "BP_TZVP_C30_1201.ctd",
        "13-normal": "BP_TZVP_C30_1301.ctd",
        "14-normal": "BP_TZVP_C30_1401.ctd",
        "15-normal": "BP_TZVP_C30_1501.ctd",
        "16-normal": "BP_TZVP_C30_1601.ctd",
        "17-normal": "BP_TZVP_C30_1701.ctd",
        "18-normal": "BP_TZVP_18.ctd",
        "19-normal": "BP_TZVP_19.ctd",
        "12-fine": "BP_TZVPD_FINE_HB2012_C30_1201.ctd",
        "13-fine": "BP_TZVPD_FINE_HB2012_C30_1301.ctd",
        "14-fine": "BP_TZVPD_FINE_C30_1401.ctd",
        "15-fine": "BP_TZVPD_FINE_C30_1501.ctd",
        "16-fine": "BP_TZVPD_FINE_C30_1601.ctd",
        "17-fine": "BP_TZVPD_FINE_C30_1701.ctd",
        "18-fine": "BP_TZVPD_FINE_18.ctd",
        "19-fine": "BP_TZVPD_FINE_19.ctd",
    }

    # load up all resources to set value options in _settings_options
    # TODO - catch error if dfa_settings cannot be created
    try:
        with open(os.path.join(ASSETS_PATH, "censo_dfa_settings.json"), "r") as dfa_file:
            dfa_settings = DfaSettings(json.load(dfa_file))

        with open(os.path.join(ASSETS_PATH, "basis_sets.json"), "r") as bs_file:
            basis_sets = tuple(json.load(bs_file))

        with open(os.path.join(ASSETS_PATH, "censo_solvents_db.json"), "r") as solv_file:
            solvents_db = json.load(solv_file)
    except FileNotFoundError:
        print(f"{'ERROR:':{WARNLEN}}Could not find DFA/basis/solvents file!")
        print("\nGoing to exit!")
        sys.exit(1)

    solv_mods = reduce(lambda x, y: x + y, SOLV_MODS.values())
    gsolv_mods = reduce(lambda x, y: x + y, GSOLV_MODS.values())

    # Available settings are defined here 
    _settings_options = {
        "paths": {
            "orcapath": {
                "default": "",
                "options": [],
            },
            "orcaversion": {
                "default": "",
                "options": []
            },
            "xtbpath": {
                "default": "",
                "options": []
            },
            "crestpath": {
                "default": "",
                "options": []
            },
            "cosmorssetup": {
                "default": "",
                "options": []
            },
            "dbpath": {
                "default": "",
                "options": []
            },
            "cosmothermversion": {
                "default": "",
                "options": []
            },
            "mpshiftpath": {
                "default": "",
                "options": []
            },
            "escfpath": {
                "default": "",
                "options": []
            },
        },
        # TODO - charge, unpaired should probably be removed from here and be given for each run specifically
        "general": {
            "maxprocs": {
                "default": 1,
                "range": [
                    1,
                    128
                ]
            },
            "omp": {
                "default": 1,
                "range": [
                    1,
                    256
                ]
            },
            "imagethr": {
                "default": -100.0,
                "range": [
                    -300.0,
                    0.0
                ]
            },
            "sthr": {
                "default": 0.0,
                "range": [
                    0.0,
                    100.0
                ]
            },
            "scale": {
                "default": 1.0,
                "range": [
                    0.0,
                    1.0
                ]
            },
            "temperature": {
                "default": 298.15,
                "range": [
                    1e-05,
                    2000.0
                ]
            },
            "solvent": {
                "default": "h2o",
                "options": [k for k in solvents_db.keys()]
            },
            "sm_rrho": {
                "default": "alpb",
                "options": [
                    "alpb",
                    "gbsa"
                ]
            },
            "cosmorsparam": {
                "default": "automatic",
                "options": [k for k in cosmors_param.keys()]
            },
            "multitemp": {
                "default": True
            },
            "evaluate_rrho": {
                "default": True
            },
            "consider_sym": {
                "default": True
            },
            "bhess": {
                "default": True
            },
            "rmsdbias": {
                "default": False
            },
            "progress": {
                "default": False
            },
            "check": {
                "default": True
            },
            "balance": {
                "default": True
            },
            "vapor_pressure": {
                "default": False
            },
            "nmrmode": {
                "default": False
            },
            "gas-phase": {
                "default": False
            },
            "trange": {
                "default": [
                    273.15,
                    378.15,
                    5
                ]
            }
        },
        "prescreening": {
            "threshold": {
                "default": 4.0,
                "range": [
                    1.0,
                    10.0
                ]
            },
            "func": {
                "default": "pbe-d4",
                "options": dfa_settings.find_func("prescreening")
            },
            "basis": {
                "default": "def2-SV(P)",
                "options": basis_sets
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "grid": {
                "default": "low",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": True
            },
            "gcp": {
                "default": True
            }
        },
        "screening": {
            "threshold": {
                "default": 3.5,
                "range": [
                    0.75,
                    7.5
                ]
            },
            "func": {
                "default": "r2scan-3c",
                "options": dfa_settings.find_func("screening")
            },
            "basis": {
                "default": "def2-TZVP",
                "options": basis_sets
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "smgsolv": {
                "default": "smd_gsolv",
                "options": gsolv_mods
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "grid": {
                "default": "low+",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": True
            },
            "gcp": {
                "default": True
            }
        },
        "optimization": {
            "optcycles": {
                "default": 8,
                "range": [
                    1,
                    100
                ]
            },
            "radsize": {
                "default": 10,
                "range": [
                    1,
                    100
                ]
            },
            "threshold": {
                "default": 2.5,
                "range": [
                    0.5,
                    5
                ]
            },
            "hlow": {
                "default": 0.01,
                "range": [
                    0.01,
                    1.0
                ]
            },
            "optimization_P_threshold": {
                "default": 99.0,
                "range": [
                    1.0,
                    10.0
                ]
            },
            "spearmanthr": {
                "default": 0.0,
                "range": [
                    0.0,
                    10.0
                ]
            },
            "func": {
                "default": "r2scan-3c",
                "options": dfa_settings.find_func("optimization")
            },
            "basis": {
                "default": "def2-TZVP",
                "options": basis_sets
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "prog2opt": {
                "default": "prog",
                "options": []
            },
            "sm": {
                "default": "smd",
                "options": solv_mods
            },
            "smgsolv": {
                "default": "smd_gsolv",
                "options": gsolv_mods
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "optlevel2": {
                "default": "automatic",
                "options": [
                    "crude",
                    "sloppy",
                    "loose",
                    "lax",
                    "normal",
                    "tight",
                    "vtight",
                    "extreme",
                    "automatic"
                ]
            },
            "grid": {
                "default": "high",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": True
            },
            "gcp": {
                "default": True
            },
            "opt_spearman": {
                "default": True
            },
            "crestcheck": {
                "default": False
            }
        },
        "refinement": {
            "threshold": {
                "default": 99.0,
                "min": [
                    1.0,
                    10.0
                ]
            },
            "prog": {
                "default": "orca",
                "options": PROGS
            },
            "func": {
                "default": "wb97x-v",
                "options": dfa_settings.find_func("refinement")
            },
            "basis": {
                "default": "def2-TZVPP",
                "options": basis_sets
            },
            "smgsolv": {
                "default": "smd_gsolv",
                "options": gsolv_mods
            },
            "gfnv": {
                "default": "gfn2",
                "options": GFNOPTIONS
            },
            "grid": {
                "default": "high+",
                "options": GRIDOPTIONS
            },
            "run": {
                "default": False
            },
            "gcp": {
                "default": True
            }
        },
        "nmr": {
            "resonance_frequency": {
                "default": 300.0,
                "range": [
                    150.0,
                    1000.0
                ]
            },
            "prog4_j": {
                "default": "tm",
                "options": PROGS
            },
            "func_j": {
                "default": "pbe0-d4",
                "options": []
            },
            "basis_j": {
                "default": "def2-TZVP",
                "options": basis_sets
            },
            "sm4_j": {
                "default": "default",
                "options": solv_mods
            },
            "prog4_s": {
                "default": "tm",
                "options": PROGS
            },
            "func_s": {
                "default": "pbe0-d4",
                "options": []
            },
            "basis_s": {
                "default": "def2-TZVP",
                "options": basis_sets
            },
            "sm4_s": {
                "default": "default",
                "options": solv_mods
            },
            "h_ref": {
                "default": "TMS",
                "options": [
                    "TMS"
                ]
            },
            "c_ref": {
                "default": "TMS",
                "options": [
                    "TMS"
                ]
            },
            "f_ref": {
                "default": "CFCl3",
                "options": [
                    "CFCl3"
                ]
            },
            "si_ref": {
                "default": "TMS",
                "options": [
                    "TMS"
                ]
            },
            "p_ref": {
                "default": "TMP",
                "options": [
                    "TMP",
                    "PH3"
                ]
            },
            "run": {
                "default": False
            },
            "couplings": {
                "default": True
            },
            "shieldings": {
                "default": True
            },
            "h_active": {
                "default": True
            },
            "c_active": {
                "default": True
            },
            "f_active": {
                "default": False
            },
            "si_active": {
                "default": False
            },
            "p_active": {
                "default": False
            }
        },
        "optrot": {
            "func": {
                "default": "pbe-d4",
                "options": dfa_settings.find_func("optrot")
            },
            "func_or_scf": {
                "default": "r2scan-3c",
                "options": []
            },
            "basis": {
                "default": "def2-SVPD",
                "options": basis_sets
            },
            "prog": {
                "default": "orca",
                "options": [
                    "orca"
                ]
            },
            "run": {
                "default": False
            },
            "freq_or": {
                "default": [
                    598.0
                ]
            }
        },
        "uvvis": {
            "nroots": {
                "default": 20,
                "range": [
                    1,
                    100
                ]
            },
            "sigma": {
                "default": 0.1,
                "range": [
                    0.1,
                    1.0
                ]
            },
            "run": {
                "default": False
            }
        },
    }
    # TODO - defaults for uvvis (func: wB97X-D3)
    # TODO - find solutions for "automatic" or "default"/"prog"/"whatever" cases
    # TODO - rename options but remember old names for later mapping (backwards compatibility)    # removed general:func and general:prog, added part1:func1/prog1, part2:func2/prog2
    # moved scale, imagthr, sthr to float type
    # removed general:nconf (no setting, rather run-specific info)
    # added general:vapor_pressure, general:gas-phase
    # renamed part-enabling settings (e.g. "part1") to "run"
    # generalized func, basis, prog (except part4 TODO)
    # contains settings sorted by type, with defaults and limits
    # mapping to immutable dict type (MappingProxyType)
    """
    settings_options: MappingProxyType
    |-int: type, MappingProxyType
    |  |-general: str, MappingProxyType
    |  |  |-> charge: str, dict
    |  |  |   |-> default: int
    |  |  |   |-> ~options
    |  |  |   
    |  |  |-> ...
    |  | 
    |  |-prescreening: str
    |  |  |-> ...
    |  | 
    | ... 
    | 
    |-float: type, MappingProxyType
    |  |-part: str, MappingProxyType
    |  |  |-> setting: str, dict
    .....
    """    
    # TODO - solvent mapping
    # still included as fallback
    """settings_options = MappingProxyType({
        int: MappingProxyType({
            "general": MappingProxyType({
                "charge": {"default": 0, "range": (-10, 10)}, # TODO - (re)move?
                "unpaired": {"default": 0, "range": (0, 14)}, # TODO - (re)move?
                "maxprocs": {"default": 1, "range": (1, 128)}, # number of processes
                "omp": {"default": 1, "range": (1, 256)}, # number of cores per process
            }),
            "prescreening": None,
            "screening": None,
            "optimization": MappingProxyType({
                "optcycles": {"default": 8, "range": (1, 100)}, # TODO - which value as min?
                "radsize": {"default": 10, "range": (1, 100)}, # TODO - same
            }),
            "refinement": None,
            "nmr": None,
            "optrot": None,
            "uvvis": MappingProxyType({
                "nroots": {"default": 20, "range": (1, 100)}, # TODO - set dynamically
            }),
        }),
        float: MappingProxyType({
            "general": MappingProxyType({
                "imagethr": {"default": -100.0, "range": (-300.0, 0.0)}, # TODO - threshold for imaginary frequencies
                "sthr": {"default": 0.0, "range": (0.0, 100.0)}, # TODO - what is this?
                "scale": {"default": 1.0, "range": (0.0, 1.0)}, # TODO - what is this?
                "temperature": {"default": 298.15, "range": (0.00001, 2000.0)}, # TODO
            }),
            "prescreening": MappingProxyType({
                "threshold": {"default": 4.0, "range": (1.0, 10.0)}, # TODO - which value as min?
            }),
            "screening": MappingProxyType({
                "threshold": {"default": 3.5, "range": (0.75, 7.5)},
            }),
            "optimization": MappingProxyType({
                "threshold": {"default": 2.5, "range": (0.5, 5)}, # TODO - rename?
                "hlow": {"default": 0.01, "range": (0.01, 1.0)}, # TODO
                "optimization_P_threshold": {"default": 99.0, "range": (1.0, 10.0)}, # TODO
                "spearmanthr": {"default": 0.0, "range": (0.0, 10.0)},
            }),
            "refinement": MappingProxyType({
                "threshold": {"default": 99.0, "range": (1.0, 10.0)}, # TODO
            }),
            "nmr": MappingProxyType({
                "resonance_frequency": {"default": 300.0, "range": (150.0, 1000.0)}, # TODO
            }),
            "optrot": None,
            "uvvis": MappingProxyType({
                "sigma": {"default": 0.1, "range": (0.1, 1.0)},
            }),
        }),
        str: MappingProxyType({
            "general": MappingProxyType({
                # TODO - removed gas here, since there is already 'gas-phase' option (also doesn't make sense)
                "solvent": {"default": "h2o", "options": tuple([k for k in solvents_db.keys()])}, 
                #"prog_rrho": {"default": "xtb", "options": ("xtb")}, # TODO - keep this here?
                "sm_rrho": {"default": "alpb", "options": ("alpb", "gbsa")}, # TODO - same
                "cosmorsparam": {"default": "automatic", "options": tuple([k for k in cosmors_param.keys()])},
            }),
            "prescreening": MappingProxyType({
                "func": {"default": "pbe-d4", "options": tuple(dfa_settings.find_func("prescreening"))},
                "basis": {"default": "def2-SV(P)", "options": ("automatic",) + tuple(dfa_settings.composite_bs) + ("def2-SV(P)", "def2-TZVP")},
                "prog": {"default": "orca", "options": PROGS},
                "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                "grid": {"default": "low", "options": GRIDOPTIONS},
            }),
            "screening": MappingProxyType({
                "func": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func1"))},
                "basis": {"default": "automatic", "options": ("automatic",) + basis_sets},
                "prog": {"default": "orca", "options": PROGS},
                "smgsolv": {"default": "smd_gsolv", "options": gsolv_mods},
                "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                "grid": {"default": "low+", "options": GRIDOPTIONS},
            }),
            "optimization": MappingProxyType({
                "func": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func2"))},
                "basis": {"default": "automatic", "options": ("automatic",) + basis_sets},
                "prog": {"default": "orca", "options": PROGS},
                "prog2opt": {"default": "prog", "options": PROGS}, # TODO - ??? prog2 ??? # FIXME
                "sm": {"default": "smd", "options": solv_mods}, # FIXME
                "smgsolv": {"default": "smd_gsolv", "options": gsolv_mods},
                "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                "optlevel2": {"default": "automatic", "options": ("crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme", "automatic")}, # TODO - what does this mean?
                "grid": {"default": "high", "options": GRIDOPTIONS},
            }),
            "refinement": MappingProxyType({
                "prog": {"default": "orca", "options": PROGS},
                "func": {"default": "wb97x-v", "options": tuple(dfa_settings.find_func("func3"))},
                "basis": {"default": "def2-TZVPP", "options": basis_sets},
                "smgsolv": {"default": "smd_gsolv", "options": gsolv_mods},
                "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                "grid": {"default": "high+", "options": GRIDOPTIONS},
            }),
            "nmr": MappingProxyType({
                "prog4_j": {"default": "tm", "options": PROGS},
                "func_j": {"default": "pbe0-d4", "options": tuple(dfa_settings.find_func("func_j"))},
                "basis_j": {"default": "def2-TZVP", "options": basis_sets},
                "sm4_j": {"default": "default", "options": solv_mods}, # FIXME
                "prog4_s": {"default": "tm", "options": PROGS},
                "func_s": {"default": "pbe0-d4", "options": tuple(dfa_settings.find_func("func_s"))},
                "basis_s": {"default": "def2-TZVP", "options": basis_sets},
                "sm4_s": {"default": "default", "options": solv_mods}, # FIXME
                "h_ref": {"default": "TMS", "options": ("TMS",)},
                "c_ref": {"default": "TMS", "options": ("TMS",)},
                "f_ref": {"default": "CFCl3", "options": ("CFCl3",)},
                "si_ref": {"default": "TMS", "options": ("TMS",)},
                "p_ref": {"default": "TMP", "options": ("TMP", "PH3")},
            }),
            "optrot": MappingProxyType({
                "func": {"default": "pbe-d4", "options": tuple(dfa_settings.find_func("func_or"))},
                "func_or_scf": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func_or_scf"))},
                "basis": {"default": "def2-SVPD", "options": basis_sets},
                "prog": {"default": "orca", "options": ("orca",)},
            }),
            "uvvis": None,
        }),
        bool: MappingProxyType({
            "general": MappingProxyType({
                "multitemp": {"default": True},
                "evaluate_rrho": {"default": True},
                "consider_sym": {"default": True},
                "bhess": {"default": True},
                "rmsdbias": {"default": False},
                "progress": {"default": False},
                "check": {"default": True},
                "balance": {"default": True},
                "vapor_pressure": {"default": False},
                "nmrmode": {"default": False},
                "gas-phase": {"default": False},
            }),
            "prescreening": MappingProxyType({
                "run": {"default": True},
                "gcp": {"default": True},
            }),
            "screening": MappingProxyType({
                "run": {"default": True},
                "gcp": {"default": True},
            }),
            "optimization": MappingProxyType({
                "run": {"default": True},
                "gcp": {"default": True},
                # "ancopt": {"default": True},
                "opt_spearman": {"default": True},
                "crestcheck": {"default": False},
            }),
            "refinement": MappingProxyType({
                "run": {"default": False},
                "gcp": {"default": True},
            }),
            "nmr": MappingProxyType({
                "run": {"default": False},
                "couplings": {"default": True},
                "shieldings": {"default": True},
                "h_active": {"default": True},
                "c_active": {"default": True},
                "f_active": {"default": False},
                "si_active": {"default": False},
                "p_active": {"default": False},
            }),
            "optrot": MappingProxyType({
                "run": {"default": False},
            }),
            "uvvis":MappingProxyType({
                "run": {"default": False},
            }),
        }),
        list: MappingProxyType({
            "general": MappingProxyType({
                "trange": {"default": [273.15, 378.15, 5]},
            }),
            "prescreening": None,
            "screening": None,
            "optimization": None,
            "refinement": None,
            "nmr": None,
            "optrot": MappingProxyType({
                "freq_or": {"default": [598.0]},
            }),
            "uvvis": None,
        }),
    })"""

    
    @classmethod
    def get_type(cls, section: str, name: str) -> Type:
        """
        Get the type of the given setting.

        Args:
            section (str): The section of the setting.
            name (str): The name of the setting.

        Returns:
            Type: The type of the setting.
        """
        return type(cls._settings_options[section][name]["default"])
    
    def __init__(self):
        """
        Handles program configuration
        On creating a new CensoSettings object it will try to read the rcfile
        If there is no rcfile it will create a new one automatically
        """ 
        # stores all settings, grouped by section (e.g. "general", "prescreening", ...)
        # every string is in lower case
        self.__settings_current: Dict[str, Dict[str, Any]]
        
        # absolute path to configuration file, try to find .censorc on construction
        self.censorc_path: str = self.__find_rcfile()

        # if no rcfile is found create a new one in home directory
        if self.censorc_path is None: 
            self.censorc_path = os.path.join(os.path.expanduser("~"), CENSORCNAME)
            self.write_rcfile(self.censorc_path)
        
        # read config file
        self.__settings_current = self.__read_rcfile()
    
    
    def print_paths(self) -> None:
        """
        print out paths of all external qm programs
        """
        lines = []
        
        lines.append("\n" + "".ljust(PLENGTH, "-") + "\n")
        lines.append("PATHS of external QM programs".center(PLENGTH, " ") + "\n")
        lines.append("".ljust(PLENGTH, "-") + "\n")
        
        for program, path in self.__settings_current["paths"].items():
            lines.append(f"{program}:".ljust(DIGILEN, " ") + f"{path}\n")
            
        for line in lines:
            print(line)      


    @property
    def settings_current(self) -> Dict[str, Dict[str, Any]]:
        """
        returns the complete __settings_current
        """
        return self.__settings_current

    
    @settings_current.setter
    def settings_current(self, settings_dict: Dict[str, Dict[str, Any]]) -> None:
        """
        Sets the __settings_current according to the settings given in settings_dict after validating them
        """
        # Validate settings_dict
        self.__complete(settings_dict)
        self.__validate(configparser.ConfigParser().read_dict(settings_dict))

        # set __settings_current
        self.__settings_current = settings_dict


    def __read_rcfile(self) -> Dict[str, Dict[str, Any]]:
        """
        Read from config data from file located at self.censorc_path 
        """
        rcdata: Dict = {}

        # read config file
        with open(self.censorc_path, "r") as file:
            parser = configparser.ConfigParser()
            parser.read_file(file)

        # Make sure that all the settings are included in the parser
        # If not, add defaults
        parser = self.__complete(parser)

        # validate parsed data
        self.__validate(parser)

        # convert parsed data to dict
        for section in parser.sections():
            rcdata[section] = {}
            for setting in parser[section]:
                rcdata[section][setting] = parser[section][setting]

        return rcdata


    def __complete(self, parser: configparser.ConfigParser) -> configparser.ConfigParser:
        """
        fill in missing settings with default values
        """
        for part, settings in self._settings_options.items():
            for setting, definition in settings.items():
                if setting not in parser[part]:
                    parser[part][setting] = definition["default"]

        return parser


    def __validate(self, parser: configparser.ConfigParser) -> None:
        """
        validate the type of each setting in the given parser
        also potentially validate if the setting is allowed by checking with CensoSettings._settings_options
        """
        # Create a mapping of data types to configparser methods
        mapping = {
            int: parser.getint,
            float: parser.getfloat,
            bool: parser.getboolean,
        }

        # go through each section and try to validate each setting's type
        for part in parser.sections():
            for setting_name in parser[part]:
                try:
                    mapping[self.get_type(part, setting_name)](part, setting_name)
                except KeyError:
                    # KeyError means that the type is not included in the mapping
                    # that means it's either a list or string
                    if self.get_type(part, setting_name) == list:
                        # try to convert to list
                        # SyntaxError not handled so it gets raised
                        ast.literal_eval(parser[part][setting_name])
                # ValueError not handled so it gets raised (happens if there is a type mismatch)

        # passed first step of validation, now check if settings are allowed for each part that should be run
        # (this works since for bools only the type needs to be checked to validate completely)
        for setting_name, setting_value in parser[part].items():
                setting_type = self.get_type(part, setting_name)
                # for strings check if string is within a list of allowed values
                if setting_type == str:
                    options = self._settings_options[part][setting_name]['options']
                    if setting_value not in options and len(options) > 0:
                        # Only check if there are options
                        # This is fatal so CENSO stops
                        raise ValueError(f"Value '{setting_value}' is not allowed for setting '{setting_name}' in part '{part}'.")
                # for numeric values check if value is within a range
                elif setting_type in (int, float):
                    interval = self._settings_options[part][setting_name]['range'] 
                    if not interval[0] <= setting_value <= interval[1]:
                        # This is fatal so CENSO stops
                        raise ValueError(f"Value '{setting_value}' is not allowed for setting '{setting_name}' in part '{part}'.")
                # setting_type is None if setting does not exist in _settings_options
                elif setting_type is None:
                    raise ValueError(f"Unknown setting type for setting '{setting_name}' in part '{part}'")


    def write_rcfile(self, path: str) -> None:
        """
        write new configuration file with default settings into file at 'path' 
        """
        # what to do if there is an existing configuration file
        external_paths = None
        if os.path.isfile(path):
            print(
                f"An existing configuration file has been found at {path}. " 
            )

            # Read program paths from the existing configuration file 
            print("Reading program paths from existing configuration file ...")
            external_paths = self.__read_program_paths(path)

        print(f"Writing new configuration file to {path} (old file will be overwritten) ...")
        with open(path, "w", newline=None) as rcfile:
            parser = configparser.ConfigParser()
            parser.read_dict({part: {setting: settings[setting]['default'] for setting in settings} for part, settings in self._settings_options.items()})

            # Try using the paths from the old file
            if not external_paths is None:
                parser['paths'] = external_paths
            else:
                # TODO - Try reading the paths from environment variables
                pass

            parser.write(rcfile)

        print(
            f"\nA new configuration file was written into {path}.\n"
            "Put the file either in your home directory or in the current working directory.\n"
            "You should adjust the settings to your needs and set the program paths.\n"
            "Right now only the default settings are used.\n"
        )

        if ".censorc" not in path:
            print(
                f"Additionally make sure that the file name is '{CENSORCNAME}'.\n"
                f"Currently it is '{os.path.split(path)[-1]}'.\n"
            )


    def __find_rcfile(self, inprcpath=None) -> Union[str, None]:
        """
        check for existing censorc
        looks for custom path and standard paths: cwd and $home dir
        if there is a configuration file in cwd and $home, it prioritizes no the one in cwd
        """

        censorc_name = CENSORCNAME

        rcpath = None
        # check for .censorc in $home if no path is given 
        if not inprcpath is None:
            if os.path.isfile(os.path.join(os.path.expanduser("~"), censorc_name)):
                rcpath = os.path.join(os.path.expanduser("~"), censorc_name)
        elif inprcpath and os.path.isfile(inprcpath):
            # if path is given and file exists, take it
            rcpath = inprcpath
        elif inprcpath and not os.path.isfile(inprcpath):
            raise FileNotFoundError(f"Configuration file {inprcpath} not found!")

        return rcpath


    def __read_program_paths(self, path: str) -> Union[Dict[str, str], None]:
        """
        Read program paths from the configuration file at 'path'
        """
        with open(path, "r") as inp:
            parser = configparser.ConfigParser()
            parser.read_file(inp)

        try:
            return parser["paths"]
        except KeyError:
            print(f"WARNING: No paths found in {path}")
            return None


    def configure(self, args: Namespace) -> None:
        """
        Overwrite the settings of CENSO with the given arguments
        """
        # TODO
        pass