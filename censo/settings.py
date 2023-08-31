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
                            WARNLEN, GRIDOPTIONS, GFNOPTIONS, SettingsDict, __version__)
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


    @property
    def functionals(self) -> Dict[str, Dict]:
        return self.__dfa_dict["functionals"]


    @property
    def composite_bs(self) -> set:
        """
        return all composite method basis sets
        """
        return set([v for v in self.__dfa_dict["composite_method_basis"].values()])


    @property
    def composites(self) -> set:
        """
        return all composite method dfas dict entries
        """
        functionals = self.__dfa_dict["functionals"]
        return set(filter(lambda x: "composite" in functionals[x]["type"], functionals))


    @property
    def ggas(self) -> set:
        """
        return all (m)gga dfas dict entries
        """
        functionals = self.__dfa_dict["functionals"]
        return set(filter(lambda x: "gga" in functionals[x]["type"] and not "composite" in functionals[x]["type"], functionals))


    @property
    def hybrids(self) -> set:
        """
        return all hybrid dfas dict entries
        """
        functionals = self.__dfa_dict["functionals"]
        return set(filter(lambda x: "hybrid" in functionals[x]["type"] and not "composite" in functionals[x]["type"], functionals))


    @property
    def doublehs(self) -> set:
        """
        return all double hybrid dfas dict entries
        """
        functionals = self.__dfa_dict["functionals"]
        return set(filter(lambda x: "double" in functionals[x]["type"] and not "composite" in functionals[x]["type"], functionals))



@dataclass
class Setting:
    type: type
    part: str
    name: str
    value: Union[int, float, str, list, bool]
    options: Union[Tuple, None]
    default: Union[int, float, str, list, bool]
    
    def __str__(self):
        return f"{self.name}: {self.value}"
    
    # TODO - add setter that automatically checks what value you try to assign to the setting
    # regarding type and options
    

class SettingsTuple(tuple):
    def __init__(self, *args, **kwargs):
        # the try-except block makes sure that the constructor only accepts 
        # sequences that contain only 'Setting' objects
        # (tuple constructor requires a sequence anyways)
        try:
            for arg in args:
                # if condition just to make sure only sequences are checked
                if hasattr(arg, "__iter__"):
                    for item in arg:
                        assert type(item) == Setting
        except AssertionError:
            raise TypeError("Cannot create a 'SettingsTuple' instance because the input sequence is not entirely made up of 'Setting' objects.")
        
        super().__init__()
    

    def __add__(self, other):
        """
        overload the '+' operator for tuple to get the correct return type
        can also accept Setting objects
        """
        try:
            return SettingsTuple([x for x in self] + [y for y in other])
        except TypeError:
            if type(other) == Setting:
                return SettingsTuple([x for x in self] + [other])
            else:
                raise(TypeError("'+' operation for SettingsTuple only supported for iterables containing exclusively Setting objects or single Setting objects."))

    
    def byname(self, name: str) -> Union[Setting, None]:
        """
        returns the first 'Setting' object in itself that has the given name
        if no matching one is found 'None' is returned
        """
        # TODO - make this more like 'bypart'
        for item in self:
            if item.name == name:
                return item
    
    
    def bypart(self, part: str) -> List[Setting]:
        """
        returns a list of all 'Setting' objects that belong to a given part
        if no matches are found an empty list is returned
        """
        
        matches = []
        
        for item in self:
            if item.part == part:
                matches.append(item)

        return matches
        
    
    def get_setting(self, type_t: type, part: str, name: str) -> Union[Setting, None]:
        """
        returns exactly the 'Setting' object defined by the args
        if it is not found 'None' is returned
        """
        for item in self:
            if item.type == type_t and item.part == part and item.name == name:
                return item
        
        return None
    
    
    def get_settings(self, types: List[type], parts: List[str], names: List[str]) -> List[Setting]:
        """
        returns a list of 'Setting' objects that match the args
        """
        matches = []
        for item in self:
            if item.type in types and item.part in parts and item.name in names:
                matches.append(item)
                
        return matches
    

class CensoSettings:
    """
    All options are saved here.
    Manages internal settings.
    CensoSettings is implemented as thread-safe singleton
    (there should only be one instance for the run)
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

    # load up all resources to set value options in settings_options
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
    
    # default settings, defines sections and default settings
    # they are used if settings are missing from the configuration file (TODO - provide warning!)
    __defaults = {
        'general': {
            # strs
            'cosmorsparam': 'automatic',
            'sm_rrho': 'alpb',
            'solvent': 'h2o',
            # ints
            'charge': 0,
            'maxprocs': 1,
            'omp': 1,
            'unpaired': 0,
            # floats
            'imagethr': -100.0,
            'scale': 1.0,
            'sthr': 0.0,
            'temperature': 298.15,
            # lists
            'trange': [273.15, 378.15, 5],
            # bools
            'balance': True,
            'bhess': True,
            'check': True,
            'consider_sym': True,
            'evaluate_rrho': True,
            'gas-phase': False,
            'multitemp': True,
            'nmrmode': False,
            'progress': False,
            'rmsdbias': False,
            'vapor_pressure': False,
        },
        'nmr': {
            # strs
            'basis_j': 'def2-TZVP',
            'basis_s': 'def2-TZVP',
            'c_ref': 'TMS',
            'f_ref': 'CFCl3',
            'func_j': 'pbe0-d4',
            'func_s': 'pbe0-d4',
            'h_ref': 'TMS',
            'p_ref': 'TMP',
            'prog4_j': 'tm',
            'prog4_s': 'tm',
            'si_ref': 'TMS',
            'sm4_j': 'default',
            'sm4_s': 'default',
            # floats
            'resonance_frequency': 300.0,
            # bools
            'c_active': True,
            'couplings': True,
            'f_active': False,
            'h_active': True,
            'p_active': False,
            'run': False,
            'shieldings': True,
            'si_active': False,
        },
        'optimization': {
            # strs
            'basis': 'automatic',
            'func': 'r2scan-3c',
            'gfnv': 'gfn2',
            'grid': 'high',
            'optlevel2': 'automatic',
            'prog': 'orca',
            'prog2opt': 'prog',
            'sm': 'smd',
            'smgsolv': 'smd_gsolv',
            # ints
            'optcycles': 8,
            'radsize': 10,
            # floats
            'hlow': 0.01,
            'optimization_P_threshold': 99.0,
            'spearmanthr': 0.0,
            'threshold': 2.5,
            # bools
            'crestcheck': False,
            'gcp': True,
            'opt_spearman': True,
            'run': True,
        },
        'optrot': {
            # strs
            'basis': 'def2-SVPD',
            'func': 'pbe-d4',
            'func_or_scf': 'r2scan-3c',
            'prog': 'orca',
            # lists
            'freq_or': [598.0],
            # bools
            'run': False,
        },
        'prescreening': {
            # strs
            'basis': 'def2-SV(P)',
            'func': 'pbe-d4',
            'gfnv': 'gfn2',
            'grid': 'low',
            'prog': 'orca',
            # floats
            'threshold': 4.0,
            # bools
            'gcp': True,
            'run': True,
        },
        'refinement': {
            # strs
            'basis': 'def2-TZVPP',
            'func': 'wb97x-v',
            'gfnv': 'gfn2',
            'grid': 'high+',
            'prog': 'orca',
            'smgsolv': 'smd_gsolv',
            # floats
            'threshold': 99.0,
            # bools
            'gcp': True,
            'run': False,
        },
        'screening': {
            # strs
            'basis': 'automatic',
            'func': 'r2scan-3c',
            'gfnv': 'gfn2',
            'grid': 'low+',
            'prog': 'orca',
            'smgsolv': 'smd_gsolv',
            # floats
            'threshold': 3.5,
            # bools
            'gcp': True,
            'run': True,
        },
        'uvvis': {
            # ints
            'nroots': 20,
            # floats
            'sigma': 0.1,
            # bools
            'run': False,
        },
    },
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
    settings_options = {
        "general": {
            "charge": {
                "default": 0,
                "range": [
                    -10,
                    10
                ]
            },
            "unpaired": {
                "default": 0,
                "range": [
                    0,
                    14
                ]
            },
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

    @property
    def settings_options(cls) -> SettingsTuple:
        """
        returns settings_options converted to SettingsTuple
        """
        default = SettingsTuple()

        # loop over all known settings (hardcoded in settings_options)
        for type_t, parts in cls.settings_options.items():
            for part, settings in parts.items():
                if settings:
                    # if settings exist create 'Setting' objects with default settings
                    for setting, definition in settings.items():
                        options = definition.get("options", None)
                        
                        default += Setting(
                            type_t, 
                            part, 
                            setting, 
                            definition["default"], 
                            options, 
                            definition["default"]
                        )
        
        return default
                
                
    @classmethod
    def get_type(cls, name: str) -> Union[type, None]:
        """
        returns the type of a given setting
        """
        for type_t, parts in cls.settings_options.items():
            for settings in parts.values():
                if settings:
                    if name in settings.keys():
                        return type_t

        return None


    @classmethod
    def byname(cls, name: str) -> Union[Dict, None]:
        """
        returns the definition of a given setting
        """
        for type_t, parts in cls.settings_options.items():
            for settings in parts.values():
                if settings:
                    if name in settings.keys():
                        return settings[name]
                    
        return None

        
    def __init__(self):
        """
        setup with default settings, all other attributes blank
        """   
        
        # stores all settings, grouped by section (e.g. "general", "prescreening", ...)
        # every string is in lower case
        self.__settings_current: Dict[str, Dict[str, Any]]
        
        # absolute path to configuration file, try to find .censorc on construction
        self.censorc_path: str = self.__find_rcfile()

        # assert that a config file is found before trying to read it
        try:
            assert self.censorc_path is not None
        except AssertionError:
            # TODO - error handling
            raise Exception("Configuration file not found")

        # read config file
        self.__settings_current = self.__read_config(self.censorc_path)

        # TODO - try to read paths from environment vars?
        self.external_paths: Dict[str, str] = {}
        self.external_paths["orcapath"] = ""
        self.external_paths["orcaversion"] = "" # TODO - maybe remove this from here
        self.external_paths["xtbpath"] = ""
        self.external_paths["crestpath"] = ""
        self.external_paths["cosmorssetup"] = ""
        self.external_paths["dbpath"] = ""
        self.external_paths["cosmothermversion"] = "" # TODO - maybe remove this from here
        self.external_paths["mpshiftpath"] = ""
        self.external_paths["escfpath"] = ""

        # read program paths from config file
        self.__read_program_paths()
    
    
    def print_paths(self) -> None:
        """
        print out paths of all external qm programs
        """
        lines = []
        
        lines.append("\n" + "".ljust(PLENGTH, "-") + "\n")
        lines.append("PATHS of external QM programs".center(PLENGTH, " ") + "\n")
        lines.append("".ljust(PLENGTH, "-") + "\n")
        
        for program, path in self.external_paths.items():
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
    def settings_current(self, args: Namespace):
        """
        iterate over all settings and set according to rcfile first, then cml args, 
        implicitily calls check_logic
        """
        rcdata = self.read_config(self.censorc_path)
        for setting in self.__settings_current:
            # update settings with rcfile
            try:
                # everything converted to lower case since there is no need for case sensitivity
                setting.value = rcdata[setting.type][setting.part][setting.name].lower()
                rcdata[setting.type][setting.part].pop(setting.name)
            except KeyError:
                print(f"Could not find {setting} in rcdata.")
            
            # overwrite settings with cml arguments (once again converting everything to lower case)
            # TODO - coverage between all actual settings and cml options?
            try:
                if getattr(args, setting.name) is not None:
                    setting.value = getattr(args, setting.name).lower()
                else:
                    print(f"{setting.name} not set via cml args.")
            # TODO - error handling
            except Exception:
                print(LogicWarning(
                    f"{setting.name}",
                    f"{setting.name} no attribute of cml args NameSpace object.",
                    f"{setting.name} might be missing in the definition of cml args in 'inputhandling.py'.",
                ))
              
        print(f"Remaining settings:\n{rcdata}\n") # TODO - where to print out?
        
        warnings = self.check_logic()
        for warning in warnings:
            print(warning)
            
        if any([warning.fatal for warning in warnings]):
            print("\nFatal error in user input while setting up run settings. Going to exit!")
            sys.exit(1)
            

    def __read_config(self) -> SettingsDict:
        """
        Read from config data from file (here .censorc),
        every part has to be in rcfile
        """
        rcdata: Dict = {}

        # read config file
        with open(self.censorc_path, "r") as file:
            parser = configparser.ConfigParser(...)
            parser.read_file(file)

        # validate parsed data
        self.__validate(parser)

        # convert parsed data to dict
        for section in parser.sections():
            rcdata[section] = {}
            for setting in parser[section]:
                rcdata[section][setting] = parser[section][setting]

        return rcdata

    
    def __validate(self, parser: configparser.ConfigParser) -> None:
        """
        validate the type of each setting in the given parser
        also potentially validate if the setting is allowed by checking with CensoSettings.settings_options
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
                    mapping[self.get_type(setting_name)](part, setting_name)
                except KeyError:
                    # KeyError means that the type is not included in the mapping
                    # that means it's either a list or string
                    if self.get_type(setting_name) == list:
                        # try to convert to list
                        # SyntaxError not handled so it gets raised
                        ast.literal_eval(parser[part][setting_name])
                # ValueError not handled so it gets raised (happens if there is a type mismatch)

        # passed first step of validation, now check if settings are allowed for each part that should be run
        # (this works since for bools only the type needs to be checked to validate completely)
        for setting_name, setting_value in parser[part].items():
                setting_type = self.get_type(setting_name)
                # for strings check if string is within a list of allowed values
                if setting_type == str:
                    options = self.settings_options[part][setting_name]['options']
                    if setting_value not in options:
                        # This is fatal so CENSO stops
                        raise ValueError(f"Value '{setting_value}' is not allowed for setting '{setting_name}' in part '{part}'.")
                # for numeric values check if value is within a range
                elif setting_type in (int, float):
                    interval = self.settings_options[part][setting_name]['range'] 
                    if not interval[0] <= setting_value <= interval[1]:
                        # This is fatal so CENSO stops
                        raise ValueError(f"Value '{setting_value}' is not allowed for setting '{setting_name}' in part '{part}'.")


    def write_config(self, args: Namespace, cwd: str) -> str:
        """
        write new configuration file with default settings into 
        either the local 'censorc_new' or ~/.censorc.
        returns the path of the new configuration file
        """
        # what to do if there is an existing configuration file
        if hasattr(self, "censorc_path") and self.censorc_path != "":
            print(
                f"An existing configuration file ({os.path.split(self.censorc_path)[-1]})" 
                f"has been found in {self.censorc_path}"
                f"Do you want to copy existing program path information to the "
                f"new remote configuration file?"
            )
            
            user_input = ""
            while user_input.strip().lower() not in ["yes", "y", "no", "n"]:
                print("Please type 'yes/y' or 'no/n':")
                user_input = input()
            
            if user_input.strip().lower() in ("y", "yes"):
                self.read_program_paths()

        # write new censorc
        if args.inprcpath:
            path = args.inprcpath
        # (path is never unbound)
        else:
            path = os.path.join(cwd, "censorc_new")
            
        print(f"Writing configuration file to {path} ...")
        with open(path, "w", newline=None) as outdata:
            # TODO - purpose?
            outdata.write(f"$CENSO global configuration file: {os.path.split(path)[-1]}\n") 
            outdata.write(f"$VERSION: {__version__} \n")
            outdata.write("\n")
            
            # if ALL program paths are set, write the paths to the new rcfile
            if all([s != "" for s in self.external_paths.values()]):
                # write stored program paths to file
                outdata.writelines([
                    f"ORCA: {self.external_paths['orcapath']}\n",
                    f"ORCA version: {self.external_paths['orcaversion']}\n",
                    f"GFN-xTB: {self.external_paths['xtbpath']}\n",
                    f"CREST: {self.external_paths['crestpath']}\n",
                    f"mpshift: {self.external_paths['mpshiftpath']}\n",
                    f"escf: {self.external_paths['escfpath']}\n",
                    "\n",
                    "#COSMO-RS\n",
                    f"{self.external_paths['cosmorssetup']}\n",
                ])
            else:
                # TODO - write some other default (e.g. "") instead of paths
                outdata.writelines([
                    "ORCA: /path/including/binary/orca-binary\n",
                    "ORCA version: 4.2.1\n",
                    "GFN-xTB: /path/including/binary/xtb-binary\n",
                    "CREST: /path/including/binary/crest-binary\n",
                    "mpshift: /path/including/binary/mpshift-binary\n",
                    "escf: /path/including/binary/escf-binary\n",
                    "\n",
                    "#COSMO-RS\n",
                ])
                outdata.write(
                    "ctd = BP_TZVP_C30_1601.ctd cdir = "
                    '"/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES" ldir = '
                    '"/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES"\n'
                )
            outdata.write("$ENDPROGRAMS\n\n")
            #outdata.write("$CRE SORTING SETTINGS:\n")

            ### loop to write settings to rcfile
            # create headers for sections in censorc and write them
            # get defaults for settings and write them
            for part in ("general",) + PARTS:
                outdata.write(f"\n${part.upper()} SETTINGS\n")
                
                # iterate over the values for all types (dicts mapping partnames to settings)
                for parts in self.settings_options.values():
                    # get only the settings for the current part
                    settings = parts[part]
                    
                    # check if settings are defined for this part
                    if settings:
                        # write defaults to rcfile
                        for setting, definition in settings.items():
                            outdata.write(f"{setting}: {definition['default']}\n")
                
            outdata.write("\n$END CENSORC\n")

        print(
            "\nA new configuration file was written into the current directory file:\n"
            "censorc_new\n"
            "You have to adjust the settings to your needs"
            " and it is mandatory to correctly set the program paths!\n"
            "Additionally move the file to the correct filename: '.censorc'\n"
            "and place it either in your /home/$USER/ or current directory.\n"
        )

        return path


    def __find_rcfile(self, cwd: str, inprcpath=None) -> Union[str, None]:
        """
        check for existing censorc
        looks for custom path and standard paths: cwd and $home dir
        if there is a configuration file in cwd and $home, it prioritizes no the one in cwd
        """

        censorc_name = ".censorc"

        # mapping the paths defined above to True/False, 
        # depending if file exists or not
        tmp = [
            os.path.join(cwd, censorc_name),
            os.path.join(os.path.expanduser("~"), censorc_name)
        ]
        check = {
            os.path.isfile(tmp[0]): tmp[0],
            os.path.isfile(tmp[1]): tmp[1],
        }

        rcpath = None

        # TODO - probably doesn't catch all cases, needs testing
        # check for .censorc in standard locations if no path is given
        if not inprcpath:
            if all(list(check.keys())) and not tmp[0] == tmp[1]:
                # if file exists in cwd, prioritize it
                # TODO - give a warning
                rcpath = tmp[0]
            elif any(list(check.keys())):
                # take the one file found
                rcpath = check[True]
        elif inprcpath and os.path.isfile(inprcpath):
            # if path is given and file exists, take it
            rcpath = inprcpath

        return rcpath


    def __read_program_paths(self):
        """
        Set absolute paths of external programs employed in censo
        Read from the configuration file .censorc
        """
        # TODO - make this nicer?
        # TODO - fix this with readrcfile decorator
        with open(self.censorc_path, "r") as inp:
            for line in inp.readlines():
                if "ctd =" in line:
                    try:
                        self.external_paths["cosmorssetup"] = str(line.rstrip(os.linesep))
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read settings for COSMO-RS from .censorc!"
                        )
                    try:
                        normal = "DATABASE-COSMO/BP-TZVP-COSMO"
                        fine = "DATABASE-COSMO/BP-TZVPD-FINE"
                        tmp_path = self.external_paths["cosmorssetup"].split()[5].strip('"')
                        if "OLDPARAM" in tmp_path:
                            tmp_path = os.path.split(tmp_path)[0]
                        tmp_path = os.path.split(tmp_path)[0]
                        self.external_paths["dbpath"] = tmp_path
                        self.external_paths["dbpath_fine"] = os.path.join(tmp_path, fine)
                        self.external_paths["dbpath_normal"] = os.path.join(
                            tmp_path, normal
                        )
                    except Exception as e:
                        print(e)
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read settings for COSMO-RS from "
                            f".censorc!\n{'':{WARNLEN}}Most probably there is a user "
                            "input error."
                        )
                if "ORCA:" in line:
                    try:
                        self.external_paths["orcapath"] = str(line.split()[1])
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read path for ORCA from .censorc!."
                        )
                if "ORCA version:" in line:
                    try:
                        tmp = line.split()[2]
                        tmp = tmp.split(".")
                        tmp.insert(1, ".")
                        tmp = "".join(tmp)
                        self.external_paths["orcaversion"] = tmp
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read ORCA version from .censorc!"
                        )
                if "GFN-xTB:" in line:
                    try:
                        self.external_paths["xtbpath"] = str(line.split()[1])
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read path for GFNn-xTB from .censorc!"
                        )
                        
                        xtbpath = shutil.which("xtb")
                        if not xtbpath:
                            raise Exception # TODO
                            
                        self.external_paths.update({"xtbpath": xtbpath})
                        print(
                            f"{'':{WARNLEN}}Going to use {self.external_paths['xtbpath']} instead."
                        )
                            
                if "CREST:" in line:
                    try:
                        self.external_paths["crestpath"] = str(line.split()[1])
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read path for CREST from .censorc!"
                        )
                        if shutil.which("crest") is not None:
                            crestpath = shutil.which("crest")
                            if not crestpath:
                                raise Exception # TODO
                            
                            self.external_paths.update({"crestpath": crestpath})
                            print(
                                f"{'':{WARNLEN}}Going to use {self.external_paths['crestpath']} instead."
                            )
                if "mpshift:" in line:
                    try:
                        self.external_paths["mpshiftpath"] = str(line.split()[1])
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read path for mpshift from .censorc!"
                        )
                if "escf:" in line:
                    try:
                        self.external_paths["escfpath"] = str(line.split()[1])
                    except Exception:
                        print(
                            f"{'WARNING:':{WARNLEN}}Could not read path for escf from .censorc!"
                        )
                if "$ENDPROGRAMS" in line:
                    break


    def check_logic(self) -> List[LogicWarning]:
        """
        Checks internal settings for impossible setting-combinations
        also checking if calculations are possible with the requested qm_codes.
        """
        # TODO - complete it for all parts
        
        print("\nCHECKING SETTINGS ...\n")
        
        # do not affect execution, although might lead to errors/unexpected behaviour down the line
        warnings = []
        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle prog_rrho
        # TODO - reimplement when there are more options available
        """ if not self.prog_rrho:
            self.save_errors.append(
                f"{'WARNING:':{WARNLEN}}Thermostatistical contribution to "
                "free energy will not be calculated, since prog_rrho ist set to 'off'!"
            )
            if shutil.which("thermo") is not None:
                # need thermo for reading thermostatistical contribution
                self.prog_rrho = "tm"
            else:
                self.prog_rrho = "xtb"
                self.save_errors.append(
                    f"{'WARNING:':{WARNLEN}}Currently are only GFNn-xTB "
                    "hessians possible and no TM hessians"
                ) """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not self.settings_current.get_setting(bool, "general", "gas-phase"):
            with open(os.path.join(ASSETS_PATH, "solvents.json"), "r") as file:
                solvents = json.load(file)
                
            with open(os.path.join(ASSETS_PATH, "solvents_dc.json"), "r") as file:
                solvents_dc = json.load(file)
            
            # FIXME - ???
            """ if self.settings_current(settings="vapor_pressure"):
                errors.append(LogicError(
                        "vapor_pressure",
                        "The vapor_pressure flag only affect settings for COSMO-RS.",
                        "Information on solvents with properties similar to the input molecule must be provided for other solvent models!"
                    )
                ) """
            
            # check availability of solvent model for given program in parts 1-3
            # check for solvent availability for given solvent model
            # TODO - what about "replacements" for solvents
            # TODO - add flexibility concerning prescreening, optrot, uvvis
            solvent = getattr(
                self.settings_current.get_setting(str, "general", "solvent"), 
                "value", 
                None
            )
            multitemp = getattr(
                self.settings_current.get_setting(bool, "general", "multitemp"),
                "value",
                None
            )
                
            if not (solvent is None or multitemp is None):
                for part in PARTS[1:4]:
                    # get 'prog', 'sm' and 'smgsolv' settings for the respective part
                    settings = SettingsTuple(
                        self.settings_current.get_settings(
                            [str], 
                            [part], 
                            ["prog", "sm", "smgsolv"]
                        )
                    )
                    
                    for model in ["sm", "smgsolv"]:

                        prog = getattr(settings.byname("prog"), "value", None)
                        model = getattr(settings.byname(model), "value", None)
                            
                        if not (prog is None or model is None):
                            # check solvent availability in solvent model (warning)
                            if solvent not in solvents[model]:
                                warnings.append(LogicWarning(
                                        "solvent",
                                        f"Solvent {solvent} was not found to be available for solvent model {model} in part {part}.",
                                        "",
                                   )
                                )

                            # check sm/smgsolv availability in program (warning)
                            if (
                                model not in SOLV_MODS.get("prog", tuple([])) + GSOLV_MODS.get("prog", tuple([]))
                            ):
                                warnings.append(LogicWarning(
                                        "sm/smgsolv",
                                        f"Solvent model {model} was not found to be available with given program {prog} in part {part}.",
                                        "",
                                    )
                                )
                    # TODO - check for dielectric constant for cosmo ->  dc only needed for dcosmors?
                        else:
                            warnings.append(LogicWarning(
                                    f"prog/{model}",
                                    f"One of the settings is not defined for part {part}.",
                                    "Cannot check for availability.",
                                )
                            )
                    
                    prog = getattr(settings.byname("prog"), "value", None)
                    
                    if not prog is None:
                        # there is no multitemp support for any orca smgsolv (error)
                        if (
                            prog == "orca"
                            and multitemp
                        ):
                            warnings.append(LogicWarning(
                                    "multitemp",
                                    "There is no multitemp support for any solvent models in ORCA.",
                                    "Disable multitemp.",
                                    fatal=True              
                                )
                            )
                    else:
                        raise LogicError(
                            "prog",
                            "Setting is not defined for part {part}.",
                            "Check initialization of run settings."
                        )
            else:
                raise LogicError(
                    "multitemp/solvent",
                    "Settings are not defined.",
                    "Check initialization of run settings."                  
                )
            
        stroptions = CensoSettings.settings_options[str]

        for part in PARTS:
            # get settings for part
            # look for func/basis
            # check if func/basis values are allowed in settings_options
            settings = SettingsTuple(
                self.settings_current.get_settings(
                    [str, bool], 
                    [part], 
                    ["run", "prog", "func", "basis"]
                )
            )
            
            run = getattr(settings.byname("run"), "value", None)
            prog = getattr(settings.byname("prog"), "value", None)
            func = getattr(settings.byname("func"), "value", None)
            basis = getattr(settings.byname("basis"), "value", None)
                
            if all(x is not None for x in [run, prog, func, basis]):
                # iterate through settings and check for func
                # TODO - handle composite method bases
                if run and not stroptions[part] is None:
                    if not func in stroptions[part]["func"]["options"]:
                        warnings.append(LogicWarning(
                                "func",
                                f"The functional {func} was not found to be available for any software package in part {part}.",
                                ""
                            )
                        )
                    else:
                        # FIXME - doesn't work properly yet
                        # check for DFA availability in prog (warning)
                        if not func in CensoSettings.dfa_settings.find_func(part, prog):
                            warnings.append(LogicWarning(
                                    "func",
                                    f"The functional {func} was not found to be available for the software package {prog} in part {part}.",
                                    ""
                                )
                            )
                        # extra check for r2scan-3c
                        elif (
                            func == "r2scan-3c" 
                            and prog == "orca"
                        ):
                            try:
                                if int(self.external_paths["orcaversion"].split(".")[0]) < 5:
                                    warnings.append(LogicWarning(
                                            "func",
                                            f"The functional r2scan-3c is only available from ORCA version 5 in part {part}.",
                                            "Choose a different functional or use a newer version of ORCA.",
                                            fatal=True    
                                        )
                                    )
                            except Exception:
                                warnings.append(LogicWarning(
                                        "func",
                                        f"The availability for {func} could not be checked because the version of Orca could not be determined.",
                                        ""
                                    )
                                )
                        # extra check for dsd-blyp
                        elif (
                            func in ("dsd-blyp", "dsd-blyp-d3")
                            and basis != "def2-TZVPP"
                        ):
                            warnings.append(LogicWarning(
                                    "func/basis",
                                    f"The functional {func} is only available with the def2-TZVPP basis set in part {part}.",
                                    "Change the basis set to def2-TZVPP or choose a different functional.",
                                    fatal=True  
                                )
                            )
                            
                        # check for basis set availability in prog (warning)
                        if not basis in stroptions[part]["basis"]["options"]:
                            warnings.append(LogicWarning(
                                    "basis",
                                    f"The basis set {basis} was not found to be available for any software package in part {part}.",
                                    ""    
                                )
                            )
            else:
                warnings.append(LogicWarning(
                        "run/prog/func/basis",
                        f"Settings at least partially not defined for part {part}",
                        "Cannot check for availability."
                    )
                )

        return warnings