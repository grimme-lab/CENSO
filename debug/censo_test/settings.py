import json
import os
import shutil
import sys
from argparse import Namespace
from dataclasses import dataclass
from functools import reduce
from multiprocessing import Lock
from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Type, Union

from censo_test.cfg import (ASSETS_PATH, GSOLV_MODS, PARTS, PROGS, SOLV_MODS,
                            WARNLEN, SettingsDict, __version__)
from censo_test.errorswarnings import LogicError, LogicWarning


class DfaSettings:
    def __init__(self, obj: Dict[str, Dict[str, Dict]]):
        self.dfa_dict = obj

    
    def find_func(self, part: str, prog=None):
        """
        return all functionals available for a certain part and (optionally) program
        """
        tmp = []
        for k, v in self.dfa_dict["functionals"].items():
            if part in v["part"]:
                if prog:
                    if v[prog] != "":
                        tmp.append(k)
                else:
                    tmp.append(k)
        
        return tmp


    def composite_bs(self) -> set:
        """
        return all composite method basis sets
        """
        return set([v for v in self.dfa_dict["composite_method_basis"].values()])


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
            raise TypeError("Cannot create a 'SettingsTuple' instance because the input sequence is not entirely made up of 'Setting' objects")
        
        super().__init__()
    
    
    def byname(self, name: str) -> Union[Setting, None]:
        """
        returns the first 'Setting' object in itself that has the given name
        if no matching one is found 'None' is returned
        """
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
    
    
    def get_settings(self, types: list[type], parts: list[str], names: list[str]) -> List[Setting]:
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

    # TODO - defaults for uvvis (func: w-B97X-D3)
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
    settings_options = MappingProxyType({
        int: MappingProxyType({
            "general": MappingProxyType({
                "charge": {"default": 0, "options": (-10, 10)}, # TODO - (re)move?
                "unpaired": {"default": 0, "options": (0, 14)}, # TODO - (re)move?
                "maxprocs": {"default": 1, "options": (1, 1024)}, # number of processes
                "omp": {"default": 1, "options": (1, 256)}, # number of cores per process
            }),
            "prescreening": None,
            "screening": None,
            "optimization": MappingProxyType({
                "optcycles": {"default": 8, "options": (1, 100)}, # TODO - which value as min?
                "radsize": {"default": 10, "options": (1, 100)}, # TODO - same
            }),
            "refinement": None,
            "nmr": None,
            "optrot": None,
            "uvvis": MappingProxyType({
                "nroots": {"default": 20, "options": (1, 100)}, # TODO - set dynamically
            }),
        }),
        float: MappingProxyType({
            "general": MappingProxyType({
                "imagethr": {"default": -100.0, "options": (-300.0, 0.0)}, # TODO - threshold for imaginary frequencies
                "sthr": {"default": 0.0, "options": (0.0, 100.0)}, # TODO - what is this?
                "scale": {"default": 1.0, "options": (0.0, 1.0)}, # TODO - what is this?
                "temperature": {"default": 298.15, "options": (0.1, 2000.0)}, # TODO - 0.0 still works?
            }),
            "prescreening": MappingProxyType({
                "threshold": {"default": 4.0, "options": (1.0, 10.0)}, # TODO - which value as min?
            }),
            "screening": MappingProxyType({
                "screening_threshold": {"default": 3.5, "options": (0.75, 7.5)},
            }),
            "optimization": MappingProxyType({
                "opt_limit": {"default": 2.5, "options": (0.5, 5)}, # TODO - rename?
                "hlow": {"default": 0.01, "options": (0.01, 1.0)}, # TODO
                "optimization_P_threshold": {"default": 99.0, "options": (1.0, 10.0)}, # TODO
                "spearmanthr": {"default": 0.0, "options": (0.0, 10.0)},
            }),
            "refinement": MappingProxyType({
                "refinement_threshold": {"default": 99.0, "min": (1.0, 10.0)}, # TODO
            }),
            "nmr": MappingProxyType({
                "resonance_frequency": {"default": 300.0, "options": (150.0, 1000.0)}, # TODO
            }),
            "optrot": None,
            "uvvis": MappingProxyType({
                "sigma": {"default": 0.1, "options": (0.1, 1.0)},
            }),
        }),
        str: MappingProxyType({
            "general": MappingProxyType({
                "solvent": {"default": "gas", "options": ("gas",) + tuple([k for k in solvents_db.keys()])},
                #"prog_rrho": {"default": "xtb", "options": ("xtb")}, # TODO - keep this here?
                "sm_rrho": {"default": "alpb", "options": ("alpb", "gbsa")}, # TODO - same
                "cosmorsparam": {"default": "automatic", "options": tuple([k for k in cosmors_param.keys()])},
            }),
            "prescreening": MappingProxyType({
                "func": {"default": "b97-d3(0)", "options": tuple(dfa_settings.find_func("prescreening"))},
                "basis": {"default": "def2-SV(P)", "options": ("automatic",) + tuple(dfa_settings.composite_bs()) + ("def2-SV(P)", "def2-TZVP")},
                "prog": {"default": "tm", "options": PROGS},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
            }),
            "screening": MappingProxyType({
                "func": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func1"))},
                "basis": {"default": "automatic", "options": ("automatic",) + basis_sets},
                "prog": {"default": "tm", "options": PROGS},
                "smgsolv": {"default": "cosmors", "options": gsolv_mods},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
            }),
            "optimization": MappingProxyType({
                "func": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func2"))},
                "basis": {"default": "automatic", "options": ("automatic",) + basis_sets},
                "prog": {"default": "tm", "options": PROGS},
                "prog2opt": {"default": "prog", "options": PROGS}, # TODO - ??? prog2 ??? # FIXME
                "sm": {"default": "default", "options": solv_mods}, # FIXME
                "smgsolv": {"default": "cosmors", "options": gsolv_mods},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
                "optlevel2": {"default": "automatic", "options": ("crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme", "automatic")}, # TODO - what does this mean?
            }),
            "refinement": MappingProxyType({
                "prog": {"default": "tm", "options": PROGS},
                "func": {"default": "pw6b95-d4", "options": tuple(dfa_settings.find_func("func3"))},
                "basis": {"default": "def2-TZVPD", "options": basis_sets},
                "smgsolv": {"default": "cosmors", "options": gsolv_mods},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
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
                "balance": {"default": False},
                "vapor_pressure": {"default": False},
                "nmrmode": {"default": False},
                "gas-phase": {"default": False},
            }),
            "prescreening": MappingProxyType({
                "run": {"default": True},
            }),
            "screening": MappingProxyType({
                "run": {"default": True},
            }),
            "optimization": MappingProxyType({
                "run": {"default": True},
                # "ancopt": {"default": True},
                "opt_spearman": {"default": True},
                "crestcheck": {"default": False},
            }),
            "refinement": MappingProxyType({
                "run": {"default": False},
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
    })
                
                
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
        
        # stores all settings in a SettingsTuple (extending built-in tuple)
        self._settings_current: SettingsTuple = self.default_settings()
        
        # absolute path to configuration file
        self.censorc_path: str
        
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
         

    @property
    def settings_current(self) -> SettingsTuple:
        """
        return the complete _settings_current
        """
        return self._settings_current

    
    @settings_current.setter
    def settings_current(self, args: Namespace):
        """
        iterate over all settings and set according to rcfile first, then cml args, 
        implicitily calls check_logic
        """
        rcdata = self.read_config(self.censorc_path)
        for setting in self._settings_current:
            # update settings with rcfile
            try:
                setting.value = rcdata[setting.type][setting.part][setting.name]
                rcdata[setting.type][setting.part].pop(setting.name)
            except KeyError:
                print(f"Could not find {setting} in rcdata.")
            
            # overwrite settings with cml arguments
            # FIXME - coverage between all actual settings and cml options?
            if hasattr(args, setting.name):
                setting.value = getattr(args, setting.name)
            else:
                print(f"Could not find {setting} in cml args.")
              
        print(f"Remaining settings:\n{rcdata}\n") # TODO - where to print out?
        
        warnings = self.check_logic()
        for warning in warnings:
            print(warning)
            
        if any([warning.fatal for warning in warnings]):
            print("\nFatal error in user input while setting up run settings. Going to exit!")
            sys.exit(1)
            

    def default_settings(self) -> SettingsTuple:
        """
        initialize settings as default
        """
        # workaround with a temporary list in order not having to redefine tuple operations
        # normal tuple operations cast result into normal tuple instead of SettingsTuple
        tmp = []
        # loop over all known settings (hardcoded in settings_options)
        for type_t, parts in CensoSettings.settings_options.items():
            for part, settings in parts.items():
                if settings:
                    # if settings exist create 'Setting' objects with default settings
                    for setting, definition in settings.items():
                        options = definition.get("options", None)
                        
                        # the typical type annotation errors
                        tmp.append(
                            Setting(
                                type_t, 
                                part, 
                                setting, 
                                definition["default"], 
                                options, 
                                definition["default"]
                            )
                        )
        
        return SettingsTuple(tmp)


    def read_config(self, censorc_path: str) -> SettingsDict:
        """
        Read from config data from file (here .censorc),
        every part has to be in rcfile
        """
        rcdata: SettingsDict = {}
        with open(censorc_path, "r") as csvfile:
            # skip header
            line = csvfile.readline()
            while not line.startswith("$GENERAL SETTINGS"):
                line = csvfile.readline()
            
            # split file into settings sections
            """
            sections begin with $
            part stated in this line is first word without $
            read lines until next $ line
            repeat until line is $END CENSORC
            """
            part = "general" # FIXME - cheeky workaround, user input may be problematic
            
            # mind the ordering of csvfile.readline(), should not lead to EOF errors
            while True:
                while not line.startswith("$"):
                    # split the line at ':' and remove leading and trailing whitespaces
                    spl = [s.strip() for s in line.split(":")]
                    sett_type = CensoSettings.get_type(spl[0]) # FIXME - eindeutig?
                    
                    # TODO - smells of copy paste
                    try:
                        # catch all cases for up until which level the dict is initialized
                        if sett_type and sett_type not in rcdata.keys():
                            # type-key not initialized
                            rcdata[sett_type] = {part: {spl[0]: sett_type(spl[1])}}
                        elif sett_type and part and part not in rcdata[sett_type].keys():
                            # part-key not initialized
                            rcdata[sett_type][part] = {spl[0]: sett_type(spl[1])}
                        elif sett_type and part and spl[0] != "" and spl[0] not in rcdata[sett_type][part].keys():
                            # setting-key not initialized
                            rcdata[sett_type][part][spl[0]] = sett_type(spl[1])
                    except Exception:
                        raise TypeError(f"Casting failed for line '{line}' while trying to read configuration file.")
                     
                    line = csvfile.readline()
                    
                if line.startswith("$END CENSORC"):
                    break
                
                # extract name of part from line
                part = line.split()[0][1:].lower()
                
                # read next line
                line = csvfile.readline()
            
        return rcdata

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
                outdata.write(f"ORCA: {self.external_paths['orcapath']}\n")
                outdata.write(f"ORCA version: {self.external_paths['orcaversion']}\n")
                outdata.write(f"GFN-xTB: {self.external_paths['xtbpath']}\n")
                outdata.write(f"CREST: {self.external_paths['crestpath']}\n")
                outdata.write(f"mpshift: {self.external_paths['mpshiftpath']}\n")
                outdata.write(f"escf: {self.external_paths['escfpath']}\n")
                outdata.write("\n")
                outdata.write("#COSMO-RS\n")
                outdata.write(f"{self.external_paths['cosmorssetup']}\n")
            else:
                # TODO - why is this set up like that (including/excluding binary)??
                # TODO - write some other default (e.g. "") instead of paths
                outdata.write("ORCA: /path/excluding/binary/\n")
                outdata.write("ORCA version: 4.2.1\n")
                outdata.write("GFN-xTB: /path/including/binary/xtb-binary\n")
                outdata.write("CREST: /path/including/binary/crest-binary\n")
                outdata.write("mpshift: /path/including/binary/mpshift-binary\n")
                outdata.write("escf: /path/including/binary/escf-binary\n")
                outdata.write("\n")
                outdata.write("#COSMO-RS\n")
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
                for parts in CensoSettings.settings_options.values():
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


    def find_rcfile(self, cwd: str, inprcpath=None) -> None:
        """check for existing censorc"""

        # looks for censorc file (global configuration file)
        # if there is no rcfile, CENSO exits
        # looks for custom path and standard paths:
        # cwd and $home dir
        # absolute path to file directly

        censorc_name = ".censorc"

        tmp = [
            os.path.join(cwd, censorc_name),
            os.path.join(os.path.expanduser("~"), censorc_name)
        ]

        # mapping the paths defined above to True/False, 
        # depending if file exists or not
        check = {
            os.path.isfile(tmp[0]): tmp[0],
            os.path.isfile(tmp[1]): tmp[1],
        }

        # FIXME - not the best solution to ensure code safety
        rcpath = ""

        # check for .censorc in standard locations if no path is given
        if not inprcpath:
            if all(list(check.keys())):
                # ask which one to use if both are found
                print(
                    f"Configuration files have been found, {tmp[0]} and "
                    f"{tmp[1]}. Which one to use? (cwd/home)"
                )
                
                user_input = ""
                while user_input.strip().lower() not in ("cwd", "home"):
                    print("Please type 'cwd' or 'home':")
                    user_input = input()
                
                if user_input.strip().lower() in ("cwd"):
                    rcpath = tmp[0]
                elif user_input.strip().lower() in ("home"):
                    rcpath = tmp[1]
                    
            elif any(list(check.keys())):
                # take the one file found
                rcpath = check[True]
        elif inprcpath and os.path.isfile(inprcpath):
            rcpath = inprcpath
        
        self.censorc_path = rcpath


    def read_program_paths(self):
        """
        Get absolute paths of external programs employed in censo
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
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # TODO - special part2
            """ if self.part2:
                # Handle sm2 --> solvent model in optimization:
                exchange_sm = {
                    "cosmo": "cpcm",
                    "cpcm": "cosmo",
                    "dcosmors": "smd",
                    "smd": "dcosmors",
                }
                if self.sm2 not in self.impsm2:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The solvent model {self.sm2} is not implemented!"
                    )
                    error_logical = True
                if self.prog2opt == "orca":
                    if self.sm2 in self.sm2_tm:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm2} is not available with "
                            f"{self.prog2opt}! Therefore {exchange_sm[self.sm2]} is used!"
                        )
                        self.sm2 = exchange_sm[self.sm2]
                    elif self.sm2 == "default":
                        self.sm2 = self.internal_defaults_orca["sm2"]["default"]
                if self.prog2opt == "tm":
                    if self.sm2 in self.sm2_orca:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm2} is not available with "
                            f"{self.prog2opt}! Therefore { exchange_sm[self.sm2]} is used!"
                        )
                        self.sm2 = exchange_sm[self.sm2]
                    elif self.sm2 == "default":
                        self.sm2 = self.internal_defaults_tm["sm2"]["default"] """
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # TODO - special part4
            """ if self.part4:
                # Handle sm4_j
                if self.prog4_j == "orca":
                    if self.sm4_j in self.sm4_j_tm:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_j} is not available with {self.prog4_j}!"
                            f" Therefore {exchange_sm[self.sm4_j]} is used!"
                        )
                        self.sm4_j = exchange_sm[self.sm4_j]
                    elif self.sm4_j == "default":
                        self.sm4_j = self.internal_defaults_orca["sm4_j"]["default"]
                if self.prog4_j == "tm":
                    if self.sm4_j in self.sm4_j_orca:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_j} is not available with {self.prog4_j}!"
                            f" Therefore {exchange_sm[self.sm4_j]} is used!"
                        )
                        self.sm4_j = exchange_sm[self.sm4_j]
                    elif self.sm4_j == "default":
                        self.sm4_j = self.internal_defaults_tm["sm4_j"]["default"]
                # Handle sm4_s
                if self.prog4_s == "orca":
                    if self.sm4_s in self.sm4_s_tm:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_s} is not available with {self.prog4_s}!"
                            f" Therefore { exchange_sm[self.sm4_s]} is used!"
                        )
                        self.sm4_s = exchange_sm[self.sm4_s]
                    elif self.sm4_s == "default":
                        self.sm4_s = self.internal_defaults_orca["sm4_s"]["default"]
                if self.prog4_s == "tm":
                    if self.sm4_s in self.sm4_s_orca:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_s} is not available with {self.prog4_s}!"
                            f" Therefore {exchange_sm[self.sm4_s]} is used!"
                        )
                        self.sm4_s = exchange_sm[self.sm4_s]
                    elif self.sm4_s == "default":
                        self.sm4_s = self.internal_defaults_tm["sm4_s"]["default"] """
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
            # FIXME - looks up solvents, if solvent not found for solvent model a replacement is chosen automatically
            """ else:
                for key, value in check_for.items():
                    if value:
                        if (
                            censo_solvent_db[self.solvent].get(key, "nothing_found")
                            == "nothing_found"
                        ):
                            self.save_errors.append(
                                f"{'ERROR:':{WARNLEN}}The solvent for solventmodel in "
                                "{key} is not found!"
                            )
                            error_logical = True
                        if key == "DC":
                            try:
                                if censo_solvent_db[self.solvent].get(key, None) is not None:
                                    _ = float(censo_solvent_db[self.solvent].get(key, None))
                                else:
                                    self.save_errors.append(
                                        f"{'ERROR:':{WARNLEN}}The dielectric constant for the solvent '{self.solvent}' "
                                        f"is not provided for the solventmodel {'cosmo / dcosmors'}!"
                                    )
                                    error_logical = True
                            except ValueError:
                                self.save_errors.append(
                                    f"{'ERROR:':{WARNLEN}}The dielectric constant can "
                                    "not be converted."
                                )
                                error_logical = True
                        elif key in ("smd", "cpcm"):
                            if (censo_solvent_db[self.solvent].get(key, ['', 'nothing_found'])[1].lower() 
                                not in getattr(self, lookup[key]) and
                                censo_solvent_db[self.solvent].get(key, ['', None])[1]):
                                self.save_errors.append(
                                    f"{'WARNING:':{WARNLEN}}The solvent "
                                    f"'{censo_solvent_db[self.solvent].get(key, 'nothing_found')[1]}'"
                                    f" for solventmodel/program {key} can not be checked "
                                    "but is used anyway."
                                )
                        else:
                            if (censo_solvent_db[self.solvent].get(key, ['', 'nothing_found'])[1]
                                not in getattr(self, lookup[key]) and
                                censo_solvent_db[self.solvent].get(key, ['', None])[1]):
                                self.save_errors.append(
                                    f"{'WARNING:':{WARNLEN}}The solvent "
                                    f"'{censo_solvent_db[self.solvent].get(key, 'nothing_found')[1]}' "
                                    f"for solventmodel/program {key} can not be checked "
                                    "but is used anyway."
                                )
                        if (key != "DC" and 
                            censo_solvent_db[self.solvent].get(key, ["",""])[0] is None and
                            censo_solvent_db[self.solvent].get(key, "nothing_found") != "nothing_found"
                            ):
                            if key == 'xtb':
                                tmp_sm = 'alpb'
                            else:
                                tmp_sm = key
                            if censo_solvent_db[self.solvent].get(key, ["",None])[1] is None:
                                self.save_errors.append(
                                    f"{'ERROR:':{WARNLEN}}The solvent '{self.solvent}' "
                                    f"is not parameterized for the solventmodel {tmp_sm}, and "
                                    f"no replacement is available!!!"
                                )
                                error_logical = True
                            else:
                                self.save_errors.append(
                                    f"{'WARNING:':{WARNLEN}}The solvent '{self.solvent}' "
                                    f"is not parameterized for the solventmodel {tmp_sm}, therefore"
                                    f" '{censo_solvent_db[self.solvent].get(key, ['',''])[1]}' is used!!!"
                                ) """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # TODO ?
        """# adjust functional-names to new naming convention (cfg.functional)
        self.func0 = self.func_info.relay_functionals.get(self.func0, self.func0)
        self.func = self.func_info.relay_functionals.get(self.func, self.func)
        self.func3 = self.func_info.relay_functionals.get(self.func3, self.func3)
        self.func_j = self.func_info.relay_functionals.get(self.func_j, self.func_j)
        self.func_s = self.func_info.relay_functionals.get(self.func_s, self.func_s)
        self.func_or = self.func_info.relay_functionals.get(self.func_or, self.func_or)
        self.func_or_scf = self.func_info.relay_functionals.get(
            self.func_or_scf, self.func_or_scf
        )
        # extracheck for r2scan-3c and ORCA
        try:
            orcaversion = int(self.external_paths['orcaversion'].split('.')[0])
        except (ValueError, AttributeError):
            orcaversion = 2"""
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
                
            if not (run is None or prog is None or func is None or basis is None):
                # iterate through settings and check for func
                # TODO - handle composite method bases
                if run and not stroptions[part] is None:
                    # FIXME - again incorrect error because pylance is too stupid
                    if not func in stroptions[part]["func"]["options"]:
                        warnings.append(LogicWarning(
                                "func",
                                f"The functional {func} was not found to be available for any software package in part {part}.",
                                ""
                            )
                        )
                    else:
                        # FIXME - doesn't work properly yet for part4_j/s
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

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # FIXME - mergable with later if
        # FIXME - should nmrmode work if either shieldings or couplings are not set?
        # TODO - fix when working on part4
        """if self.part4 and (self.couplings or self.shieldings):
            self.nmrmode = True

        if self.part4 and not self.couplings and not self.shieldings:
            self.part4 = False
            self.save_errors.append(
                f"{'INFORMATION:':{WARNLEN}}Neither coupling nor "
                "shielding constants are activated! Part4 is not executed."
            )
        elif self.part4 and not any(
            [
                getattr(self, flag)
                for flag in (
                    "h_active",
                    "c_active",
                    "f_active",
                    "si_active",
                    "p_active",
                )
            ]
        ):
            self.save_errors.append(
                f"{'INFORMATION:':{WARNLEN}}No active element for the calculation of NMR is "
                "activated in the .censorc! Therefore all nuclei are calculated!"
            )

        # no unpaired electrons in coupling or shielding calculations!
        if self.unpaired > 0:
            if self.part4 and (self.couplings or self.shieldings):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Coupling and shift calculations "
                    "(part4) are only available for closed-shell systems!"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """# Handle optlevel2:
        # sm2 needs to be set (not default!)
        if self.optlevel2 in ("None", None, "automatic"):
            if self.sm2 in ("smd", "dcosmors") and self.solvent != "gas":
                self.optlevel2 = "lax"
            else:
                # gas phase
                self.optlevel2 = "normal" """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        return warnings