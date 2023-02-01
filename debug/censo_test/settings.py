from argparse import Namespace
import os
import sys
import json
from types import MappingProxyType
from typing import Any, Callable, Dict, List, Union
from functools import reduce

from censo_test.cfg import (
    WARNLEN, 
    PARTS, 
    ASSETS_PATH, 
    PROGS, 
    SOLV_MODS, 
    SOLV_GMODS,
    Settings,
    __version__,
)
from censo_test.errors import LogicError
from censo_test.storage import CensoStorage

class DfaSettings:
    def __init__(self, obj: Dict[str, Dict[str, Dict]]):
        self.dfa_dict = obj

    
    def find_func(self, part):
        tmp = []
        for k, v in self.dfa_dict["functionals"].items():
            if part in v["part"]:
                tmp.append(k)
        
        return tmp


    def composite_bs(self) -> set:
        return set([v for v in self.dfa_dict["composite_method_basis"].values()])


class InternalSettings:
    """
    All options are saved here.
    Manages internal settings
    """
    
    # TODO - restart
    # must not be changed if restart(concerning optimization)
    restart_unchangeable = [
        "unpaired",
        "charge",
        "solvent",
        "prog",
        "prog2opt",
        "ancopt",
        "opt_spearman",
        "optlevel2",
        "func",
        "basis",
        "sm2",
        "nat",
        "radsize",
        "cosmorsparam",
    ]
    # may be changed but data may be lost/overwritten
    restart_changeable = {
        "multitemp": False,
        # "temperature": False, # should not be changeable all solvent and
        # rrho values depend on this
        "trange": False,
        "bhess": False,
        "part1_gfnv": False,
        "optimization_gfnv": False,
        "refinement_gfnv": False,
        "smgsolv1": False,
        "smgsolv2": False,
        "smgsolv3": False,
        "func_or": False,
        "basis_or": False,
        "func_or_scf": False,
        "freq_or": False,
        "func3": False,
        "basis3": False,
        "func_j": False,
        "basis_j": False,
        "sm4_j": False,
        "func_s": False,
        "basis_s": False,
        "sm4_s": False,
        # "consider_sym": calculated on the fly
    }

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
    solv_gmods = reduce(lambda x, y: x + y, SOLV_GMODS.values())

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
    # TODO - change min values to ranges?
    # TODO - find solution for None cases
    # TODO - solvent mapping
    settings_options = MappingProxyType({
        int: MappingProxyType({
            "general": MappingProxyType({
                "charge": {"default": 0, "min": -256}, # TODO - (re)move?
                "unpaired": {"default": 0, "min": 0}, # TODO - (re)move?
                "maxthreads": {"default": 1, "min": 1},
                "omp": {"default": 1, "min": 1}, # number of cores per thread
            }),
            "prescreening": None,
            "screening": None,
            "optimization": MappingProxyType({
                "optcycles": {"default": 8, "min": 2}, # TODO - which value as min?
                "radsize": {"default": 10, "min": 5}, # TODO - same
            }),
            "refinement": None,
            "nmr": None,
            "optrot": None,
            "uvvis": MappingProxyType({
                "nroots": {"default": 20, "min": 1}, # TODO - set dynamically
            }),
        }),
        float: MappingProxyType({
            "general": MappingProxyType({
                "imagethr": {"default": -100.0, "min": -300.0}, # TODO - threshold for imaginary frequencies
                "sthr": {"default": 0.0, "min": 0.0}, # TODO - what is this?
                "scale": {"default": 1.0, "min": 0.0}, # TODO - what is this?
                "temperature": {"default": 298.15, "min": 0.1}, # TODO - 0.0 still works?
            }),
            "prescreening": MappingProxyType({
                "prescreening_threshold": {"default": 4.0, "min": 1.0}, # TODO - which value as min?
            }),
            "screening": MappingProxyType({
                "screening_threshold": {"default": 3.5, "min": 0.75},
            }),
            "optimization": MappingProxyType({
                "opt_limit": {"default": 2.5, "min": 0.5}, # TODO - rename?
                "hlow": {"default": 0.01, "min": 0.01}, # TODO
                "optimization_P_threshold": {"default": 99.0, "min": 1.0}, # TODO
                "spearmanthr": {"default": 0.0, "min": 0.0},
            }),
            "refinement": MappingProxyType({
                "refinement_threshold": {"default": 99.0, "min": 1.0}, # TODO
            }),
            "nmr": MappingProxyType({
                "resonance_frequency": {"default": 300.0, "min": 150.0}, # TODO
            }),
            "optrot": None,
            "uvvis": MappingProxyType({
                "sigma": {"default": 0.1, "min": 0.1},
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
                "func": {"default": "b97-d", "options": tuple(dfa_settings.find_func("func0"))},
                "basis": {"default": "def2-SV(P)", "options": ("automatic",) + tuple(dfa_settings.composite_bs()) + ("def2-SV(P)", "def2-TZVP")},
                "prog": {"default": "tm", "options": PROGS},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
            }),
            "screening": MappingProxyType({
                "func": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func1"))},
                "basis": {"default": "automatic", "options": ("automatic",) + basis_sets},
                "prog": {"default": "tm", "options": PROGS},
                "smgsolv": {"default": "cosmors", "options": solv_gmods},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
            }),
            "optimization": MappingProxyType({
                "func": {"default": "r2scan-3c", "options": tuple(dfa_settings.find_func("func2"))},
                "basis": {"default": "automatic", "options": ("automatic",) + basis_sets},
                "prog": {"default": "tm", "options": PROGS},
                "prog2opt": {"default": "prog", "options": PROGS}, # TODO - ??? prog2 ??? # FIXME
                "sm": {"default": "default", "options": solv_mods}, # FIXME
                "smgsolv": {"default": "cosmors", "options": solv_gmods},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
                "optlevel2": {"default": "automatic", "options": ("crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme", "automatic")}, # TODO - what does this mean?
            }),
            "refinement": MappingProxyType({
                "prog": {"default": "tm", "options": PROGS},
                "func": {"default": "pw6b95", "options": tuple(dfa_settings.find_func("func3"))},
                "basis": {"default": "def2-TZVPD", "options": basis_sets},
                "smgsolv": {"default": "cosmors", "options": solv_gmods},
                "gfnv": {"default": "gfn2", "options": ("gfn1", "gfn2", "gfnff")},
            }),
            "nmr": MappingProxyType({
                "prog4_j": {"default": "tm", "options": PROGS},
                "func_j": {"default": "pbe0", "options": tuple(dfa_settings.find_func("func_j"))},
                "basis_j": {"default": "def2-TZVP", "options": basis_sets},
                "sm4_j": {"default": "default", "options": solv_mods}, # FIXME
                "prog4_s": {"default": "tm", "options": PROGS},
                "func_s": {"default": "pbe0", "options": tuple(dfa_settings.find_func("func_s"))},
                "basis_s": {"default": "def2-TZVP", "options": basis_sets},
                "sm4_s": {"default": "default", "options": solv_mods}, # FIXME
                "h_ref": {"default": "TMS", "options": ("TMS",)},
                "c_ref": {"default": "TMS", "options": ("TMS",)},
                "f_ref": {"default": "CFCl3", "options": ("CFCl3",)},
                "si_ref": {"default": "TMS", "options": ("TMS",)},
                "p_ref": {"default": "TMP", "options": ("TMP", "PH3")},
            }),
            "optrot": MappingProxyType({
                "func": {"default": "pbe", "options": tuple(dfa_settings.find_func("func_or"))},
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


    @staticmethod
    def get_type(setting) -> Union[type, None]:
        """
        returns the type of a given setting
        """
        for type_t, parts in InternalSettings.settings_options.items():
            for settings in parts.values():
                if settings:
                    if setting in settings.keys():
                        return type_t

        return None

        
    def __init__(self, storage):

        self._settings_current: Settings

        self.storage: CensoStorage = storage

        # no censorc found at standard dest./given dest.
        if self.storage.censorc_path == "":
            print(
                f"No rcfile has been found. Do you want to create a new one?\n"
            )

            user_input = ""
            while user_input.strip().lower() not in ["yes", "y", "no", "n"]:
                print("Please type 'yes/y' or 'no/n':")
                user_input = input()
            
            if user_input.strip().lower() in ("y", "yes"):
                self.storage.censorc_path = self.write_config()
                self.storage.censorc_name = "censorc_new"
            elif user_input.strip().lower() in ("n", "no"):
                print(
                    "Configuration file needed to run CENSO!\n"
                    "Going to exit!"
                )
                sys.exit(1)

        self.settings_current = self.storage.args
        
        # print errors and exit if there are any conflicting settings
        errors = self.check_logic()
        for error in errors:
            print(error)

        if len(errors) != 0:
            sys.exit(1) 
            

    # TODO - write special timer decorator for this function
    @property
    def settings_current(self) -> Callable:
        """
        returns settings for given types and/or parts structured as self.settings_current
        can also find value of single setting
        property used as proxy for actual filter functionality in order to still use decorator functionality
        usage: <instance>.settings_current(types=..., parts=..., settings=...) or omit some/all keywords
        inclusion of more redundant keywords avoids unnecessary recursions
        """
        
        # TODO - check runtime intensity and optimize if problematic
        
        def f(types=set(self._settings_current.keys()), parts=set(PARTS), settings=set([])) -> Any:

            # declare filter mapping with correct typing
            filters: dict[str, set[str]] = {}
            
            # if the keyword value is iterable, assert that it is a set
            for keyword, keyname in zip((parts, types, settings), ("parts", "types", "settings")):
                try:
                    if iter(keyword) and type(keyword) != set:
                        # convert str to a set of itself in order not to get a set of chars
                        if type(keyword) == str:
                            keyword = set([keyword])
                        else:
                            keyword = set(keyword)
                # keyword not iterable, make it iterable by putting it into a set
                except TypeError:
                    keyword = set([keyword])
                
                # use the strings from the zip as keys for the filters dict
                filters[keyname] = keyword
                    
            # recursively filter input dictionary
            def filter_f(in_dict: dict = self._settings_current, i: int = 0) -> Any:
                fmap = ["types", "parts", "settings"]
                out_dict = {}

                # go through in_dict and filter out values for keys, which are
                # defined in filters
                for key, value in in_dict.items():
                    if key in filters[fmap[i]] and i < 2 and value:
                        out_dict[key] = filter_f(value, i+1)
                        
                        # remove key if there is nothing inside
                        try:
                            if len(out_dict[key]) == 0:
                                out_dict.pop(key)
                        except TypeError:
                            pass
                    # if we are at the level of settings
                    elif key in filters[fmap[i]] and i == 2:
                        out_dict[key] = value
                 
                # collapse to smallest possible nesting    
                # partX: only one setting -> returns value of that setting
                # accessing: settingsdict[part] = value of that setting   
                if len(out_dict) == 1:
                    return reduce(lambda x, y: x + y, iter(out_dict.items()).__next__()[::-1])
                else:
                    return out_dict, None
                       
            return filter_f()

        return f
    
    # sets _settings_current as opposed to name, in order to be able
    # to declare it's type beforehand
    @settings_current.setter
    def settings_current(self, args: Namespace):
        """
        iterate over all settings and set according to rcfile first, then cml args, 
        else set to default 
        throw error if setting not allowed in options
        """
        # TODO - wait for fix for compatibility
        """ ### if restart read all settings from previous run (enso.json)
        if args.restart and os.path.isfile(os.path.join(self.parent.cwd, "enso.json")):
            tmp = config.read_json(os.path.join(self.parent.cwd, "enso.json"), silent=True)
            previous_settings = tmp.get("settings")
            for key, value in previous_settings.items():
                if vars(self.args).get(key, "unKn_own") == "unKn_own":
                    # print(key, 'not_known')
                    continue
                if getattr(self.args, key, "unKn_own") is None:
                    setattr(self.args, key, value)
        ### END if restart """
        
        # get settings from rcfile first
        self._settings_current = self.read_config()
        
        # overwrite settings
        for type_t, partsd in InternalSettings.settings_options.items():
            for part, settings in partsd.items():
                if settings:
                    # set the value of the settings according to cml if given
                    for setting in settings.keys():
                        # FIXME - coverage between all actual settings and cml options?
                        try:
                            if getattr(args, setting):
                                self._settings_current[type_t][part][setting] = getattr(args, setting)
                        except AttributeError:
                            pass


    def read_config(self) -> Settings:
        """
        Read from config data from file (here .censorc),
        every part has to be in rcfile
        """
        rcdata: Settings = {}
        with open(self.storage.censorc_path, "r") as csvfile:
            # skip header
            line = csvfile.readline()
            while not line.startswith("$PRESCREENING SETTINGS"):
                line = csvfile.readline()
            
            # split file into settings sections
            """
            sections begin with $
            part stated in this line is first word without $
            read lines until next $ line
            repeat until line is $END CENSORC
            """
            part = "prescreening" # FIXME - cheeky workaround, user input may be problematic
            # mind the ordering of csvfile.readline(), should not lead to EOF errors
            while True:
                while not line.startswith("$"):
                    spl = line.strip().split(":")
                    sett_type = InternalSettings.get_type(spl[0])
                    if sett_type and sett_type not in rcdata.keys():
                        try:
                            rcdata[sett_type] = {}
                            rcdata[sett_type][part] = {}
                            rcdata[sett_type][part][spl[0]] = sett_type(spl[1])
                        except Exception:
                            """casting failed"""
                            # TODO
                     
                    line = csvfile.readline()
                    
                if line.startswith("$END CENSORC"):
                    break
                
                part = line.split()[0][1:].lower()
                line = csvfile.readline()
            
        return rcdata


    def write_config(self) -> str:
        """
        write new configuration file with default settings into 
        either the local or ~/.censorc directory.
        """

        # FIXME - is this legal
        """ if self.args.copyinput:
            self.write_config(global_config=False, usepaths=True)
            print(
                "The file censo.inp with the current settings has been written to "
                "the current working directory."
            )
            print("\nGoing to exit!")
            sys.exit(0) """

        # what to do if there is an existing configuration file
        if not (self.storage.censorc_path or self.storage.censorc_path == ""):
            print(
                f"An existing configuration file ({self.storage.censorc_name})" 
                f"has been found in {self.storage.censorc_path}"
                f"Do you want to copy existing program path information to the "
                f"new remote configuration file?"
            )
            
            user_input = ""
            while user_input.strip().lower() not in ["yes", "y", "no", "n"]:
                print("Please type 'yes/y' or 'no/n':")
                user_input = input()
            
            if user_input.strip().lower() in ("y", "yes"):
                self.storage.read_program_paths()

        # write new censorc
        if self.storage.args.inprcpath:
            path = self.storage.args.inprcpath
        # (path is never unbound)
        else:
            path = os.path.join(self.storage.cwd, "censorc_new")
            
        print(f"Writing configuration file to {path} ...")
        with open(path, "w", newline=None) as outdata:
            outdata.write(f"$CENSO global configuration file: {self.storage.censorc_name}\n")
            outdata.write(f"$VERSION:{__version__} \n")
            outdata.write("\n")
            
            # if ALL program paths are set, write the paths to the new rcfile
            if all([s != "" for s in self.storage.external_paths.values()]):
                # write stored program paths to file
                outdata.write(f"ORCA: {self.storage.external_paths['orcapath']}\n")
                outdata.write(f"ORCA version: {self.storage.external_paths['orcaversion']}\n")
                outdata.write(f"GFN-xTB: {self.storage.external_paths['xtbpath']}\n")
                outdata.write(f"CREST: {self.storage.external_paths['crestpath']}\n")
                outdata.write(f"mpshift: {self.storage.external_paths['mpshiftpath']}\n")
                outdata.write(f"escf: {self.storage.external_paths['escfpath']}\n")
                outdata.write("\n")
                outdata.write("#COSMO-RS\n")
                outdata.write(f"{self.storage.external_paths['cosmorssetup']}\n")
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
                
                # get all parts (partsd) for all types
                for partsd in InternalSettings.settings_options.values():
                    # only do this loop if the current part has settings at all
                    if partsd[part]:
                        # get defaults
                        for setting, definition in partsd[part].items():
                            outdata.write(f"{setting}: {definition['default']}\n")
                
            outdata.write("\n$END CENSORC\n")

        print(
            "\nA new configuration file was written into the current directory file:\n"
            "censorc_new\n"
            "You have to adjust the settings to your needs"
            " and it is mandatory to correctly set the program paths!\n"
            "Additionally move the file to the correct filename: '.censorc'\n"
            "and place it either in your /home/$USER/ or current directory.\n"
            "\nAll done!"
        )

        return path


    def check_logic(self) -> List[LogicError]:
        """
        Checks internal settings for impossible setting-combinations
        also checking if calculations are possible with the requested qm_codes.
        """
        # TODO - move checks concerning options given by commandline
        # TODO - check for uv/vis args
        # TODO - pass parts to be checked via arg
        # TODO - divide into sections
        # TODO - collect all errors and return them
        
        errors = []
        
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
        if not self.settings_current(settings="gas-phase")[0]:
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
            tmp_settings, _ = self.settings_current(
                types=str, 
                parts=PARTS, 
                settings=["prog", "sm", "smgsolv"]
            )
            for part in PARTS:
                try:
                    prog: str = tmp_settings[part]["prog"]
                    sm: str = tmp_settings[part]["sm"]
                    smgsolv: str = tmp_settings[part]["smgsolv"]
                except KeyError:
                    continue
                
                # availability in program
                # FIXME - calling the attributes of the job types should work?
                if (
                    sm not in SOLV_MODS[prog]
                    or smgsolv not in SOLV_GMODS[prog]
                ):
                    errors.append(LogicError(
                            "sm/smgsolv",
                            f"{sm} or {smgsolv} not available with given program {prog} in part {part}.",
                            "Choose a different solvent model!"
                        )
                    )
                # availability in solvent model
                elif self.settings_current(settings="solvent") not in solvents[sm]:
                    errors.append(LogicError(
                            "solvent",
                            f"Solvent {self.settings_current(settings='solvent')} not available for solvent model {sm} in part {part}.",
                            "Choose a different solvent!"
                        )
                    )
                    # TODO - check for dielectric constant for cosmo ->  dc only needed for dcosmors?
                    
                # there is no multitemp support for any orca smgsolv
                if (
                    self.settings_current(parts=part, settings="prog") == "orca"
                    and self.settings_current(settings="multitemp")
                ):
                    errors.append(LogicError(
                            "multitemp",
                            "There is no multitemp support for any solvent models in ORCA.",
                            "Disable multitemp!"                  
                        )
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
            # (why would you want that)
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
        tmp_stroptions = InternalSettings.settings_options[str]

        for part in PARTS:
            """
            get settings for part
            look for funcX/basisX
            check if funcX/basisX values are allowed in settings_options
            """
            # get settings for the part as dict
            tmp_settings = self.settings_current(parts=part, settings=["run", "func", "basis"])

            # iterate through settings and check for funcX
            # TODO - handle composite method bases
            if tmp_settings[bool] and tmp_stroptions[part]:
                try:
                    func = tmp_settings["func"]
                except KeyError:
                    func = tmp_settings[str]
                
                # FIXME - again incorrect error because pylance is too stupid
                if not func in tmp_stroptions[part]["func"]["options"]:
                    errors.append(LogicError(
                            "func",
                            f"The functional {func} is not available for any software package in part {part}.",
                            "Choose a different functional!"
                        )
                    )
                else:
                    # FIXME - doesn't work properly yet for part4_j/s
                    tmp_prog_funcname = InternalSettings.dfa_settings.dfa_dict["functionals"][tmp_settings["func"]][tmp_settings["prog"]]
                    if not tmp_prog_funcname:
                        errors.append(LogicError(
                                "func",
                                f"The functional {func} is not available for the software package {tmp_settings['prog']} in part {part}.",
                                "Choose a different functional or change the software package!"
                            )
                        )
                    # extra check for r2scan-3c
                    elif (
                        func == "r2scan-3c" 
                        and tmp_settings["prog"] == "orca"
                        and int(self.storage.external_paths["orcaversion"].split(".")[0]) < 5
                    ):
                        errors.append(LogicError(
                                "func",
                                f"The functional r2scan-3c is only available from ORCA version 5 in part {part}.",
                                "Choose a different functional or use a newer version of ORCA!"    
                            )
                        )
                    # extra check for dsd-blyp
                    elif (
                        func in ("dsd-blyp", "dsd-blyp-d3")
                        and tmp_settings["basis"] != "def2-TZVPP"
                    ):
                        errors.append(LogicError(
                                "func",
                                f"The functional dsd-blyp is only available with the def2-TZVPP basis set in part {part}.",
                                "Change the basis set to def2-TZVPP or choose a different functional!"    
                            )
                        )
                        
                    if not tmp_settings["basis"] in tmp_stroptions[part]["basis"]["options"]:
                        errors.append(LogicError(
                                "basis",
                                f"The basis set {tmp_settings['basis']} is not available for any software package in part {part}.",
                                "Choose a different basis set!"    
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
            
        return errors