"""
Wrapper to manage external program calls as well as setting up CENSO enviroment for runtime
store paths here
"""

import os
import sys
import json
import shutil
from argparse import Namespace
from collections import OrderedDict
from typing import Dict, Union
import weakref
from numpy import exp


from censo.cfg import (
    CODING,
    WARNLEN,
    __version__,
)
from censo.orca_job import OrcaJob
from censo.settings import InternalSettings
from censo.inputhandling import config_setup
from censo.tm_job import TmJob
from censo.utilities import (
    mkdir_p,
    do_md5,
    t2x,
    print,
)

# TODO - how do the assets files get into ~/.censo_assets?
class CensoCore:
    _instance_ref = weakref.WeakValueDictionary()
    
    prog_job: dict[str, type] = {"tm": TmJob, "orca": OrcaJob}

    @staticmethod
    def factory(cwd, args: Namespace):
        """keep track of the CensoCore instance (should only be one)"""
        instance = CensoCore(cwd, args)
        CensoCore._instance_ref[id(instance)] = instance
        return instance


    @staticmethod
    def check_instance(f):
        """check number of CensoCore instances, exit if 0 or more than 1 to assure CENSO runs properly"""
        def wrap(*args, **kwargs):
            n_inst = len(CensoCore._instance_ref)
            if n_inst > 1:
                print(
                    f"{'ERROR:':{WARNLEN}}There should be only one instance of censo_core"
                    "Current instances are:"
                )
                for ref in CensoCore._instance_ref.values():
                    print(ref)
                sys.exit(1)
            elif n_inst == 0:
                print(
                    f"{'ERROR:':{WARNLEN}}There is currently no instance of censo_core"
                )
                sys.exit(1)
            else:
                return f(*args, **kwargs)

        return wrap


    @staticmethod 
    def core():
        try:
            core: CensoCore = next(CensoCore._instance_ref.values())
            return core
        except StopIteration:
            print(
                f"{'ERROR:':{WARNLEN}}There is currently no instance of censo_core"
            )
            sys.exit(1)


    def __init__(self, cwd, args):
        """initialize core"""
        self.args: Namespace = args
        self.cwd: str = cwd
        self.censorc_name = ".censorc"
        
        # looks for censorc file (global configuration file, might be omitted)
        self.censorc_path: Union[str, None] = self.find_rcfile()
        
        # if no path is found, CENSO exits (assets are essential for functionality)
        self.assets_path: str = self.find_assets()
        
        # if no input ensemble is found, CENSO exits
        self.ensemble_path: str = self.find_ensemble()

        self.internal_settings = InternalSettings(self.args)

        # pathsdefaults: --> read_program_paths
        self.external_paths: dict[str, str] = {}
        self.external_paths["orcapath"] = ""
        self.external_paths["orcaversion"] = ""
        self.external_paths["xtbpath"] = ""
        self.external_paths["crestpath"] = ""
        self.external_paths["cosmorssetup"] = ""
        self.external_paths["dbpath"] = ""
        self.external_paths["cosmothermversion"] = ""
        self.external_paths["mpshiftpath"] = ""
        self.external_paths["escfpath"] = ""
        

    def setup_censo(self) -> None:
        """setup external and internal settings"""
        try:
            self.args
            self.cwd
        except AttributeError:
            print(f"{'ERROR:':{WARNLEN}}Error while accessing cml arguments!")
            print("\nGoing to exit!")
            sys.exit(1)

        # TODO - collapse if conditions to function calls by iteration over args attributes

        if self.args.cleanup:
            self.cleanup_run()
            print("Removed files and going to exit!")
            sys.exit(0)
        elif self.args.cleanup_all:
            self.cleanup_run(True)
            print("Removed files and going to exit!")
            sys.exit(0)

        if self.args.writeconfig:
            self.write_config(True)
            sys.exit(0)

        if self.args.inprcpath:
            try:
                self.find_rcfile(os.path.abspath(os.path.expanduser(self.args.inprcpath)))
            except FileNotFoundError:
                print(
                    f"{'ERROR:':{WARNLEN}}Could not find the configuration file: {self.censorc_name}."
                )
                print("\nGoing to exit!")
                sys.exit(1)
        # TODO - leave until compatibility is fixed
        elif (
            self.args.restart
            and os.path.isfile(os.path.join(config.cwd, "enso.json"))
            and os.path.isfile(
                config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
                .get("settings", {})
                .get("configpath", "")
            )
        ):
            # read configpath from previous enso.json if restart is requested
            # additionally check for vapor pressure
            tmp = config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
            config.configpath = tmp.get("settings", {}).get("configpath", "")
            if not self.args.vapor_pressure:
                self.args.vapor_pressure = tmp.get("settings", {}).get("vapor_pressure", False)
        elif os.path.isfile(os.path.join(config.cwd, configfname)):
            # local configuration file before remote configuration file
            config.configpath = os.path.join(config.cwd, configfname)
        elif os.path.isfile(os.path.join(os.path.expanduser("~"), configfname)):
            # remote configuration file
            config.configpath = os.path.join(os.path.expanduser("~"), configfname)
        else:
            print(
                f"{'ERROR:':{WARNLEN}}Could not find the configuration file: {configfname}.\n"
                f"{'':{WARNLEN}}The file has to be either in /home/$USER/ or the current "
                "working directory!\n"
                "The configuration file can be otherwise directly referenced using: "
                "'censo -inprc /path/to/.censorc'"
            )
            print("\nGoing to exit!")
            sys.exit(1)
        
        # TODO - how to manage user customizable solvent db?
        """
        solvent_user_path = os.path.expanduser(
            os.path.join("~/.censo_assets/", "censo_solvents.json")
        )

        if os.path.isfile(solvent_user_path):
            config.save_infos.append(
                "Reading file: {}\n".format(os.path.basename(solvent_user_path))
            )
            try:
                with open(solvent_user_path, "r", encoding=CODING, newline=None) as inp:
                    censo_solvent_db.update(json.load(inp, object_pairs_hook=OrderedDict))
            except (ValueError, TypeError, FileNotFoundError):
                print(
                    f"{'ERROR:':{WARNLEN}}Your censo_solvents.json file in {solvent_user_path} is corrupted!\n"
                    f"{'':{WARNLEN}}You can delete your corrupted file and a new "
                    f"censo_solvents.json will be created on the next start of CENSO."
                )
                print("\nGoing to exit!")
                sys.exit(1)

            #vapor_pressure
            if self.args.vapor_pressure == "on" or self.args.vapor_pressure is True:
                config.vapor_pressure = True
                if os.path.isfile(os.path.join(config.cwd, "same_solvent.json")):
                    try:
                        with open(os.path.join(config.cwd, "same_solvent.json"), "r", encoding=CODING, newline=None) as inp:
                            censo_solvent_db.update(json.load(inp, object_pairs_hook=OrderedDict))
                    except (ValueError, TypeError, FileNotFoundError):
                        print(f"{'ERROR:':{WARNLEN}}Your file same_solvents.json is corrupted!")
                        print("\nGoing to exit!")
                        sys.exit(1)
        else:
            try:
                with open(solvent_user_path, "w") as out:
                    json.dump(censo_solvent_db, out, indent=4, sort_keys=True)
                config.save_infos.append(
                    "Creating file: {}".format(os.path.basename(solvent_user_path))
                )
            except OSError as error:
                print(f"{'WARNING:':{WARNLEN}}The file {solvent_user_path} could not be "
                        "created and internal defaults will be applied!"
                )
                print(f"{'INFORMATION:':{WARNLEN}}The corresponding python error is: {error}.")
        ### END solvent database adjustable by user"""
        # TODO - same here
        """### NMR reference shielding constant database adjustable by user
        nmr_ref_user_path = os.path.expanduser(
            os.path.join("~/.censo_assets/", "censo_nmr_ref.json")
        )
        if not os.path.isfile(nmr_ref_user_path):
            try:
                with open(nmr_ref_user_path, "w") as out:
                    tmp = NmrRef()
                    json.dump(tmp, out, default=NmrRef.NMRRef_to_dict, indent=4, sort_keys=True)
                config.save_infos.append(
                    "Creating file: {}".format(os.path.basename(nmr_ref_user_path))
                )
            except OSError as error:
                print(f"{'WARNING:':{WARNLEN}}The file {nmr_ref_user_path} could not be "
                        "created and internal defaults will be applied!"
                )
                print(f"{'INFORMATION:':{WARNLEN}}The corresponding python error is: {error}.")
        ### END NMR reference shielding constant database adjustable by user"""

        # TODO - again similar thing
        """### ORCA user editable input
        orcaeditablepath = os.path.expanduser(
            os.path.join("~/.censo_assets/", "censo_orca_editable.dat")
        )
        if not os.path.isfile(orcaeditablepath):
            try:
                with open(orcaeditablepath, "w", newline=None) as out:
                    for line in editable_ORCA_input.get("default", []):
                        out.write(line+'\n')
                config.save_infos.append(
                    "Creating file: {}".format(os.path.basename(orcaeditablepath))
                )
            except OSError as error:
                print(f"{'INFORMATION:':{WARNLEN}}The file {orcaeditablepath} could not be "
                        "created and internal defaults will be applied!"
                )
                print(f"{'INFORMATION:':{WARNLEN}}The corresponding python error is: {error}.")
        else: 
            if os.path.isfile(orcaeditablepath):
                try:
                    config.save_infos.append(
                    "Reading file: {}\n".format(os.path.basename(orcaeditablepath))
                    )
                    tmp_orca = []
                    with open(orcaeditablepath, "r", encoding=CODING) as inp:
                        for line in inp:
                            if '#' in line:
                                tmp_orca.append(line.split("#")[0].rstrip())
                            else:
                                tmp_orca.append(line.rstrip())
                        editable_ORCA_input.update({'default':tmp_orca})
                except (ValueError, TypeError, FileNotFoundError, OSError):
                    print(
                        f"{'INFORMATION:':{WARNLEN}}Your {os.path.basename(orcaeditablepath)} file in "
                        f"{orcaeditablepath} is corrupted and internal defaults will be applied!\n"
                        f"{'':{WARNLEN}}You can delete your corrupted file and a new "
                        f"one will be created on the next start of CENSO."
                    )
        ### END ORCA user editable input"""

        # TODO - wait for fix for compatibility
        ### if restart read all settings from previous run (enso.json)
        if self.args.restart and os.path.isfile(os.path.join(config.cwd, "enso.json")):
            tmp = config.read_json(os.path.join(config.cwd, "enso.json"), silent=True)
            previous_settings = tmp.get("settings")
            for key, value in previous_settings.items():
                if vars(self.args).get(key, "unKn_own") == "unKn_own":
                    # print(key, 'not_known')
                    continue
                if getattr(self.args, key, "unKn_own") is None:
                    setattr(self.args, key, value)
        ### END if restart

        # TODO - fix compatibility here
        # read censorc for settings
        if config.configpath:
            # combine args und commandline
            # check if startread in file:
            startread = "$CRE SORTING SETTINGS:"
            with open(config.configpath, "r") as myfile:
                try:
                    data = myfile.readlines()
                    censorc_version = "0.0.0"
                    for line in data:
                        if "$VERSION" in line:
                            censorc_version = line.split(":")[1]
                    if int(censorc_version.split(".")[1]) < int(
                        __version__.split(".")[1]
                    ) or int(censorc_version.split(".")[0]) < int(
                        __version__.split(".")[0]
                    ):
                        print(
                            f"{'ERROR:':{WARNLEN}}There has been an API break and you have to "
                            f"create a new .censorc.\n{'':{WARNLEN}}E.g. 'censo -newconfig'"
                        )
                        print(
                            f"{'INFORMATION:':{WARNLEN}}Due to the API break data in ~/.censo_assets/ might need updating.\n"
                            f"{'':{WARNLEN}}Therefore delete the files {'censo_nmr_ref.json'} and {'censo_solvents.json'} so "
                            "that they can be re-created upon the next censo call."
                        )
                        sys.exit(1)
                    myfile.seek(0)  # reset reader
                except (ValueError, KeyError, AttributeError) as e:
                    print(e)
                    print(
                        f"{'ERROR:':{WARNLEN}}Please create a new .censorc --> 'censo -newconfig'"
                    )
                    sys.exit(1)
                if not startread in myfile.read():
                    print(
                        f"{'ERROR:':{WARNLEN}}You are using a corrupted .censorc. Create a new one!"
                    )
                    sys.exit(1)
            # combine commandline and .censorc
            config.read_config(config.configpath, startread, args)

        if self.args.copyinput:
            self.read_program_paths()
            self.write_config(global_config=False, usepaths=True)
            print(
                "The file censo.inp with the current settings has been written to "
                "the current working directory."
            )
            print("\nGoing to exit!")
            sys.exit(0)

        self.read_program_paths(silent=True)
        self.read_input()

        # FIXME - temporary place for remaining settings
        if not self.args.spearmanthr:
            # set spearmanthr by number of atoms:
            self.spearmanthr = 1 / (exp(0.03 * (self.internal_settings.runinfo["nat"] ** (1 / 4))))

        self.internal_settings.runinfo["consider_unconverged"] = False


    def find_rcfile(self, path=None) -> Union[str, None]:
        """check for existing censorc"""

        if not path:
            # check if .censorc in local or home dir
            if os.path.isfile(os.path.join(self.cwd, self.censorc_name)):
                censorc_path = os.path.join(self.cwd, self.censorc_name)
            elif os.path.isfile(os.path.join(os.path.expanduser("~"), self.censorc_name)):
                censorc_path = os.path.join(os.path.expanduser("~"), self.censorc_name)
            else:
                censorc_path = None
        else:
            censorc_path = path
            if not os.path.isfile(censorc_path):
                raise FileNotFoundError

        return censorc_path

    
    def find_assets(self) -> str:
        """
        look for assets folder
        if it is not found, CENSO exits
        """
        assets_path = os.path.expanduser("~/.censo_assets")
        if not os.path.isdir(assets_path):
            print(f"{'WARNING:':{WARNLEN}}The folder '~/.censo_assets/' designed for additional "
                "remote configuration files can not be found!\n"
            )
            sys.exit(1)

        return assets_path


    def find_ensemble(self) -> str:
        """check for ensemble input file"""
        # if input file given via args use this path, otherwise set path to a default value
        if self.args.inp:
            ensemble_path = os.path.join(self.cwd, self.args.inp)
        else:
            ensemble_path = os.path.join(self.cwd, "crest_conformers.xyz")

        if os.path.isfile(ensemble_path):
            return ensemble_path
        else:
            print(
                f"{'ERROR:':{WARNLEN}}The input ensemble cannot be found!"
            )
            sys.exit(1)
        

    def orca_version(self):
        """get the version of orca given via the program paths"""


    # TODO - optimize
    def read_input(self) -> None:
        """
        read ensemble input file (e.g. crest_conformers.xyz)
        look for enso.json and load if found
        """
        
        if os.path.isfile(os.path.join(self.cwd, "enso.json")):
            tmp = config.read_json(os.path.join(self.cwd, "enso.json"), silent=True)
            if "ensemble_info" in tmp and self.args.inp is None:
                inpfile = os.path.basename(tmp["ensemble_info"].get("filename"))

                if os.path.isfile(inpfile):
                    self.args.inp = inpfile

                if self.args.debug:
                    print(f"Using Input file from: {inpfile}")

        if self.ensemble_path:
            # store md5 hash for quick comparison of inputs later
            self.internal_settings.runinfo["md5"] = do_md5(self.ensemble_path)
            with open(self.ensemble_path, "r", encoding=CODING, newline=None) as inp:
                foundcoord = False
                for line in inp:
                    if "$coord" in line:
                        foundcoord = True
                        break
            # if $coord in file => tm format
            if foundcoord:
                _, self.internal_settings.runinfo["nat"] = t2x(
                    self.ensemble_path, writexyz=True, outfile="converted.xyz"
                )
                self.ensemble_path = os.path.join(self.cwd, "converted.xyz")
                self.internal_settings.runinfo["maxconf"] = 1
                self.internal_settings.runinfo["nconf"] = 1
            # else xyz format (hopefully)
            else:
                with open(
                    self.ensemble_path, "r", encoding=CODING, newline=None
                ) as infile:
                    try:
                        self.internal_settings.runinfo["nat"] = int(infile.readline().split()[0])
                        filelen = 1
                    except (ValueError, TypeError, IndexError) as e:
                        print(
                            f"{'ERROR:':{WARNLEN}}Could not get the number of atoms or the "
                            "number of conformers from the inputfile "
                            f"{os.path.basename(self.args.inp)}"
                        )
                        sys.exit(1)
                    for line in infile:
                        filelen += 1

                    # check if coords and nat are given in proper xyz format for all conformers
                    try:
                        self.internal_settings.runinfo["maxconf"] = int(filelen / (self.internal_settings.runinfo["nat"] + 2))
                        if filelen % (self.internal_settings.runinfo["nat"] + 2) != 0:
                            raise ValueError
                    except ValueError:
                        print(
                            f"{'ERROR:':{WARNLEN}}Could not get the number of atoms or the "
                            "number of conformers from the inputfile "
                            f"{os.path.basename(self.args.inp)}"
                        )
                        sys.exit(1)
        else:
            print(f"{'ERROR:':{WARNLEN}}The input file can not be found!")
            sys.exit(1)


    # FIXME - new arg needed?
    def write_config(self, new=False, global_config=True, usepaths=False, update=False) -> None:
        """
        write new local or global configuration file into either the local or .censorc directory.
        """

        if new and global_config:
            if self.censorc_path is not None:
                print(f"An existing .censorc has been found in {self.censorc_path}")
                print(
                    f"Do you want to copy existing program path information to the "
                    f"new remote configuration file?"
                )

                user_input = ""
                while user_input.strip().lower() not in ["yes", "y", "no", "n"]:
                    print("Please type 'yes/y' or 'no/n':")
                    user_input = input()
                
                if user_input.strip().lower() in ("y", "yes"):
                    usepaths=True
                elif user_input.strip().lower() in ("n", "no"):
                    usepaths=False

            print("\nPlease chose your QM code applied in parts 0-2 either TURBOMOLE (TM) or ORCA (ORCA) or decide later (later):") # FIXME - decide later option where?
            user_input = ""
            while user_input.strip().lower() not in ["tm", "orca"]:
                user_input = input()

            if user_input.strip().lower() in ('tm', 'orca'):
                # TODO - manage this via internal_settings instance
                if user_input.strip().lower() in ('tm'):
                    config = config_setup(path=os.path.abspath(self.cwd), **{'prog':'tm'})
                elif user_input.strip().lower() in ('orca'):
                    config = config_setup(path=os.path.abspath(self.cwd), **{'prog':'orca'})

            if usepaths:
                self.read_program_paths()
            # write new censorc
            new_censorc_name = "censorc_new"
            # TODO - merge with write_rcfile
            config.write_rcfile(os.path.join(config.cwd, new_censorc_name), usepaths=usepaths, update=True)
            print(
                "\nA new ensorc was written into the current directory file: "
                f"{new_censorc_name}!\nYou have to adjust the settings to your needs"
                " and it is mandatory to correctly set the program paths!\n"
                "Additionally move the file to the correct filename: '.censorc'\n"
                "and place it either in your /home/$USER/ or current directory.\n"
                "\nAll done!"
            )

        # TODO - avoid passing censorc_path = None

        args_key = {v: k for k, v in self.internal_settings.key_args_dict.items()}
        if update:
            data = self.internal_settings.provide_runinfo(extend=False) # TODO - fix provide_runinfo
        else:
            data = {}

        if global_config:
            pathtofile = self.censorc_path
        else:
            pathtofile = os.path.join(self.cwd, "censo.inp")

        with open(pathtofile, "w", newline=None) as outdata:
            outdata.write("$CENSO global configuration file: .censorc\n")
            outdata.write(f"$VERSION:{__version__} \n")
            outdata.write("\n")
            if usepaths:
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
                # outdata.write("cosmothermversion: 16\n")
            else:
                # TODO - why is this set up like that (including/excluding binary)??
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
                # outdata.write("cosmothermversion: 16\n")
            outdata.write("$ENDPROGRAMS\n\n")
            outdata.write("$CRE SORTING SETTINGS:\n")
            outdata.write("$GENERAL SETTINGS:\n")

            # FIXME - fix copy/paste code, fix function calls
            for key in OrderedDict(self.internal_settings.defaults_refine_ensemble_general):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.internal_settings.defaults_refine_ensemble_general)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                if key == "nconf" and value is None:
                    value = "all"
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART0 - CHEAP-PRESCREENING - SETTINGS:\n")
            for key in OrderedDict(self.internal_settings.defaults_refine_ensemble_part0):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.internal_settings.defaults_refine_ensemble_part0)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART1 - PRESCREENING - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.internal_settings.defaults_refine_ensemble_part1):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.internal_settings.defaults_refine_ensemble_part1)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART2 - OPTIMIZATION - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.internal_settings.defaults_refine_ensemble_part2):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.internal_settings.defaults_refine_ensemble_part2)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART3 - REFINEMENT - SETTINGS:\n")
            for key in OrderedDict(self.internal_settings.defaults_refine_ensemble_part3):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.internal_settings.defaults_refine_ensemble_part3)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$NMR PROPERTY SETTINGS:\n")
            outdata.write("$PART4 SETTINGS:\n")
            for key in OrderedDict(self.internal_settings.defaults_nmrprop_part4):
                value = self._exchange_onoff(
                    data.get(
                        key, OrderedDict(self.internal_settings.defaults_nmrprop_part4)[key]["default"]
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$OPTICAL ROTATION PROPERTY SETTINGS:\n")
            outdata.write("$PART5 SETTINGS:\n")
            for key in OrderedDict(self.internal_settings.defaults_optical_rotation_part5):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.internal_settings.defaults_optical_rotation_part5)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("$END CENSORC\n")

    def write_censo_inp(self, path=None):
        """
         Write file "censo.inp" which stores current settings in the .censorc format
         """
        if path is None:
            path = self.cwd
        args_key = {v: k for k, v in self.key_args_dict.items()}
        data = self.provide_runinfo(extend=False)
        with open(os.path.join(path, "censo.inp"), "w", newline=None) as outdata:
            outdata.write("$File: censo.inp settings of current calculation\n")
            outdata.write(f"$VERSION:{__version__} \n")
            outdata.write("\n")
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
            # outdata.write("cosmothermversion: 16\n")
            outdata.write("$ENDPROGRAMS\n\n")
            outdata.write("$CRE SORTING SETTINGS:\n")
            outdata.write("$GENERAL SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_general):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_general)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                if key == "nconf" and value is None:
                    value = "all"
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART0 - CHEAP-PRESCREENING - SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part0):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part0)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART1 - PRESCREENING - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part1):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part1)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART2 - OPTIMIZATION - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part2):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part2)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART3 - REFINEMENT - SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part3):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part3)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$NMR PROPERTY SETTINGS:\n")
            outdata.write("$PART4 SETTINGS:\n")
            for key in OrderedDict(self.defaults_nmrprop_part4):
                value = self._exchange_onoff(
                    data.get(
                        key, OrderedDict(self.defaults_nmrprop_part4)[key]["default"]
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$OPTICAL ROTATION PROPERTY SETTINGS:\n")
            outdata.write("$PART5 SETTINGS:\n")
            for key in OrderedDict(self.defaults_optical_rotation_part5):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_optical_rotation_part5)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$END censo.inp\n")


    def read_json(self, path, silent=False):
        """
        Reading stored data on conformers and information on settings of
        previous run.
        """
        if os.path.isfile(path):
            if not silent:
                print("Reading file: {}\n".format(os.path.basename(path)))
            try:
                with open(path, "r", encoding=CODING, newline=None) as inp:
                    save_data = json.load(inp, object_pairs_hook=OrderedDict)
            except (ValueError, TypeError, FileNotFoundError):
                print(f"{'ERROR:':{WARNLEN}}Your Jsonfile (enso.json) is corrupted!\n")
                time.sleep(0.02)
                raise
        return save_data


    def write_json(self, conformers, outfile="enso.json", path=None):
        """
        Dump conformer data and settings information of current run to json file
        """
        if not path:
            path = self.cwd
        
        data = {}
        data["settings_current"] = self.internal_settings.settings_current()
        
        data["ensemble"] = {}        
        for conf in conformers:
            data["ensemble"][conf.id] = conf
            
        with open(os.path.join(path, outfile), "w") as out:
            json.dump(data, out, indent=4, sort_keys=False)


    def read_program_paths(self, silent=False):
        """
        Get absolute paths of external programs employed in censo
        Read from the configuration file .censorc
        """
        # FIXME
        if silent:
            keep = self.save_errors

        # TODO - fix this with readrcfile decorator
        with open(self.censorc_path, "r") as inp:
            stor = inp.readlines()
        for line in stor:
            if "ctd =" in line:
                try:
                    self.external_paths["cosmorssetup"] = str(line.rstrip(os.linesep))
                except Exception:
                    self.save_errors.append(
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
                    self.save_errors.append(e)
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read settings for COSMO-RS from "
                        f".censorc!\n{'':{WARNLEN}}Most probably there is a user "
                        "input error."
                    )
            if "ORCA:" in line:
                try:
                    self.external_paths["orcapath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
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
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read ORCA version from .censorc!"
                    )
            if "GFN-xTB:" in line:
                try:
                    self.external_paths["xtbpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for GFNn-xTB from .censorc!"
                    )
                    if shutil.which("xtb") is not None:
                        self.external_paths["xtbpath"] = shutil.which("xtb")
                        self.save_errors.append(
                            f"{'':{WARNLEN}}Going to use {self.external_paths['xtbpath']} instead."
                        )
            if "CREST:" in line:
                try:
                    self.external_paths["crestpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for CREST from .censorc!"
                    )
                    if shutil.which("crest") is not None:
                        self.external_paths["crestpath"] = shutil.which("crest")
                        self.save_errors.append(
                            f"{'':{WARNLEN}}Going to use {self.external_paths['crestpath']} instead."
                        )
            if "mpshift:" in line:
                try:
                    self.external_paths["mpshiftpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for mpshift from .censorc!"
                    )
            if "escf:" in line:
                try:
                    self.external_paths["escfpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for escf from .censorc!"
                    )
            if "$ENDPROGRAMS" in line:
                break

            if silent:
                self.save_errors = keep


    def cleanup_run(self, complete=False):
            """
            Delete all unneeded files.
            """

            # files containing these patterns are deleted
            to_delete = [
                "enso.json.", 
                "enso_ensemble_part1.xyz.", 
                "enso_ensemble_part2.xyz.", 
                "enso_ensemble_part3.xyz.",
                "a3mat.tmp",
                "a2mat.tmp",
                "amat.tmp",
            ]

            
            # remove conformer_rotamer_check folder if complete cleanup
            if complete:
                print("Cleaning up the directory from ALL unneeded files!")
                to_delete[0] = "enso.json"
                if os.path.isdir(os.path.join(self.cwd, "conformer_rotamer_check")):
                    print("Removing conformer_rotamer_check")
                    shutil.rmtree(os.path.join(self.cwd, "conformer_rotamer_check"))
            else:
                print("Cleaning up the directory from unneeded files!")

            print(f"Be aware that files in {self.cwd} and subdirectories with names containing the following substrings will be deleted:")
            for sub in to_delete:
                print(sub)

            print("Do you wish to continue?")
            print("Please type 'yes' or 'no':")

            ui = input()
            if ui.strip().lower() not in ["yes", "y"]:
                print("Aborting cleanup!")
                sys.exit(0)

            # iterate over files in cwd and subdirs recursively and remove them if to delete
            deleted = 0
            for subdir, dirs, files in os.walk(self.cwd):
                for file in files:
                    if any([str in file for str in to_delete]) and int(file.split(".")[2]) > 1:
                        print(f"Removing: {file}")
                        os.remove(os.path.join(subdir, file))
                        deleted += os.path.getsize(os.path.join(subdir, file))

            print(f"Removed {deleted / (1024 * 1024): .2f} MB")

            """ # deprecated
            files_in_cwd = [
                f for f in os.listdir(self.cwd) if os.path.isfile(os.path.join(self.cwd, f))
            ]
                
            # get files like amat.tmp from COSMO calculation (can be several Mb)
            files_in_sub_cwd = glob.glob(self.cwd + "/**/**/**/*mat.tmp", recursive=True)
            size = 0
            for tmpfile in files_in_sub_cwd:
                if any(x in tmpfile for x in ["a3mat.tmp", "a2mat.tmp", "amat.tmp"]):
                    if os.path.isfile(tmpfile):
                        size += os.path.getsize(tmpfile)
                        # print(f"Removing {tmpfile} {os.path.getsize(tmpfile)} byte")
                        os.remove(tmpfile)
            print(f"Removed {size/(1024*1024): .2f} Mb")
                # ask if CONF folders should be removed """