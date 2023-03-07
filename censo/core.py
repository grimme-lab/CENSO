"""
Wrapper to set up CENSO enviroment for runtime
store paths here
stores ensembledata and conformers
"""

import os
import sys
import json
import shutil
from argparse import Namespace
from collections import OrderedDict
from typing import Dict
import weakref
import functools
from numpy import exp

from censo.cfg import (
    CODING,
    WARNLEN,
    __version__,
)
from censo.orca_job import OrcaProc
from censo.settings import CensoSettings, PARTS, Settings
from censo.tm_job import TmProc
from censo.datastructure import MoleculeData
from censo.ensembledata import EnsembleData
from censo.utilities import (
    check_for_float,
    mkdir_p,
    do_md5,
    t2x,
    print,
)

# TODO - how do the assets files get into ~/.censo_assets?
class CensoCore:
    _instance_ref = weakref.WeakValueDictionary()
    
    prog_job: Dict[str, type] = {"tm": TmProc, "orca": OrcaProc}

    @staticmethod
    def factory(cwd, args: Namespace):
        """keep track of the CensoCore instance (should only be one)"""
        instance = CensoCore(cwd, args)
        CensoCore._instance_ref[id(instance)] = instance
        return instance


    @staticmethod
    def check_instance(f):
        """
        check number of CensoCore instances,
        exit if 0 or more than 1 to assure CENSO runs properly
        """
        @functools.wraps(f)
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
        # TODO - remove hardcoding here, also remove os.path.join uses with this
        self.censorc_name = ".censorc"
        
        # looks for censorc file (global configuration file)
        # if there is no rcfile, CENSO exits
        # looks for custom path and standard paths:
        # cwd and home dir
        # path to file directly
        self.censorc_path: str = self.find_rcfile()
        
        # if no path is found, CENSO exits (assets are essential for functionality)
        # checks standard path first:
        # "~/.censo_assets"
        # TODO - add option for cml input
        # path to folder
        self.assets_path: str = self.find_assets()
        
        # if no input ensemble is found, CENSO exits
        # path has to be given via cml or the default path will be used:
        # "{cwd}/crest_conformers.xyz"
        # path to file directly
        self.ensemble_path: str = self.find_ensemble()

        self.internal_settings = CensoSettings(self, self.args)

        # pathsdefaults: --> read_program_paths
        self.external_paths: Dict[str, str] = {}
        self.external_paths["orcapath"] = ""
        self.external_paths["orcaversion"] = ""
        self.external_paths["xtbpath"] = ""
        self.external_paths["crestpath"] = ""
        self.external_paths["cosmorssetup"] = ""
        self.external_paths["dbpath"] = ""
        self.external_paths["cosmothermversion"] = ""
        self.external_paths["mpshiftpath"] = ""
        self.external_paths["escfpath"] = ""
        
        self.conformers: list[MoleculeData] = []
        self.ensembledata: EnsembleData = EnsembleData()
        
        # run actions for which no complete setup is needed (FIXME - works?)
        if self.args.cleanup:
            self.cleanup_run()
            print("Removed files and going to exit!")
            sys.exit(0)
        elif self.args.cleanup_all:
            self.cleanup_run(True)
            print("Removed files and going to exit!")
            sys.exit(0)

        if self.args.writeconfig:
            self.write_config()
            sys.exit(0)
            
        # FIXME - temporary place for remaining settings
        if not self.args.spearmanthr:
            # set spearmanthr by number of atoms:
            self.spearmanthr = 1 / (exp(0.03 * (self.internal_settings.runinfo["nat"] ** (1 / 4))))

        self.internal_settings.runinfo["consider_unconverged"] = False
        

    """ def run_args(self) -> None:
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
        ### END solvent database adjustable by user
        # TODO - same here
        ### NMR reference shielding constant database adjustable by user
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
        ### END NMR reference shielding constant database adjustable by user

        # TODO - again similar thing
        ### ORCA user editable input
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
        ### END ORCA user editable input """

        
    def find_rcfile(self) -> str:
        """check for existing censorc"""

        tmp = [
            os.path.join(self.cwd, self.censorc_name),
            os.path.join(os.path.expanduser("~"), self.censorc_name)
        ]

        check = {
            os.path.isfile(tmp[0]): tmp[0],
            os.path.isfile(tmp[1]): tmp[1],
        }

        # FIXME - not the best solution to ensure code safety
        rcpath = ""

        # check for .censorc in standard locations if no path is given
        if not self.args.inprcpath:
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
        else:
            if os.path.isfile(self.args.inprcpath):
                rcpath = self.args.inprcpath
        
        # no censorc found at standard dest./given dest.
        if rcpath == "":
            print(
                f"No rcfile has been found. Do you want to create a new one?\n"
            )

            user_input = ""
            while user_input.strip().lower() not in ["yes", "y", "no", "n"]:
                print("Please type 'yes/y' or 'no/n':")
                user_input = input()
            
            if user_input.strip().lower() in ("y", "yes"):
                rcpath = self.write_config()
                self.censorc_name = "censorc_new"
            elif user_input.strip().lower() in ("n", "no"):
                print(
                    "Configuration file needed to run CENSO!\n"
                    "Going to exit!"
                )
                sys.exit(1)
            
        return rcpath

    
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


    def read_input(self) -> None:
        """
        read ensemble input file (e.g. crest_conformers.xyz)
        """
        
        # restart capability
        """ if self.args.restart and os.path.isfile(os.path.join(self.cwd, "enso.json")):
            tmp = config.read_json(os.path.join(self.cwd, "enso.json"), silent=True)
            if "ensemble_info" in tmp and self.args.inp is None:
                inpfile = os.path.basename(tmp["ensemble_info"].get("filename"))

                if os.path.isfile(inpfile):
                    self.args.inp = inpfile

                if self.args.debug:
                    print(f"Using Input file from: {inpfile}") """

        # store md5 hash for quick comparison of inputs later
        self.internal_settings.runinfo["md5"] = do_md5(self.ensemble_path)
        
        # if $coord in file => tm format, needs to be converted to xyz
        with open(self.ensemble_path, "r", encoding=CODING, newline=None) as inp:
            for line in inp:
                if "$coord" in line:
                    _, self.internal_settings.runinfo["nat"], self.ensemble_path = t2x(
                        self.ensemble_path, writexyz=True, outfile="converted.xyz"
                    )
                    break
        
        try:
            self.setup_conformers()
        except Exception as error: # TODO
            print(error)
            sys.exit(1)
        

    def setup_conformers(self) -> None:
        """
        open ensemble input
        split into conformers
        create MoleculeData objects out of coord input
        read out energy from xyz file if possible
        """
        # open ensemble input
        with open(self.ensemble_path, "r") as file:
            lines = file.readlines()
            nat = self.internal_settings.runinfo["nat"]
            
            # check for correct line count in input 
            # assuming consecutive xyz-format coordinates
            if len(lines) % (nat + 2):
                nconf = min(self.args.nconf, len(lines) / (nat + 2))
                if self.args.nconf > nconf:
                    print(
                        f"{'WARNING:':{WARNLEN}}Given nconf is larger than max. number"
                        "of conformers in input file. Setting to the max. amount automatically."
                    )
                    
                self.internal_settings.runinfo["nconf"] = nconf
            else:
                raise Exception # TODO
            
            # get precalculated energies if possible
            for i in range(nconf):
                self.conformers.append(MoleculeData(i, lines[i*nat+2:(i+1)*nat+2]))
                self.conformers[i].xtb_energy = check_for_float(lines[i*nat+1])
            
            # also works, if xtb_energy is None (None is put first)    
            self.conformers.sort(key=lambda x: x.xtb_energy)


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
        if not (self.censorc_path or self.censorc_path == ""):
            print(
                f"An existing configuration file ({self.censorc_name})" 
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
        if self.args.inprcpath:
            path = self.args.inprcpath
        # (path is never unbound)
        else:
            path = os.path.join(self.cwd, "censorc_new")
            
        print(f"Writing configuration file to {path} ...")
        with open(path, "w", newline=None) as outdata:
            outdata.write(f"$CENSO global configuration file: {self.censorc_name}\n")
            outdata.write(f"$VERSION:{__version__} \n")
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

            ### three loops to write settings to rcfile (TODO - make this into one?)
            # first, create headers and setup temporary storage of settings
            headers = []
            tmpsettings = {}
            for tmppart in PARTS:
                headers.append(f"${tmppart.upper()} SETTINGS\n")
                tmpsettings[tmppart] = {}
                
            # second, get default for each setting and store it
            for parts in CensoSettings.settings_options.values():
                for part, settings in parts.items():
                    if settings:
                        for setting, definition in settings.items():
                            tmpsettings[part][setting] = definition["default"]
              
            # third, write everything to outdata
            for header in headers:
                outdata.write(header)
                outdata.write("\n")
                for tmppart, tmpsett in tmpsettings.items():
                    for s, v in tmpsett.items():
                        outdata.write(f"{s}: {v}\n")
                        
            outdata.write("$END CENSORC\n")

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


    def read_config(self) -> Settings:
        """
        Read from config data from file (here enso.inp or .censorc),
        cml > .censorc > internal defaults
        every part has to be in rcfile
        """
        rcdata: Settings = {}
        with open(self.censorc_path, "r") as csvfile:
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
                    sett_type = CensoSettings.get_type(spl[0])
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


    """ def read_json(self, path, silent=False):
        Reading stored data on conformers and information on settings of
        previous run.
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
        Dump conformer data and settings information of current run to json file
        if not path:
            path = self.cwd
        
        data = {}
        data["settings_current"] = self.internal_settings.settings_current()
        
        data["ensemble"] = {}        
        for conf in conformers:
            data["ensemble"][conf.id] = conf
            
        with open(os.path.join(path, outfile), "w") as out:
            json.dump(data, out, indent=4, sort_keys=False)
 """

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