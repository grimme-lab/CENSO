"""
stores all variables which have to be accessed by multiple classes
to avoid circular imports

also avoid to add functionality other than finding and storing information

these are primarily paths and the cml arguments
"""

from argparse import Namespace
import os
import shutil
import sys
from typing import Dict

class CensoStorage:
    def __init__(self, args: Namespace, cwd: str):
        self.args: Namespace = args
        self.cwd: str = cwd
        
        # TODO - remove hardcoding
        self.censorc_name: str = ".censorc"
        
        # looks for censorc file (global configuration file)
        # if there is no rcfile, CENSO exits
        # looks for custom path and standard paths:
        # cwd and home dir
        # path to file directly
        self.censorc_path: str = self.find_rcfile()
        
        # if no input ensemble is found, CENSO exits
        # path has to be given via cml or the default path will be used:
        # "{cwd}/crest_conformers.xyz"
        # path to file directly
        self.ensemble_path: str = self.find_ensemble()

        # if no path is found, CENSO exits (assets are essential for functionality)
        # checks standard path first:
        # "~/.censo_assets"
        # TODO - add option for cml input
        # assets_path = path to folder
        self.assets_path: str = self.find_assets()
        
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
            """ print(
                f"{'ERROR:':{WARNLEN}}The input ensemble cannot be found!"
            ) """
            sys.exit(1)
            
        
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
        
        return rcpath
    
    
    def find_assets(self) -> str:
        """
        look for assets folder
        if it is not found, CENSO exits
        """
        assets_path = os.path.expanduser("~/.censo_assets")
        if not os.path.isdir(assets_path):
            """ print(f"{'WARNING:':{WARNLEN}}The folder '~/.censo_assets/' designed for additional "
                "remote configuration files can not be found!\n"
            ) """
            sys.exit(1)

        return assets_path
    
    
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