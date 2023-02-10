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
        
        # looks for censorc file (global configuration file)
        # if there is no rcfile, CENSO exits
        # looks for custom path and standard paths:
        # cwd and $home dir
        # absolute path to file directly
        self.censorc_path: str = self.find_rcfile()
        
        # if no input ensemble is found, CENSO exits
        # path has to be given via cml or the default path will be used:
        # "{cwd}/crest_conformers.xyz"
        # absolute path to file directly
        self.ensemble_path: str = self.find_ensemble()

        # TODO - remove (defined in cfg)
        # if no path is found, CENSO exits (assets are essential for functionality)
        # checks standard path first:
        # "~/.censo_assets"
        # TODO - add option for cml input
        # assets_path = path to folder
        self.assets_path: str = self.find_assets()
        
    
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

        censorc_name = ".censorc"

        tmp = [
            os.path.join(self.cwd, censorc_name),
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