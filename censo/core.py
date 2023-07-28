"""
stores ensembledata and conformers
functionality for program setup
"""

from argparse import Namespace
import os
import sys
from typing import Callable, Dict, List
from multiprocessing import Lock
from math import exp

from censo.cfg import (
    CODING,
    WARNLEN,
    __version__,
    AU2J,
    KB
)
from censo.datastructure import MoleculeData
from censo.ensembledata import EnsembleData
from censo.utilities import (
    check_for_float,
    do_md5,
    t2x,
    print,
)


# TODO - how do the assets files get into ~/.censo_assets?
class CensoCore:
    """
    """
    
    def __init__(self, args: Namespace, cwd: str):
        """
        store cwd and args, setup blank information storages
        """

        # current working directory
        self.cwd: str = cwd

        # parsed commandline arguments
        self.args: Nampespace = args
        
        # contains run-specific info that may change during runtime
        # initialized in CensoCore.read_input
        self.runinfo = {
            "nconf": int,
            "nat": int,
            "maxconf": int,
            "md5": str,
            "consider_unconverged": bool,
        }
            
        # stores the conformers with all info
        self.conformers: List[MoleculeData] = []
 
        # stores the conformers which were sorted out
        self.rem: List[MoleculeData] = []       
        
        # absolute path to ensemble input file
        self.ensemble_path: str

        # TODO - should this be done here?
        self.find_ensemble()

    
    def find_ensemble(self) -> None:
        """check for ensemble input file"""
        # if input file given via args use this path, otherwise set path to a default value
        if self.args.inp:
            ensemble_path = os.path.join(self.cwd, self.args.inp)
        else:
            ensemble_path = os.path.join(self.cwd, "crest_conformers.xyz")

        if os.path.isfile(ensemble_path):
            self.ensemble_path = ensemble_path
        else:
            """ print(
                f"{'ERROR:':{WARNLEN}}The input ensemble cannot be found!"
            ) """
            sys.exit(1)
            
        
    def read_input(self) -> None:
        """
        read ensemble input file (e.g. crest_conformers.xyz)
        """

        # store md5 hash for quick comparison of inputs later
        self.runinfo["md5"] = do_md5(self.ensemble_path)
        
        # if $coord in file => tm format, needs to be converted to xyz
        with open(self.ensemble_path, "r", encoding=CODING, newline=None) as inp:
            lines = inp.readlines()
            if any(["$coord" in line for line in lines]):
                _, self.runinfo["nat"], self.ensemble_path = t2x(
                        self.ensemble_path, writexyz=True, outfile="converted.xyz"
                    )
            else:
                self.runinfo["nat"] = int(lines[0].split()[0])
        
        # FIXME - temporary place for remaining settings
        # FIXME - where to put spearmanthr???
        if not self.args.spearmanthr:
            # set spearmanthr by number of atoms:
            self.spearmanthr = 1 / (exp(0.03 * (self.runinfo["nat"] ** (1 / 4))))

        self.runinfo["consider_unconverged"] = False
        # self.onlyread = False # FIXME - basically only used to print error???
        
        try:
            self.setup_conformers()
        except Exception as error: # TODO
            print(error.with_traceback)
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
            nat = self.runinfo["nat"]
            
            # check for correct line count in input 
            # assuming consecutive xyz-format coordinates
            if len(lines) % (nat + 2) == 0:
                if self.args.nconf:
                    nconf = int(min(self.args.nconf, len(lines) / (nat + 2)))
                    if self.args.nconf > nconf:
                        print(
                            f"{'WARNING:':{WARNLEN}}Given nconf is larger than max. number"
                            "of conformers in input file. Setting to the max. amount automatically."
                        )   
                else:
                    nconf = int(len(lines) / (nat + 2))
                
                self.runinfo["nconf"] = nconf
            else:
                raise Exception # TODO
            
            # get precalculated energies if possible
            for i in range(nconf):
                self.conformers.append(MoleculeData(f"CONF{i}", lines[i*nat+2:(i+1)*nat+2]))
                
                # precalculated energy set to 0.0 if it cannot be found
                self.conformers[i].xtb_energy = check_for_float(lines[i*nat+1]) or 0.0
            
            # also works if xtb_energy is None for some reason (None is put first)    
            self.conformers.sort(key=lambda x: x.xtb_energy, reverse=True)


    def update_conformers(self, filtered: List[MoleculeData]):
        """
        removing conformers from further consideration
        """

        # move the sorted out conformers to rem list
        for conf in filtered:
            # pop item from conformers and insert this item at index 0 in rem
            self.rem.insert(0, self.conformers.pop(self.conformers.index(conf)))


    def calc_boltzmannweights(self, temp: float, part: str) -> None:
        """
        Calculate weights for boltzmann distribution of ensemble at given temperature
        and given values for free enthalpy
        """
        # find lowest gtot value
        minfree: float = min([conf.results[part]["gtot"] for conf in self.conformers])
        
        # calculate boltzmann factors
        bmfactors = {
            id(conf):   conf.results[part]["gtot"] 
                        * exp(-(conf.results[part]["gtot"] - minfree) * AU2J / (KB * temp)) 
            for conf in self.conformers
        }
       
        # calculate boltzmann sum from bmfactors 
        bsum: float = sum(bmfactors.values())
        
        for conf in self.conformers:
            conf.results[part]["bmw"] = bmfactors[id(conf)] / bsum