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
from pprint import pprint

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
    
    def __init__(self, workdir: str, args: Namespace = None):
        """
        Setup a CensoCore object using the args from the command line
        workdir is the directory where the CENSO run should be executed in
        """

        # current working directory
        self.workdir: str = workdir

        # if args are given set accordingly, otherwise assume CENSO is used without commandline
        self.args: Nampespace = args
        
        # contains run-specific info that may change during runtime
        # initialized in CensoCore.read_input
        self.runinfo = {
            "nconf": int,
            "nat": int,
            "maxconf": int,
            "md5": str,
            "consider_unconverged": bool,
            "charge": int,
            "unpaired": int,
        }
            
        # stores the conformers with all info
        # NOTE: this is deliberately chosen to be a list since lists are ordered
        self.conformers: List[MoleculeData] = []
 
        # stores the conformers which were sorted out
        self.rem: List[MoleculeData] = []       
        
        # absolute path to ensemble input file
        self.ensemble_path: str

        
    def read_input(self, ensemble_path: str, charge: int = None, unpaired: int = None) -> None:
        """
        read ensemble input file (e.g. crest_conformers.xyz)
        """

        self.ensemble_path = ensemble_path

        # store md5 hash for quick comparison of inputs later
        self.runinfo["md5"] = do_md5(self.ensemble_path)
        
        # if $coord in file => tm format, needs to be converted to xyz
        with open(self.ensemble_path, "r") as inp:
            lines = inp.readlines()
            if any(["$coord" in line for line in lines]):
                _, self.runinfo["nat"], self.ensemble_path = t2x(
                        self.ensemble_path, writexyz=True, outfile="converted.xyz"
                    )
            else:
                self.runinfo["nat"] = int(lines[0].split()[0])

        # set charge and unpaired via funtion args or cml args
        if not self.args is None:
            self.runinfo["charge"] = charge or self.args.charge
            self.runinfo["unpaired"] = unpaired or self.args.unpaired
            
        if self.runinfo["charge"] is None or self.runinfo["unpaired"] is None:
            # TODO - error handling
            raise Exception("Charge or number of unpaired electrons not defined.") 

        # FIXME - temporary place for remaining settings
        # FIXME - where to put spearmanthr???
        if not self.args.spearmanthr:
            # set spearmanthr by number of atoms:
            self.spearmanthr = 1 / (exp(0.03 * (self.runinfo["nat"] ** (1 / 4))))

        self.runinfo["consider_unconverged"] = False
        # self.onlyread = False # FIXME - basically only used to print error???
        
        try:
            self.setup_conformers()
        except Exception as error: 
            raise error
        

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
                self.conformers.append(MoleculeData(f"CONF{i+1}", lines[2+i*(nat+2):(i+1)*(nat+2)]))
                
                # precalculated energy set to 0.0 if it cannot be found
                self.conformers[i].xtb_energy = check_for_float(lines[i*(nat+2)+1]) or 0.0
            
            # also works if xtb_energy is None for some reason (None is put first)    
            self.conformers.sort(key=lambda x: x.xtb_energy)


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