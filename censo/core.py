"""
stores ensembledata and conformers
functionality for program setup
"""

from argparse import Namespace
import os
import sys
from typing import Callable, Dict, List
from numpy import exp
from multiprocessing import Lock

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