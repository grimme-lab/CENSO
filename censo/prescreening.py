"""

"""
import os
from typing import List
from censo.cfg import PLENGTH, DIGILEN, AU2KCAL
from censo.utilities import (
    print,
    timeit,
)
from censo.part import CensoPart
from censo.core import CensoCore
from censo.settings import CensoSettings
from censo.parallel import ProcessHandler
from censo.datastructure import MoleculeData

class Prescreening(CensoPart):
    
    def __init__(self, core: CensoCore, settings: CensoSettings):
        super().__init__(core, settings)
        
        self.core = core

        # grabs the settings required for this part from the passed 'CensoSettings' instance
        self.settings = settings
        self.run_settings = settings.settings_current.bypart("general") + settings.settings_current.bypart("prescreening")

        # transfers settings into a dict of instructions to be passed to the process handler
        self._instructions = {}
        for setting in self.run_settings:
            self._instructions[setting.name] = setting.value
    
    @timeit
    def run(self) -> None:
        """
        prescreening == part0, calculate cheap free energy on GFNn-xTB input geometry
        The idea is to improve on the description of E with a very fast DFT method.

        Consists of two parts:
        - DFT single-point
        - calculation of gsolv with xTB (optional, only if run not in gas-phase)
        """
        
        """
        1. setup job instructions (either gas phase sp or computation with solvent model)
        2. print instructions TODO - print only when 'verbose' (for debugging basically)?
        3. run calculations in parallel via helper
        4. print results
        """
        
        # initialize process handler for current program with conformer geometries
        handler = ProcessHandler(self.settings, [conf.geom for conf in self.core.conformers])
        
        # set jobtype to pass to handler
        jobtype: List[str] = []
        if self._instructions.get("gas-phase", None):
            jobtype = ["sp", "xtb_gsolv"]
        else:
            jobtype = ["sp"]
        
        # print instructions
        self.print_info()
        
        # set folder to do the calculations in for prescreening
        folder = os.path.join(self.core.cwd, 'prescreening')
        if os.path.isdir(folder):
            # TODO - warning? stderr?
            print(f"Folder {folder} already exists. Potentially overwriting files.")
        elif os.system(f"mkdir {folder}") != 0 and not os.path.isdir(folder):
            # TODO - stderr case?
            print(f"Could not create directory for {self.__class__.__name__.lower()}. Executing calculations in {self.core.cwd}.")
            folder = self.core.cwd
        
        # compute results
        # for structure of results from handler.execute look there
        results = handler.execute(jobtype, self._instructions, folder)

        # update results for each conformer
        for conf in self.core.conformers:
            # store results of single jobs for each conformer
            conf.results[self.__class__.__name__.lower()] = results[id(conf)]
            
            # calculate free enthalpy values for every conformer
            conf.results[self.__class__.__name__.lower()]["gtot"] = self.key(conf)
        
        # sort conformers list with prescreening key (gtot)
        self.core.conformers.sort(
            key=lambda conf: conf.results[self.__class__.__name__.lower()]["gtot"],
        )  
        
        # update conformers with threshold
        threshold = self._instructions.get("threshold", None)
        if not threshold is None:
            # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
            limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
            
            # filter out all conformers below threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.key(x) > limit + threshold, 
                    self.core.conformers
                )
            ]
            
            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(filtered)  
            
            # calculate boltzmann weights from gtot values calculated here
            # trying to get temperature from instructions, set it to room temperature if that fails for some reason
            self.core.calc_boltzmannweights(
                self._instructions.get("temperature", 298.15),
                self.__class__.__name__.lower()
            )
            
            # if no conformers are filtered basically nothing happens
        else:
            """
            TODO
            print warning that no threshold could be determined
            (should not happen but doesn't necessarily brake the program)
            """
            print("...")
        
        # write results (analogous to deprecated print)
        self.write_results()
                
        # DONE


    def key(self, conf: MoleculeData) -> float:
        """
        prescreening key for conformer sorting
        calculates Gtot = E (DFT) + Gsolv (xtb) for a given conformer
        """
        
        # Gtot = E (DFT) + Gsolv (xtb)
        gtot: float = conf.results[self.__class__.__name__.lower()]["sp"]["energy"]
        
        # Gsolv is only calculated if prescreening is not in the gas-phase
        if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys():
            gtot += conf.results[self.__class__.__name__.lower()]["xtb_gsolv"]["gsolv"]
        
        return gtot


    def print_info(self) -> None:
        """
        formatted write for prescreening
        write out settings used
        """
        
        lines = []
        lines.append("\n" + "".ljust(PLENGTH, "-") + "\n")
        lines.append("PRESCREENING - PART0".center(PLENGTH, " ") + "\n")
        lines.append("".ljust(PLENGTH, "-") + "\n")
        lines.append("\n")
        
        for instruction, val in self._instructions.items():
            lines.append(f"{instruction}:".ljust(DIGILEN, " ") + f"{val}\n")
            
        # print everything to console TODO - to stderr instead?
        for line in lines:
            print(line)


    def write_results(self) -> None:
        """
        writes: 
            E (xtb), 
            δE (xtb), 
            G_solv (xtb), 
            δG_solv,
            
            E(DFT), 
            δE(DFT), 
            
            E(DFT) + G_solv, 
            δ(E(DFT) + G_solv) 
            
        TODO also writes data in easily digestible format
        """
        # column headers
        headers = [
            "CONF",
            "E (xtb)",
            "ΔE (xtb)",
            "E (DFT)",
            "δGsolv (xtb)",
            "Gtot",
            "ΔE (DFT)",
            "ΔδGsolv",
            "ΔGtot",
        ]
        
        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
            "[kcal/mol]",
            "[kcal/mol]",
        ]
        
        # variables for printmap
        # minimal xtb single-point energy
        xtbmin = min(
            conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas'] 
            for conf in self.core.conformers
        )
        
        # minimal dft single-point energy
        dftmin = min(
            conf.results[self.__class__.__name__.lower()]['sp']['energy'] 
            for conf in self.core.conformers
        )
        
        # minimal solvation free enthalpy
        gsolvmin = min(
            conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv'] 
            if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() 
            else 0.0 
            for conf in self.core.conformers
        ) 
        
        # minimal total free enthalpy
        gtotmin = min(self.key(conf) for conf in self.core.conformers)
        
        # determines what to print for each conformer in each column
        # TODO - remaining float accuracies
        printmap = {
            "CONF": lambda conf: conf.name,
            "E (xtb)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas']}",
            "ΔE (xtb)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}",
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']}",
            "δGsolv (xtb)": lambda conf: 
                f"{conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv']}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys()
                else "---",
            "Gtot": lambda conf: f"{self.key(conf)}",
            "ΔE (DFT)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['sp']['energy'] - dftmin) * AU2KCAL:.2f}",
            "ΔδGsolv": lambda conf: 
                f"{(conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv'] - gsolvmin) * AU2KCAL:.2f}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() 
                else "---",
            "ΔGtot": lambda conf: f"{(self.key(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        # determine column width 'collen' of column with header 'header' 
        # by finding the length of the maximum width entry
        # for each column (header)
        collens = {
            header: collen for header, collen in zip(
                headers,
                (max(
                    len(header), max(len(printmap[header](conf)) 
                    for conf in self.core.conformers)
                ) for header in headers)
            )
        }
        
        lines = []
        
        # add table header
        # note: needs this amount of {}s because of one {} indicates fstring,
        #       two {} indicates to print {}, three to not print {collen},
        #       four {} to use collen within the header fstring variable
        lines.append(" ".join(f"{{header:^{collen}}}" for header, collen in collens.items()) + "\n")
        lines.append(" ".join(f"{{unit:^{collen}}}" for unit, collen in zip(units, collens.values())) + "\n")
        
        # add a row for every conformer (this time sorted by name)
        for conf in sorted(self.core.conformers, key=lambda conf: conf.name):
            # print floats with 2 digits accuracy if it is a difference, else with 6 digits
            lines.append(
                " ".join(
                        f"{{printmap[header](conf):^{collen}}}" 
                        if "Δ" in header 
                        else f"{{printmap[header](conf):^{collen}}}" 
                        for header, collen in collens.items()
                    ) 
                # draw an arrow if conformer is the best in current ranking
                + ("    <------\n" if self.key(conf) == self.key(self.core.conformers[0]) else "\n")
            )
            
        # list all conformers still considered with their boltzmann weights
        # as well as the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):\n"
        )
        lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n")
        print("".ljust(int(PLENGTH), "-") + "\n")

        # calculate averaged free enthalpy
        avG = sum([
            conf.results[self.__class__.__name__.lower()]["bmw"] 
            * conf.results[self.__class__.__name__.lower()]["gtot"] 
            for conf in self.core.conformers
        ])
        
        # calculate averaged free energy
        avE = sum([
            conf.results[self.__class__.__name__.lower()]["bmw"]
            * conf.results[self.__class__.__name__.lower()]["sp"]["energy"]
            for conf in self.core.conformers
        ])

        # append the lines for the free energy/enthalpy
        lines.append(f"{self._instructions.get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0==")
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")
        
        lines.append(">>> END of Prescreening <<<".center(PLENGTH, " ") + "\n")
            
        # write everything to a file
        with open(os.path.join(self.core.cwd, "prescreening.out"), "w") as outfile:
            outfile.writelines(lines)