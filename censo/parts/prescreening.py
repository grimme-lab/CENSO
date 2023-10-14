"""

"""
import os
from typing import List
import csv

from censo.cfg import PLENGTH, DIGILEN, AU2KCAL
from censo.utilities import (
    print,
    timeit,
    format_data,
)
from censo.parts.part import CensoPart
from censo.core import CensoCore
from censo.settings import CensoSettings
from censo.parallel import ProcessHandler
from censo.datastructure import MoleculeData

class Prescreening(CensoPart):
    
    alt_name = "part0"
    
    def __init__(self, core: CensoCore, settings: CensoSettings):
        super().__init__(core, settings, "prescreening")
        

    @timeit
    def run(self) -> None:
        """
        first screening of the ensemble by doing single-point calculation on the input geometries,
        using a (cheap) DFT method. if the ensemble optimization is not taking place in the gas-phase,
        the gsolv contribution is calculated using xtb.

        the list of conformers is then updated using Gtot (only DFT single-point energy if in gas-phase).
        """
        
        # initialize process handler for current program with conformer geometries
        handler = ProcessHandler(self._instructions, [conf.geom for conf in self.core.conformers])
        
        # set jobtype to pass to handler
        jobtype: List[str] = []
        if not self._instructions.get("gas-phase", None):
            if self._instructions.get("implicit", None):
                jobtype = ["sp", "gsolv"]
            else:
                jobtype = ["sp", "xtb_gsolv"]
        else:
            jobtype = ["sp"]
        
        # print instructions
        self.print_info()
        
        # set folder to do the calculations in for prescreening
        folder = os.path.join(self.core.workdir, self.__class__.__name__.lower())
        if os.path.isdir(folder):
            # TODO - warning? stderr?
            print(f"Folder {folder} already exists. Potentially overwriting files.")
        elif os.system(f"mkdir {folder}") != 0 and not os.path.isdir(folder):
            # TODO - error handling stderr case? is this reasonable?
            raise RuntimeError(f"Could not create directory for {self.__class__.__name__.lower()}")
        
        # compute results
        # for structure of results from handler.execute look there
        results = handler.execute(jobtype, folder)

        # update results for each conformer
        for conf in self.core.conformers:
            # if 'copy_mo' is enabled, get the mo_path from the results and store it in the respective GeometryData object
            if self._instructions.get("copy_mo", None):
                conf.geom.mo_path = results[id(conf)]["sp"]["mo_path"]

            # store results of single jobs for each conformer
            conf.results.setdefault(self.__class__.__name__.lower(), {}).update(results[id(conf)])
            
            # calculate free enthalpy values for every conformer
            conf.results[self.__class__.__name__.lower()]["gtot"] = self.key(conf)
        
        # sort conformers list with prescreening key (gtot)
        self.core.conformers.sort(
            key=lambda conf: conf.results[self.__class__.__name__.lower()]["gtot"],
        )  
        
        # calculate boltzmann weights from gtot values calculated here
        # trying to get temperature from instructions, set it to room temperature if that fails for some reason
        self.core.calc_boltzmannweights(
            self._instructions.get("temperature", 298.15),
            self.__class__.__name__.lower()
        )
            
        # write results (analogous to deprecated print)
        self.write_results()
                
        # update conformers with threshold
        threshold = self._instructions.get("threshold", None)
        if not threshold is None:
            # pick the free enthalpy of the lowest conformer
            limit = min([conf.results[self.__class__.__name__.lower()]["gtot"] for conf in self.core.conformers])
            
            # filter out all conformers above threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.key(x) > limit + threshold, 
                    self.core.conformers
                )
            ]
            
            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(filtered)  

            # TODO - print out which conformers are no longer considered
        else:
            """
            TODO
            print warning that no threshold could be determined
            (should not happen but doesn't necessarily brake the program)
            """
            print("...")
           
        # DONE


    def key(self, conf: MoleculeData) -> float:
        """
        prescreening key for conformer sorting
        calculates Gtot = E (DFT) + Gsolv (xtb) or Gsolv (DFT) for a given conformer
        """
        
        # Gtot = E (DFT) + Gsolv (xtb) or Gsolv (DFT)
        gtot: float = conf.results[self.__class__.__name__.lower()]["sp"]["energy"]
        
        # Gsolv is only calculated if prescreening is not in the gas-phase
        if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys():
            gtot += conf.results[self.__class__.__name__.lower()]["xtb_gsolv"]["gsolv"]
        elif "gsolv" in conf.results[self.__class__.__name__.lower()].keys():
            gtot += conf.results[self.__class__.__name__.lower()]["gsolv"]["gsolv"]
        
        return gtot


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
            
        also writes data in easily digestible format
        """
        
        # column headers
        headers = [
            "CONF#",
            "E (xTB)",
            "ΔE (xTB)",
            "E (DFT)",
            "ΔGsolv (xTB)",
            "Gtot",
            "ΔE (DFT)",
            "δΔGsolv",
            "ΔGtot",
            "Boltzmann weight",
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
            f"\% at {self._instructions.get('temperature', 298.15)} K",
        ]
        
        # variables for printmap
        # minimal xtb single-point energy
        # TODO - where do prescreening and screening get xtb single-point from?
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
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas']:.6f}",
            "ΔE (xTB)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}",
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}",
            "ΔGsolv (xTB)": lambda conf: 
                f"{conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv']:.6f}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys()
                else "---",
            "Gtot": lambda conf: f"{self.key(conf):.6f}",
            "ΔE (DFT)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['sp']['energy'] - dftmin) * AU2KCAL:.2f}",
            "δΔGsolv": lambda conf: 
                f"{(conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv'] - gsolvmin) * AU2KCAL:.2f}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() 
                else "---",
            "ΔGtot": lambda conf: f"{(self.key(conf) - gtotmin) * AU2KCAL:.2f}",
            "Boltzmann weight": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['bmw'] * 100:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # list the averaged free enthalpy of the ensemble
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
        lines.append(f"{self._instructions.get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part0==\n")
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")
        
        # lines.append(f">>> END of {self.__class__.__name__} <<<".center(PLENGTH, " ") + "\n")
            
        # write everything to a file
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.out"), "w", newline=None) as outfile:
            outfile.writelines(lines)

        # additionally, write data in csv format
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.csv"), "w", newline=None) as outfile:
            writer = csv.DictWriter(outfile, headers, delimiter=" ")
            writer.writeheader()
            rows = [{header: printmap[header](conf) for header in headers} for conf in self.core.conformers]
            writer.writerows(rows)