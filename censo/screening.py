"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv implicitly and include the RRHO contribution.
"""
from typing import List
import os
from math import isclose

from censo.prescreening import Prescreening
from censo.part import CensoPart
from censo.utilities import print, timeit, format_data
from censo.parallel import ProcessHandler
from censo.core import CensoCore
from censo.settings import CensoSettings
from censo.datastructure import MoleculeData
from censo.cfg import (
    PLENGTH,
    AU2KCAL,
)


class Screening(Prescreening):
    
    alt_name = "part1"
    
    def __init__(self, core: CensoCore, settings: CensoSettings):
        CensoPart.__init__(self, core, settings, "screening")


    @timeit
    def run(self) -> None:
        """
        Advanced screening of the ensemble by doing single-point calculations on the input geometries,
        but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

        Basically consists of two parts:
            - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
            - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
        """
        # TODO - test this
        # TODO - maybe put folder/handler as instance variable such that it can be reused later instead of reinstantiating
        super().run()

        # PART (2)
        # TODO - overwrite 'gtot'?
        threshold = self._instructions.get("threshold", None)
        
        # folder should already exist if previous part didn't raise runtime error
        folder = os.path.join(self.core.workdir, self.__class__.__name__.lower())
        if self._instructions["evaluate_rrho"]:
            # initialize process handler for current program with conformer geometries
            handler = ProcessHandler(self._instructions, [conf.geom for conf in self.core.conformers])
            
            jobtype = ["xtb_rrho"]

            # append results to previous results
            results = handler.execute(jobtype, folder)
            for conf in self.core.conformers:
                # update results for each conformer
                conf.results[self.__class__.__name__.lower()].update(results[id(conf)])
                
                # calculate new gtot including RRHO contribution
                conf.results[self.__class__.__name__.lower()]["gtot"] = self.key2(conf)

            # sort conformers list
            self.core.conformers.sort(key=lambda conf: conf.results[self.__class__.__name__.lower()]["gtot"])

            # update conformers with threshold
            # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
            limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
            
            # filter out all conformers above threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.key2(x) > limit + threshold, 
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

            # TODO - write results for second part
            self.write_results2()
                
        # DONE


    def key2(self, conf: MoleculeData) -> float:
        """
        Calculate the total Gibbs free energy (Gtot) of a given molecule using DFT single-point and gsolv (if included) and RRHO contributions.
        
        Parameters:
            conf (MoleculeData): The MoleculeData object containing the necessary information for the calculation.
        
        Returns:
            float: The total Gibbs free energy (Gtot) of the molecule.
        """
        # Gtot = E(DFT) + Gsolv + Grrho
        # NOTE: key2 should only be called if evaluate_rrho is True
        return self.key(conf) + conf.results[self.__class__.__name__.lower()]["xtb_rrho"]["energy"]


    def write_results(self) -> None:
        """
        Overrides the write_results function of Prescreening.
        Write the results to a file in formatted way.
        writes (1):
            E (xtb),
            δE (xtb),
            E (DFT),
            δGsolv (DFT),
            Gtot,
            δGtot
            
        Generates NO csv file. All info is included in the file written in write_results2.
        """
        # PART (1) of writing
        # column headers
        headers = [
            "CONF#",
            "E (xTB)",
            "ΔE (xTB)",
            "E (DFT)",
            "ΔGsolv",
            "Gtot",
            "ΔGtot",
        ]
        
        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[kcal/mol]",
        ]
        
        # variables for printmap
        # minimal xtb single-point energy (taken from prescreening)
        # TODO - where do prescreening and screening get xtb single-point from? - if prescreening is done grab from there, otherwise xtb_sp should be run
        xtbmin = min(
            conf.results["prescreening"]['xtb_gsolv']['energy_xtb_gas'] 
            for conf in self.core.conformers
        )
        
        # minimal total free enthalpy (single-point and potentially gsolv)
        gtotmin = min(self.key(conf) for conf in self.core.conformers)
        
        # determines what to print for each conformer in each column
        # TODO - remaining float accuracies
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: f"{conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas']:.6f}", # TODO
            "ΔE (xTB)": lambda conf: f"{(conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}", # TODO
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}",
            "ΔGsolv": lambda conf: 
                f"{self.key(conf) - conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() or "gsolv" in conf.results[self.__class__.__name__.lower()].keys()
                else "---", 
            "Gtot": lambda conf: f"{self.key(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.key(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # write everything to a file
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.out"), "w", newline=None) as outfile:
            outfile.writelines(lines)
        
        
    def write_results2(self) -> None:
        """
        Additional write function in case RRHO is used.
        Write the results to a file in formatted way. This is appended to the first file.
        writes (2):
            G (xtb),
            δG (xtb),
            E (DFT),
            δGsolv (DFT),
            Grrho,
            Gtot,
            δGtot
        
        Also writes them into an easily digestible format.
        """
        # column headers
        headers = [
            "CONF#",
            "G (xTB)",
            "ΔG (xTB)",
            "E (DFT)",
            "ΔGsolv",
            "GmRRHO",
            "Gtot",
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
            "[Eh]",
            "[kcal/mol]",
        ]

        # minimal xtb energy from single-point (and mRRHO)
        gxtbmin = min(
            conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions["temperature"]] # TODO?
            if self._instructions["evaluate_rrho"] else conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas']
            for conf in self.core.conformers
        )

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            self.key2(conf)
            for conf in self.core.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "G (xTB)": lambda conf: f"{conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}", # TODO
            "ΔG (xTB)": lambda conf: f"{(conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']] - gxtbmin) * AU2KCAL:.2f}", # TODO
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}",
            "ΔGsolv": lambda conf: 
                f"{self.key(conf) - conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}"
                if not isclose(self.key(conf), conf.results[self.__class__.__name__.lower()]['sp']['energy'])
                else "---", 
            "GmRRHO": lambda conf: 
                f"{conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}"
                if self._instructions["evaluate_rrho"]
                else "---", 
            "Gtot": lambda conf: f"{self.key2(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.key2(conf) - gtotmin) * AU2KCAL:.2f}",
        }
        
        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)
        
        # append lines to already existing file
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.out"), "a", newline=None) as outfile:
            outfile.writelines(lines)