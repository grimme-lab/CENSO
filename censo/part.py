"""
Part class as parent class for all parts of the calculation to
implement complete OOP approach.
"""

"""
- print format
- data created for each part (only the data needed for the part which is not global data
such as conformers or ensembledata)
- jobs
- main functions with @timeit
"""
from typing import Dict, Any

from censo.utilities import timeit
from censo.core import CensoCore
from censo.settings import CensoSettings, SettingsTuple


class CensoPart:
    
    def __init__(self, core: CensoCore, settings: CensoSettings, part: str):
        self.core: CensoCore = core
        self.settings: CensoSettings = settings
        
        # contains settings grabbed from CensoSettings instance, such as general settings etc.
        self._instructions: Dict[str, Any]

        # grabs the settings required for this part from the passed 'CensoSettings' instance
        paths = settings.settings_current.get("paths")
        general = settings.settings_current.get("general")
        specific = settings.settings_current.get(part)
        self._instructions = paths + general + specific

        # add some additional settings to _instructions so that the processors don't have to do any lookups
        # note: [1] auto-selects replacement solvent (TODO - print warning!)
        self._instructions["solvent_key_xtb"] = settings.solvents_db.get(self._instructions["solvent"])["xtb"][1]
        if 'sm' in self._instructions.keys():
            self._instructions[f"solvent_key_prog"] = settings.solvents_db.get(self._instructions["solvent"])[self._instructions["sm"]][1]    
        # TODO - doesn't work yet for parts where 'func' keyword doesn't exist or there are multiple functionals
        self._instructions["func_type"] = settings.dfa_settings.get_type(self._instructions["func"]) 

        # add 'charge' and 'unpaired' to instructions
        self._instructions["charge"] = core.runinfo.get("charge")
        self._instructions["unpaired"] = core.runinfo.get("unpaired")

        # set the correct name for 'func'
        self._instructions["func_name"] = settings.dfa_settings.get_name(self._instructions["func"], self._instructions["prog"])
        self._instructions["disp"] = settings.dfa_settings.get_disp(self._instructions["func"], self._instructions["prog"])
    
    @timeit
    def run(self) -> None:
        """
        what gets executed if the part is run
        should be implemented for every part respectively
        """
        pass    
    
    def key(self):
        """
        key to sort conformers list (optional)
        """
        pass
    
    
    def write_info(self):
        """
        formatted write of part instructions
        """
        pass
    
    
    def write_results(self):
        """
        formatted write of part results (optional)
        """
        pass