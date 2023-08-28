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
    
    def __init__(self, core: CensoCore, settings: CensoSettings):
        self.core: CensoCore
        self.settings: CensoSettings
        
        # contains settings grabbed from CensoSettings instance, such as general settings etc.
        self._instructions: Dict[str, Any]

    
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