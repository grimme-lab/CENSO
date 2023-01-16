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

from censo.utilities import timeit
from censo.core import CensoCore

class CensoPart:
    def __init__(self, core: CensoCore):
        self.core: CensoCore = core
    
    
    @timeit
    def run(self) -> None:
        pass
    
    
    def print(self):
        pass