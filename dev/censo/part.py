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
from censo.settings import CensoSettings, SettingsTuple


class CensoPart:
    
    def __init__(self):
        self.core: CensoCore = CensoCore()
        self.settings: SettingsTuple = CensoSettings().settings_current.get_settings(
            ,
            ,
            
        )
    
    
    @timeit
    def run(self) -> None:
        """
        setup jobs for all conformers
        run jobs via caller
        update ranking (decorator?)
        
        => basically the same for all parts, except for optrot, nmr, uvvis
           because ranking is not affected there
        
        difference between parts:
        job type
        settings for jobs
        print format?
        """
        settings = self.core.internal_settings.settings_current(parts=self.__class__.name)
        jobs = []
        
        for conformer in self.core.conformers:
            jobs.append(self.core.prog_job[settings[str]["prog"]]())
    
    
    def print(self):
        """
        print format for every part
        """
        pass