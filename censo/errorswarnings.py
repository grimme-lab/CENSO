"""
Contains custom error types
"""

class LogicError(Exception):
    """
    error type for check logic
    compared to normal exception:
    - info about which settings are problematic
    - info about what the problem is
    """
    def __init__(self, setting: str, problem: str, fix: str, *args):
        super().__init__(args)
        self._setting = setting
        self._problem = problem
        self._fix = fix
        
    def __str__(self) -> str:
        return f"LogicError: Error for setting {self._setting}\n{self._problem}\n{self._fix}\n"
    
    
class LogicWarning(Warning):
    """
    warning for potentially fatal errors in user input
    """
    def __init__(self, setting: str, problem: str, hint: str, fatal: bool = False, *args):
        super().__init__(args)
        self._setting = setting
        self._problem = problem
        self._hint = hint
        self.fatal = fatal
        
    
    def __str__(self) -> str:
        return  f"""{'Fatal ' if self.fatal else ''}LogicWarning: Possibly invalid value for setting {self._setting}
              {self._problem}
              {self._hint}
              {'Trying to continue.' if not self.fatal else ''}\n"""
                
                
class JobFailedWarning(Warning):
    """
    warning for failed job, specifying reasons TODO
    """