"""
Contains custom error types
"""

class LogicError(Exception):
    """
    error type for check logic
    compared to normal exception:
    - info about which settings are problematic
    - info about what the problem is
    - info about how to fix the error
    """
    def __init__(self, setting: str, problem: str, fix: str, *args):
        super().__init__(args)
        self._setting = setting
        self._problem = problem
        self._fix = fix
        
    def __str__(self) -> str:
        return f"LogicError: Invalid value for setting {self._setting}\n{self._problem}\n{self._fix}\n"