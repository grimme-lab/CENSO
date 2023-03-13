"""
"""

from typing import Dict


class EnsembleData:        
    """
    contains metadata about the ensemble
    """
    
    def __init__(self):
        # stores the most current sorted list of conformers

        # stores the ids of the best conformers for each part
        self.bestconf: Dict[str, int]
        
        # stores the number of considered conformers for each part
        self.partsnconf: Dict[str, int]