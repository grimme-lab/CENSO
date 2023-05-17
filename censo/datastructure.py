from typing import Dict, List


class GeometryData:
    """
    Geometry contains geometry information as well as identifier to match it to a MoleculeData object
    in order to keep the object small, since it has to be pickled for multiprocessing
    """

    def __init__(self, id: int, xyz):
        """
        takes an identifier and the geometry lines from the xyz-file as input
        """
        
        # identifier linking it to a MoleculeData object
        self.id: int = id
        
        # dict with element symbols as keys and lists of three-item lists as values
        self.xyz: Dict[str, List[List[float]]] = {}
        
        # set up xyz dict from the input lines
        for line in xyz:
            spl = [s.strip() for s in line.split(":")]
            element = spl[0].capitalize()
            tmp = spl[1:]
            if element not in self.xyz.keys():
                self.xyz[element] = []
                
            self.xyz[element].append([float(i) for i in tmp])
        
        
class MoleculeData:
    """
    MoleculeData contains identifier, a GeometryData object, 
    as well as the sorting keys
    """
    
    def __init__(self, name: str, xyz):
        """
        takes geometry lines from the xyz-file as input to pass it to the GeometryData constructor
        """
        
        # stores a name for printing and (limited) between-run comparisons
        self.name: str = name
        
        # stores the geometry info to have a small object to be used for multiprocessing
        self.geom: GeometryData = GeometryData(id(self), xyz)
        
        # stores the initial xtb energy from CREST (or whatever was used before)
        self.xtb_energy: float
        
        # stores the results of the calculations
        # should be structured like:
        # 'part': <results from part jobs/in-part-calculations>
        # => e.g. self.results["prescreening"]["gtot"] 
        #    would return the free enthalpy of the conformer calculated in prescreening
        #    (not calculated with an external program)
        #
        #    self.results["prescreening"]["sp"] 
        #    returns the energy calculated in the DFT single point in prescreening
        #    (calculated by external program)
        self.results = {}