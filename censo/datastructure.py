from typing import Dict, List
from functools import reduce

from censo.cfg import BOHR2ANG

class GeometryData:
    """
    Geometry contains geometry information as well as identifier to match it to a MoleculeData object
    in order to keep the object small, since it has to be pickled for multiprocessing
    """

    def __init__(self, id: int, name: str, xyz):
        """
        takes an identifier and the geometry lines from the xyz-file as input
        """
        
        # identifier linking it to a MoleculeData object
        # NOTE: this is runtime specific since it is set via built-in id-function
        self.id: int = id

        # name of the linked MoleculeData
        self.name: str = name
        
        # dict with element symbols as keys and lists of three-item lists as values
        # the coordinates should be given in Angstrom
        # self.xyz = {"H": [[0.0, 0.0, 0.0], [...], ...], "C": [[0.0, 0.0, 0.0], ...], ...}
        self.xyz: Dict[str, List[List[float]]] = {}
        
        # set up xyz dict from the input lines
        for line in xyz:
            spl = [s.strip() for s in line.split()]
            element = spl[0].capitalize()
            tmp = spl[1:]
            if element not in self.xyz.keys():
                self.xyz[element] = []
                
            self.xyz[element].append([float(i) for i in tmp])

        # compute number of atoms
        self.nat: int = sum(len(i) for i in self.xyz.values())

        # stores path to the most recent DFT MO-file (.gbw file in orca) 
        # (FIXME - quick and dirty solution, it would seem to make more sense to put this into a MoleculeData object)
        self.mo_path: str = None

    
    def toorca(self) -> List:
        """
        method to convert the internal cartesian coordinates to a data format usable by the OrcaParser
        """
        coord = []
        for element, allcoords in self.xyz.items():
            for atom in allcoords:
                coord.append([element] + atom)
        
        return coord


    def tocoord(self, path: str) -> None:
        """
        method to convert the internal cartesian coordinates (self.xyz) to a coord file (for tm or xtb)
        """
        coord = []
        for element, allcoords in self.xyz.items():
            for atom in allcoords:
                atom = list(map(lambda x: float(x) / BOHR2ANG, atom))
                coord.append(reduce(lambda x, y: f"{x} {y}", atom + [f"{element}\n"]))

        # write coord file
        with open(path, "w") as file:
            file.write("$coord\n")
            file.writelines(coord)
            file.write("$end\n")


    def fromcoord(self, path: str) -> None:
        """
        method to convert the content of a coord file to cartesian coordinates for the 'xyz' attribute
        """
        with open(path, "r") as file:
            lines = file.readlines()

        tmp = {}
        for line in lines:
            if not line.startswith("$"):
                tmp.setdefault(line.split()[-1], []).append([float(x) * BOHR2ANG for x in line.split()[:-1]])
            elif line.startswith("$end"):
                break


    def toxyz(self, path: str) -> None:
        """
        method to convert self.xyz to a xyz-file
        """
        with open(path, "w") as file:
            file.writelines([
                f"{self.nat}\n",
                f"{self.name}\n",
            ])
            for element, allcoords in self.xyz.items():
                for atom in allcoords:
                    file.write(f"{element} {atom[0]:.10f} {atom[1]:.10f} {atom[2]:.10f}\n")


class MoleculeData:
    """
    MoleculeData contains identifier, a GeometryData object, 
    as well as the sorting keys

    The confomers' MoleculeDatas are set up in censo.core.CensoCore.setup_conformers
    """
    
    def __init__(self, name: str, xyz):
        """
        takes geometry lines from the xyz-file as input to pass it to the GeometryData constructor
        """
        
        # stores a name for printing and (limited) between-run comparisons
        self.name: str = name
        
        # stores the geometry info to have a small object to be used for multiprocessing
        self.geom: GeometryData = GeometryData(id(self), self.name, xyz)
        
        # stores the initial xtb energy from CREST (or whatever was used before)
        self.xtb_energy: float = None
        
        # stores the results of the calculations
        self.results = {}
        # should be structured like the following:
        # 'part': <results from part jobs/in-part-calculations>
        # => e.g. self.results["prescreening"]["gtot"] 
        #    would return the free enthalpy of the conformer calculated in prescreening
        #    (not calculated with an external program)
        #
        #    self.results["prescreening"]["sp"] 
        #    returns the 'result' of the DFT single point in prescreening
        #    (calculated by external program)
        #    to get the single-point energy: self.results["prescreening"]["sp"]["energy"] 
        #    (confer to the results for each jobtype)