from pydantic import BaseModel

from .params import BOHR2ANG, QmProg


class Atom(BaseModel):
    element: str
    xyz: list[float]


class GeometryData:
    """
    Geometry contains geometry information as well as identifier to match it to a MoleculeData object
    in order to keep the object small, since it has to be pickled for multiprocessing
    """

    def __init__(self, name: str, xyz: list[str]):
        """
        takes an identifier and the geometry lines from the xyz-file as input
        """
        # name of the linked MoleculeData
        self.name: str = name

        # list of dicts preserving the order of the input file for easy mapping
        # the coordinates should be given in Angstrom
        # self.xyz = [{"element": "H", "xyz": [0.0, 0.0, 0.0]}, {"element": "C", "xyz": [0.0, 0.0 0.7]}, ...]
        self.xyz: list[Atom] = []

        # set up xyz dict from the input lines
        for line in xyz:
            spl = [s.strip() for s in line.split()]
            element = spl[0].capitalize()
            self.xyz.append(Atom(element=element, xyz=[float(i) for i in spl[1:]]))

        # Count atoms
        self.nat: int = len(self.xyz)

    def toorca(self) -> list[str]:
        """
        method to convert the internal cartesian coordinates to a data format usable by the OrcaParser
        """
        coord = []
        for atom in self.xyz:
            coord.append(" ".join([atom.element] + [str(c) for c in atom.xyz]))

        return coord

    def tocoord(self) -> list[str]:
        """
        method to convert the internal cartesian coordinates (self.xyz) to coord file format (for tm or xtb)
        """
        coord = ["$coord\n"]
        for atom in self.xyz:
            coord.append(
                " ".join(
                    list(map(lambda x: str(float(x) / BOHR2ANG), atom.xyz))
                    + [f"{atom.element}\n"]
                )
            )

        coord.append("$end\n")

        return coord

    def fromcoord(self, path: str) -> None:
        """
        method to convert the content of a coord file to cartesian coordinates for the 'xyz' attribute
        """
        with open(path, "r") as file:
            lines = file.readlines()

        self.xyz = []
        for line in lines:
            if not line.startswith("$"):
                coords = line.split()
                element = coords[-1]
                cartesian_coords = [float(x) * BOHR2ANG for x in coords[:-1]]
                self.xyz.append(Atom(element=element, xyz=cartesian_coords))
            elif line.startswith("$end"):
                break

    def fromxyz(self, path: str) -> None:
        """
        Method to convert the content of an xyz file to cartesian coordinates for the 'xyz' attribute
        """
        with open(path, "r") as file:
            lines = file.readlines()

        self.xyz = []
        # Just skip the first two lines
        for line in lines[2:]:
            split = line.split()
            element = split[0]
            coords = [float(x) for x in split[1:]]
            self.xyz.append(Atom(element=element, xyz=coords))

    def toxyz(self) -> list[str]:
        """
        method to convert self.xyz to xyz-file format
        """
        lines = [
            f"{self.nat}\n",
            f"{self.name}\n",
        ]
        for atom in self.xyz:
            lines.append(
                f"{atom.element} {atom.xyz[0]:.10f} {atom.xyz[1]:.10f} {atom.xyz[2]:.10f}\n"
            )

        return lines


class MoleculeData:
    """
    The confomers' MoleculeData are set up in censo.ensembledata.EnsembleData.setup_conformers
    """

    def __init__(self, name: str, xyz: list[str], charge: int = 0, unpaired: int = 0):
        """
        takes geometry lines from the xyz-file as input to pass it to the GeometryData constructor
        """

        # stores a name for printing and (limited) between-run comparisons
        self.name: str = name

        # stores the geometry info to have a small object to be used for multiprocessing
        self.geom: GeometryData = GeometryData(self.name, xyz)

        # stores the degeneration factor of the conformer
        self.degen: int = 1

        self.charge: int = charge
        self.unpaired: int = unpaired

        # stores the initial (biased) xtb energy from CREST (or whatever was used before)
        self.xtb_energy: float | None = None

        # self.__energies: list[float] = []
        # self.__gsolvs: list[float] = []
        # self.__grrhos: list[float] = []
        self.energy: float = 0.0
        self.gsolv: float = 0.0
        self.grrho: float = 0.0

        # Stores Boltzmann populations in order of calculation
        # self.__bmws: list[float] = []
        self.bmw: float = 0.0

        # list to store the paths to all MO-files from the jobs run for this conformer
        # might also include tuples if open shell and tm is used
        self.mo_paths: dict[str, list[str | tuple[str, str]]] = {
            QmProg.ORCA: [],
            QmProg.TM: [],
        }

    # @property
    # def energy(self) -> float:
    #     """Most recent electronic gas-phase energy."""
    #     # TODO: is this really a good solution?
    #     # NOTE: this is the only time we don't check for the list length, since this really shouldn't be the case
    #     # when you want to access it
    #     return self.__energies[-1]
    #
    # @energy.setter
    # def energy(self, energy: float):
    #     """Set most recent gas-phase/implicitly solvated energy."""
    #     self.__energies.append(energy)
    #
    # @property
    # def gsolv(self) -> float:
    #     """Most recent solvation Gibbs free energy."""
    #     if len(self.__gsolvs) > 0:
    #         return self.__gsolvs[-1]
    #     else:
    #         return 0.0
    #
    # @gsolv.setter
    # def gsolv(self, gsolv: float):
    #     """Set most recent solvation Gibbs free energy."""
    #     self.__gsolvs.append(gsolv)
    #
    # @property
    # def grrho(self) -> float:
    #     """Most recent mRRHO contribution."""
    #     if len(self.__grrhos) > 0:
    #         return self.__grrhos[-1]
    #     else:
    #         return 0.0
    #
    # @grrho.setter
    # def grrho(self, grrho: float):
    #     """Set most recent mRRHO contribution."""
    #     self.__grrhos.append(grrho)
    #
    #
    # @property
    # def bmw(self) -> float:
    #     """Most recent Boltzmann population."""
    #     return self.__bmws[-1]
    #
    # @bmw.setter
    # def bmw(self, bmw: float):
    #     """Set most recent Boltzmann weight."""
    #     self.__bmws.append(bmw)
    @property
    def gtot(self) -> float:
        """Most recent energy+gsolv+grrho."""
        return self.energy + self.gsolv + self.grrho
