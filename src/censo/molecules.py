from dataclasses import dataclass
from pydantic import BaseModel

from .params import BOHR2ANG, QmProg


@dataclass
class Contributions:
    """
    Energy contributions for a molecule.
    """

    energy: float = 0.0
    gsolv: float = 0.0
    grrho: float = 0.0


class Atom(BaseModel):
    """
    Represents an atom with element and coordinates.
    """

    element: str
    xyz: list[float]


class GeometryData:
    """
    Geometry contains geometry information as well as identifier to match it to a MoleculeData object
    in order to keep the object small, since it has to be pickled for multiprocessing
    """

    def __init__(self, name: str, xyz: list[str]):
        """
        Takes an identifier and the geometry lines from the xyz-file as input.

        :param name: Name of the geometry.
        :type name: str
        :param xyz: List of xyz lines.
        :type xyz: list[str]
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
        Method to convert the internal cartesian coordinates to a data format usable by the OrcaParser.

        :return: List of coordinate strings.
        :rtype: list[str]
        """
        coord = []
        for atom in self.xyz:
            coord.append(" ".join([atom.element] + [str(c) for c in atom.xyz]))

        return coord

    def tocoord(self) -> list[str]:
        """
        Method to convert the internal cartesian coordinates (self.xyz) to coord file format (for tm or xtb).

        :return: List of coord lines.
        :rtype: list[str]
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
        Method to convert the content of a coord file to cartesian coordinates for the 'xyz' attribute.

        :param path: Path to the coord file.
        :type path: str
        :return: None
        :rtype: None
        """
        with open(path) as file:
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
        Method to convert the content of an xyz file to cartesian coordinates for the 'xyz' attribute.

        :param path: Path to the xyz file.
        :type path: str
        :return: None
        :rtype: None
        """
        with open(path) as file:
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
        Method to convert self.xyz to xyz-file format.

        :return: List of xyz lines.
        :rtype: list[str]
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
    The confomers' MoleculeData are set up in censo.ensemble.EnsembleData.setup_conformers
    """

    def __init__(self, name: str, xyz: list[str], charge: int = 0, unpaired: int = 0):
        """
        Takes geometry lines from the xyz-file as input to pass it to the GeometryData constructor.

        :param name: Name of the molecule.
        :type name: str
        :param xyz: List of xyz lines.
        :type xyz: list[str]
        :param charge: Charge of the molecule.
        :type charge: int
        :param unpaired: Number of unpaired electrons.
        :type unpaired: int
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

        self.__energy: float = 0.0
        self.__gsolv: float = 0.0
        self.__grrho: float = 0.0

        # list to store the paths to all MO-files from the jobs run for this conformer
        # might also include tuples if open shell and tm is used
        self.mo_paths: dict[str, list[str | tuple[str, str]]] = {
            QmProg.ORCA: [],
            QmProg.TM: [],
        }

    @property
    def energy(self) -> float:
        """
        Current energy.

        :return: Energy value.
        :rtype: float
        """
        return self.__energy

    @property
    def gsolv(self) -> float:
        """
        Current gsolv.

        :return: Gsolv value.
        :rtype: float
        """
        return self.__gsolv

    @property
    def grrho(self) -> float:
        """
        Current grrho.

        :return: Grrho value.
        :rtype: float
        """
        return self.__grrho

    @property
    def gtot(self) -> float:
        """
        Current energy+gsolv+grrho.

        :return: Total free energy.
        :rtype: float
        """
        return self.__energy + self.__gsolv + self.__grrho

    def update(self, contributions: Contributions):
        """
        Update contributions.

        :param contributions: Contributions to update with.
        :type contributions: Contributions
        :return: None
        :rtype: None
        """
        self.__energy = contributions.energy
        self.__gsolv = contributions.gsolv
        self.__grrho = contributions.grrho
