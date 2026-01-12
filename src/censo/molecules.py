from typing import cast
from dataclasses import dataclass
from pydantic import BaseModel
import warnings

from .params import BOHR2ANG, QmProg


@dataclass
class Contributions:
    """
    Energy contributions for a molecule.

    :ivar energy: Electronic energy contribution in Hartree.
    :vartype energy: float
    :ivar gsolv: Solvation free energy contribution in Hartree.
    :vartype gsolv: float
    :ivar grrho: Rigid-rotor-harmonic-oscillator (RRHO) free energy contribution in Hartree.
    :vartype grrho: float
    """

    energy: float = 0.0
    gsolv: float = 0.0
    grrho: float = 0.0


class Atom(BaseModel):
    """
    Represents an atom with element and coordinates.

    :ivar element: Chemical element symbol (e.g., 'C', 'H', 'O').
    :vartype element: str
    :ivar xyz: Cartesian coordinates (x, y, z) in Angstrom.
    :vartype xyz: tuple[float, float, float]
    """

    element: str
    xyz: tuple[float, float, float]


class GeometryData:
    """
    Container for molecular geometry information.

    Stores atomic positions and provides conversion methods between different
    coordinate formats (xyz, coord). Designed to be lightweight for efficient
    pickling in multiprocessing contexts.

    :ivar name: Identifier name of the geometry, used to match with MoleculeData objects.
    :vartype name: str
    :ivar xyz: List of Atom objects containing element symbols and Cartesian coordinates.
    :vartype xyz: list[Atom]
    :ivar nat: Number of atoms in the geometry.
    :vartype nat: int
    """

    def __init__(self, name: str, atoms: list[Atom]):
        """
        Initialize a GeometryData object with atoms.

        :param name: Name identifier for the geometry.
        :type name: str
        :param atoms: List of Atom objects. Also accepts list[str] for backwards
            compatibility (deprecated).
        :type atoms: list[Atom]
        :raises ValueError: If xyz line format is invalid when parsing from strings.
        """
        # name of the linked MoleculeData
        self.name: str = name

        self.xyz: list[Atom] = []

        if isinstance(atoms, list) and all(isinstance(x, Atom) for x in atoms):
            self.xyz = atoms
        elif isinstance(atoms, list) and all(isinstance(x, str) for x in atoms):
            # ensure backwards compatibility (would expect second arg to be list[str])
            warnings.warn(
                "Passing xyz file content as list of strings is deprecated, use GeometryData.from_xyz instead",
                DeprecationWarning,
            )
            # set up xyz dict from the input lines
            for line in atoms:
                spl = [s.strip() for s in cast(str, cast(object, line)).split()]
                if len(spl) != 4:
                    raise ValueError(f"Unexpected xyz line: {line}")
                element = spl[0].capitalize()
                x, y, z = (float(i) for i in spl[1:])
                self.xyz.append(Atom(element=element, xyz=(x, y, z)))

        # Count atoms
        self.nat: int = len(self.xyz)

    @classmethod
    def from_xyz(cls, name: str, xyz: list[str]):
        """
        Create a GeometryData object from xyz-file format lines.

        :param name: Name identifier for the geometry.
        :type name: str
        :param xyz: List of strings in xyz format (element x y z).
        :type xyz: list[str]
        :return: New GeometryData instance.
        :rtype: GeometryData
        :raises ValueError: If line format is invalid (not 4 space-separated values).
        """
        atoms = []
        # set up xyz dict from the input lines
        for line in xyz:
            spl = [s.strip() for s in line.split()]
            if len(spl) != 4:
                raise ValueError(f"Unexpected xyz line: {line}")
            element = spl[0].capitalize()
            x, y, z = (float(i) for i in spl[1:])
            atoms.append(Atom(element=element, xyz=(x, y, z)))
        return cls(name, atoms)

    @classmethod
    def from_mol(
        cls,
        name: str,
        atomic_numbers: list[int],
        coords: list[tuple[float, float, float]],
    ):
        """
        Create a GeometryData object from atomic numbers and coordinates.

        :param name: Name identifier for the geometry.
        :type name: str
        :param atomic_numbers: List of atomic numbers (e.g., [6, 1, 1] for CH2).
        :type atomic_numbers: list[int]
        :param coords: List of (x, y, z) coordinate tuples in Angstrom.
        :type coords: list[tuple[float, float, float]]
        :return: New GeometryData instance.
        :rtype: GeometryData
        """
        from .params import PSE

        atoms = [Atom(element=PSE[i], xyz=c) for i, c in zip(atomic_numbers, coords)]
        return cls(name, atoms)

    def toorca(self) -> list[str]:
        """
        Convert internal cartesian coordinates to ORCA format.

        :return: List of coordinate strings in ORCA format (element x y z).
        :rtype: list[str]
        """
        coord = []
        for atom in self.xyz:
            coord.append(" ".join([atom.element] + [str(c) for c in atom.xyz]))

        return coord

    def tocoord(self) -> list[str]:
        """
        Convert internal cartesian coordinates to TURBOMOLE coord file format.

        Coordinates are converted from Angstrom to Bohr for TURBOMOLE/XTB compatibility.

        :return: List of strings in coord file format, including $coord and $end markers.
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

    def fromcoord(self, *args):
        """
        Deprecated wrapper for update_from_coord_file.

        .. deprecated::
            Use update_from_coord_file instead.
        """
        warnings.warn(
            "fromcoord is deprecated, use update_from_coord_file instead",
            DeprecationWarning,
        )
        self.update_from_coord_file(*args)

    def update_from_coord_file(self, path: str) -> None:
        """
        Update geometry from a TURBOMOLE coord file.

        Reads a coord file and updates the xyz attribute with parsed coordinates.
        Coordinates are automatically converted from Bohr to Angstrom.

        :param path: Path to the coord file.
        :type path: str
        :raises ValueError: If coord line format is invalid.
        """
        with open(path) as file:
            lines = file.readlines()

        self.xyz = []
        for line in lines:
            if not line.startswith("$"):
                coords = line.split()
                if len(coords) != 4:
                    raise ValueError(f"Unexpected coord line: {line}")
                element = coords[-1]
                x, y, z = (float(x) * BOHR2ANG for x in coords[:-1])
                self.xyz.append(Atom(element=element, xyz=(x, y, z)))
            elif line.startswith("$end"):
                break

    def fromxyz(self, *args):
        """
        Deprecated wrapper for update_from_xyz_file.

        .. deprecated::
            Use update_from_xyz_file instead.
        """
        warnings.warn(
            "fromxyz is deprecated, use update_from_xyz_file instead",
            DeprecationWarning,
        )
        self.update_from_xyz_file(*args)

    def update_from_xyz_file(self, path: str) -> None:
        """
        Update geometry from an xyz file.

        Reads an xyz file and updates the xyz attribute with parsed coordinates.
        Skips the first two header lines (atom count and comment).

        :param path: Path to the xyz file.
        :type path: str
        """
        with open(path) as file:
            lines = file.readlines()

        self.xyz = []
        # Just skip the first two lines
        for line in lines[2:]:
            split = line.split()
            element = split[0]
            x, y, z = (float(x) for x in split[1:])
            self.xyz.append(Atom(element=element, xyz=(x, y, z)))

    def toxyz(self) -> list[str]:
        """
        Convert geometry to xyz-file format.

        :return: List of strings in xyz format:
            - Line 1: number of atoms
            - Line 2: molecule name
            - Remaining lines: element x y z (coordinates in Angstrom)
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
    Container for molecular data including geometry and energy contributions.

    Stores conformer information including geometry, charge, spin state, and
    energy contributions. Conformers are typically set up in
    censo.ensemble.EnsembleData.setup_conformers.

    :ivar name: Identifier name of the molecule/conformer.
    :vartype name: str
    :ivar geom: GeometryData object containing atomic coordinates.
    :vartype geom: GeometryData
    :ivar degen: Degeneracy factor of the conformer (default: 1).
    :vartype degen: int
    :ivar charge: Molecular charge.
    :vartype charge: int
    :ivar unpaired: Number of unpaired electrons.
    :vartype unpaired: int
    :ivar xtb_energy: Initial (biased) xTB energy from CREST or similar (optional).
    :vartype xtb_energy: float | None
    :ivar mo_paths: Dictionary storing paths to molecular orbital files by QM program.
    :vartype mo_paths: dict[str, list[str | tuple[str, str]]]
    """

    def __init__(
        self, name: str, geom: GeometryData, charge: int = 0, unpaired: int = 0
    ):
        """
        Initialize a MoleculeData object.

        :param name: Name identifier for the molecule.
        :type name: str
        :param geom: GeometryData object. Also accepts list[str] for backwards
            compatibility (deprecated).
        :type geom: GeometryData
        :param charge: Molecular charge (default: 0).
        :type charge: int
        :param unpaired: Number of unpaired electrons (default: 0).
        :type unpaired: int
        """

        # stores a name for printing and (limited) between-run comparisons
        self.name: str = name

        self.geom: GeometryData
        if isinstance(geom, GeometryData):
            # stores the geometry info to have a small object to be used for multiprocessing
            self.geom = geom
        elif isinstance(geom, list) and all(isinstance(x, str) for x in geom):
            # ensure backwards compatibility (would expect second arg to be list[str])
            warnings.warn(
                "Passing xyz file content as list of strings is deprecated, use MoleculeData.from_xyz instead",
                DeprecationWarning,
            )
            self.geom = GeometryData.from_xyz(name, geom)

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

    @classmethod
    def from_mol(
        cls,
        name: str,
        atomic_numbers: list[int],
        coords: list[tuple[float, float, float]],
        charge: int = 0,
        unpaired: int = 0,
    ):
        """
        Create a MoleculeData object from atomic numbers and coordinates.

        :param name: Name identifier for the molecule.
        :type name: str
        :param atomic_numbers: List of atomic numbers (e.g., [6, 1, 1] for CH2).
        :type atomic_numbers: list[int]
        :param coords: List of (x, y, z) coordinate tuples in Angstrom.
        :type coords: list[tuple[float, float, float]]
        :param charge: Molecular charge (default: 0).
        :type charge: int
        :param unpaired: Number of unpaired electrons (default: 0).
        :type unpaired: int
        :return: New MoleculeData instance.
        :rtype: MoleculeData
        """
        geom = GeometryData.from_mol(name, atomic_numbers, coords)
        return cls(name, geom, charge=charge, unpaired=unpaired)

    @classmethod
    def from_xyz(cls, name: str, xyz: list[str], charge: int = 0, unpaired: int = 0):
        """
        Create a MoleculeData object from xyz-file format lines.

        :param name: Name identifier for the molecule.
        :type name: str
        :param xyz: List of strings in xyz format (element x y z).
        :type xyz: list[str]
        :param charge: Molecular charge (default: 0).
        :type charge: int
        :param unpaired: Number of unpaired electrons (default: 0).
        :type unpaired: int
        :return: New MoleculeData instance.
        :rtype: MoleculeData
        """
        geom = GeometryData.from_xyz(name, xyz)
        return cls(name, geom, charge=charge, unpaired=unpaired)

    @property
    def energy(self) -> float:
        """
        Electronic energy of the molecule.

        :return: Electronic energy in Hartree.
        :rtype: float
        """
        return self.__energy

    @property
    def gsolv(self) -> float:
        """
        Solvation free energy contribution.

        :return: Solvation free energy (Gsolv) in Hartree.
        :rtype: float
        """
        return self.__gsolv

    @property
    def grrho(self) -> float:
        """
        Rigid-rotor-harmonic-oscillator (RRHO) free energy contribution.

        :return: RRHO free energy (Grrho) in Hartree.
        :rtype: float
        """
        return self.__grrho

    @property
    def gtot(self) -> float:
        """
        Total free energy (sum of all contributions).

        :return: Total free energy (energy + gsolv + grrho) in Hartree.
        :rtype: float
        """
        return self.__energy + self.__gsolv + self.__grrho

    def update(self, contributions: Contributions):
        """
        Update energy contributions for this molecule.

        :param contributions: Contributions object containing energy, gsolv, and grrho values.
        :type contributions: Contributions
        """
        self.__energy = contributions.energy
        self.__gsolv = contributions.gsolv
        self.__grrho = contributions.grrho
