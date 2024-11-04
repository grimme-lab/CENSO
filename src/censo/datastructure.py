from functools import reduce
from typing import TypedDict

from .params import BOHR2ANG, Config


class Atom(TypedDict):
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
            self.xyz.append({"element": element, "xyz": [float(i) for i in spl[1:]]})

        # Count atoms
        self.nat: int = len(self.xyz)

    def toorca(self) -> list:
        """
        method to convert the internal cartesian coordinates to a data format usable by the OrcaParser
        """
        coord = []
        for atom in self.xyz:
            coord.append([atom["element"]] + atom["xyz"])

        return coord

    def tocoord(self) -> list[str]:
        """
        method to convert the internal cartesian coordinates (self.xyz) to coord file format (for tm or xtb)
        """
        coord = ["$coord\n"]
        for atom in self.xyz:
            coord.append(
                reduce(
                    lambda x, y: f"{x} {y}",
                    list(map(lambda x: float(x) / BOHR2ANG, atom["xyz"]))
                    + [f"{atom['element']}\n"],
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
                self.xyz.append({"element": element, "xyz": cartesian_coords})
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
            self.xyz.append({"element": element, "xyz": coords})

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
                f"{atom['element']} {atom['xyz'][0]:.10f} {atom['xyz'][1]:.10f} {atom['xyz'][2]:.10f}\n"
            )

        return lines


class MoleculeData:
    """
    The confomers' MoleculeData are set up in censo.ensembledata.EnsembleData.setup_conformers
    """

    def __init__(self, name: str, xyz: list[str]):
        """
        takes geometry lines from the xyz-file as input to pass it to the GeometryData constructor
        """

        # stores a name for printing and (limited) between-run comparisons
        self.name: str = name

        # stores the geometry info to have a small object to be used for multiprocessing
        self.geom: GeometryData = GeometryData(self.name, xyz)

        # stores the degeneration factor of the conformer
        self.degen: int = 1

        # stores the initial (biased) xtb energy from CREST (or whatever was used before)
        self.xtb_energy: float = None

        # list to store the paths to all MO-files from the jobs run for this conformer
        # might also include tuples if open shell and tm is used
        self.mo_paths: list[str, tuple] = []


class ParallelJob:

    def __init__(self, conf: GeometryData, jobtype: list[str]):
        # conformer for the job
        self.conf = conf

        # list of jobtypes to execute for the processor
        self.jobtype = jobtype

        # number of cores to use
        self.omp = Config.OMPMIN

        # stores path to an mo file which is supposed to be used as a guess
        # In case of open shell tm calculation this can be a tuple of files
        self.mo_guess = None

        # Stores all the important information for preparation of the input files for every jobtype
        # Always contains the 'general' key, which basically stores settings from the general section
        # that are supposed to be applied for every job
        # Also should always contain the name of the part where the job is launched from, as well as charge and
        # number of unpaired electrons
        # NOTE: prepinfo.keys() and items in jobtype are not necessarily the same! E.g. for NMR
        # jobtype = ["nmr"], prepinfo.keys() = ["nmr_s"], or prepinfo.keys() = ["nmr_s", "nmr_j"], ...
        self.prepinfo: dict[str, dict[str, any]] = {
            "general": {},
            "partname": "",
            "charge": 0,
            "unpaired": 0,
        }

        # store metadata, is updated by the processor
        # structure e.g.: {"sp": {"success": True, "error": None}, "xtb_rrho": {"success": False, ...}, ...}
        # always contains the "mo_path" key
        self.meta: dict[str, any] = {"mo_path": None}

        # store the results of the job
        self.results: dict[str, any] = {}

        # stores all flags for the jobtypes
        self.flags: dict[str, any] = {}
