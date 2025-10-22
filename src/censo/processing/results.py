from dataclasses import dataclass, field

from ..molecules import Atom


@dataclass
class Result:
    """Base class for results returned by calculations."""

    mo_path: str | tuple[str, str] = ""


@dataclass
class SPResult(Result):
    """Results class for single-point calculations."""

    energy: float = 0.0


@dataclass
class GsolvResult(Result):
    """Results class for Gsolv calculations."""

    gsolv: float = 0.0
    energy_gas: float = 0.0
    energy_solv: float = 0.0


@dataclass
class RRHOResult(Result):
    """Results class for mRRHO calculations."""

    energy: float = 0.0
    rmsd: float = 0.0
    gibbs: dict[float, float] = field(default_factory=dict)
    enthalpy: dict[float, float] = field(default_factory=dict)
    entropy: dict[float, float] = field(default_factory=dict)
    symmetry: str = "c1"
    symnum: int = 1
    linear: bool = False


@dataclass
class OptResult(Result):
    """Results class for geometry optimizations."""

    # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
    ecyc: list[float] = field(default_factory=list)
    # 'gncyc' contains the gradient norms for all cycles
    gncyc: list[float] = field(default_factory=list)
    # 'energy' contains the final energy of the optimization (converged or unconverged)
    energy: float = 0.0
    # 'geom' stores the optimized geometry in GeometryData.xyz format
    geom: list[Atom] = field(default_factory=list)
    grad_norm: float = 0.0
    cycles: int = 0
    converged: bool = False


@dataclass
class NMRResult(Result):
    """Results class for NMR calculations."""

    # atom index mapped onto shielding value, starting from 0
    shieldings: list[tuple[int, float]] = field(default_factory=list)
    # two atom indices mapped onto coupling value, starting from 0
    couplings: list[tuple[tuple[int, int], float]] = field(default_factory=list)


@dataclass
class RotResult(Result):
    """Results class for optical rotation calculations."""

    rotations_velocity: list[tuple[float, float]] = field(default_factory=list)
    rotations_length: list[tuple[float, float]] = field(default_factory=list)


@dataclass
class UVVisResult(Result):
    """Results class for UVVis calculations."""

    excitations: list[dict[str, float]] = field(default_factory=list)


@dataclass
class MetaData:
    """Metadata for job execution results."""

    conf_name: str
    success: bool = False
    error: str = ""
