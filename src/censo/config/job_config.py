from pydantic import ValidationInfo, field_validator
from typing import Literal

from ..params import (
    GenericConfig,
    OrcaSolvMod,
    TmSolvMod,
    XtbSolvMod,
    GfnVersion,
    GridLevel,
)


class XTBJobConfig(GenericConfig):
    gfnv: GfnVersion
    """Version of GFN to use (GFN1-xTB, GFN2-xTB, GFN-FF)."""

    gas_phase: bool
    """Whether to ignore solvation related settings."""

    solvent: str | None = None
    """Solvent name to apply. Naming in solvation model will be automatically mapped."""

    sm_rrho: XtbSolvMod | None = None
    """Solvation model to use in xtb for the mRRHO calculation."""

    @field_validator("solvent", "sm_rrho")
    @classmethod
    def need_solvent_when_no_gp(cls, v: str | None, info: ValidationInfo):
        if not info.data["gas_phase"] and v is None:
            raise ValueError(
                "Solvent and solvation model need to be provided if not doing a gas-phase calculation."
            )

        return v


class SPJobConfig(GenericConfig):
    # DFT related
    copy_mo: bool
    """Whether to copy MO files for better guesses."""

    func: str
    """Functional name to use. Naming in program will be automatically mapped."""

    basis: str
    """Basis set to use."""

    grid: GridLevel
    """Grid level to be applied to this calculation."""

    template: bool
    """Whether to use template file for this calculation."""

    gas_phase: bool
    """Whether to ignore solvation related settings."""

    sm: OrcaSolvMod | TmSolvMod | None = None
    """Solvent model to use."""

    solvent: str | None = None
    """Solvent name to apply. Naming in solvation model will be automatically mapped."""

    # Only used for COSMORS calculations
    multitemp: bool | None = None
    """Equivalent to general settings variant. Whether to calculate temperature range (only used in COSMORS calculations)."""

    temperature: float | None = None
    """Equivalent to general settings variant. If not multitemp then use this temperature (only used in COSMORS calculations)."""

    trange: tuple[float, float, float] | None = None
    """Equivalent to general settings variant. (min, max, step) for the temperature range (only used in COSMORS calculations)."""

    @field_validator("solvent", "sm")
    @classmethod
    def need_solvent_when_no_gp(cls, v: str | None, info: ValidationInfo):
        if not info.data["gas_phase"] and v is None:
            raise ValueError(
                "Solvent and solvation model need to be provided if not doing a gas-phase calculation."
            )
        return v


class OptJobConfig(SPJobConfig):
    # Optimization
    macrocycles: bool
    """Whether to use macrocycle optimization."""

    optlevel: Literal[
        "crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"
    ]
    """Optlevel for the ANCOPT optimizer. Will get mapped to roughly equivalent levels in native optimizers."""

    optcycles: int
    """Steps per macrocycles if macrocycles are used."""

    maxcyc: int | None = None
    """Maximum number of optimization steps if doing a full opt."""


class XTBOptJobConfig(OptJobConfig):
    hlow: float
    """hlow setting to be used in xtb (refer to xtb documentation)."""


class RRHOJobConfig(XTBJobConfig):
    multitemp: bool
    """Equivalent to general settings variant. Whether to calculate temperature range (only used in COSMORS calculations)."""

    sthr: float

    imagthr: float

    bhess: bool
    """Whether to utilize a single-point Hessian calculation for frequencies. Otherwise geometry will be optimized (equivalent to `--bhess` vs `--ohess` for xtb)."""

    temperature: float

    trange: tuple[float, float, float]

    # rmsdbias: bool


class NMRJobConfig(SPJobConfig):
    # NMR
    couplings: bool

    shieldings: bool

    active_nuclei: str

    fc_only: bool

    ss_cutoff: float


class UVVisJobConfig(SPJobConfig):
    # UVVis
    nroots: int
