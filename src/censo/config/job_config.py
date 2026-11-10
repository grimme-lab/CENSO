from pydantic import ValidationInfo, field_validator
from typing import Literal


from ..config.paths import PathsConfig
from .generic import GenericConfig
from ..params import (
    OrcaSolvMod,
    TmSolvMod,
    XtbSolvMod,
    GfnVersion,
    GridLevel,
)


class XTBJobConfig(GenericConfig):
    """Configuration for XTB jobs."""

    paths: PathsConfig
    """Paths to external programs."""

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
        """
        Validate that solvent and solvation model are provided if not gas-phase.

        :param v: The value to validate.
        :param info: Validation info.
        :return: The validated value.
        """
        if not info.data["gas_phase"] and v is None:
            raise ValueError(
                "Solvent and solvation model need to be provided if not doing a gas-phase calculation."
            )

        return v


class SPJobConfig(GenericConfig):
    """Configuration for single-point jobs."""

    paths: PathsConfig
    """Paths to external programs."""

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

    temperatures: list[float] | None = None
    """Temperatures in K for COSMO-RS evaluation, main temperature last (only used in COSMO-RS calculations)."""

    @field_validator("solvent", "sm")
    @classmethod
    def need_solvent_when_no_gp(cls, v: str | None, info: ValidationInfo):
        """
        Validate that solvent and solvation model are provided if not gas-phase.

        :param v: The value to validate.
        :param info: Validation info.
        :return: The validated value.
        """
        if not info.data["gas_phase"] and v is None:
            raise ValueError(
                "Solvent and solvation model need to be provided if not doing a gas-phase calculation."
            )
        return v


class OptJobConfig(SPJobConfig):
    """Configuration for optimization jobs."""

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
    """Configuration for XTB optimization jobs."""

    hlow: float
    """hlow setting to be used in xtb (refer to xtb documentation)."""

    constraints: str | None = None
    """Constraints in xtb format as string."""


class RRHOJobConfig(XTBJobConfig):
    """Configuration for RRHO jobs."""

    sthr: float

    imagthr: float

    # bhess: bool
    # """Whether to utilize a single-point Hessian calculation for frequencies. Otherwise geometry will be optimized (equivalent to `--bhess` vs `--ohess` for xtb)."""

    temperatures: list[float]
    """Temperatures in K for the thermostatistical evaluation, main temperature last (xtb reports G(T) of the last entry in xtb_enso.json)."""

    # rmsdbias: bool


class NMRJobConfig(SPJobConfig):
    """Configuration for NMR jobs."""

    # NMR
    couplings: bool

    shieldings: bool

    active_nuclei: str

    fc_only: bool

    ss_cutoff: float


class RotJobConfig(SPJobConfig):
    """Configuration for optical rotation jobs."""

    freq: list[float]


class UVVisJobConfig(SPJobConfig):
    """Configuration for UVVis jobs."""

    # UVVis
    nroots: int
