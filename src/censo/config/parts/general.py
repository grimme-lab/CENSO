from pydantic import Field, model_validator


from .base import BasePartConfig
from ...params import XtbSolvMod
from ...assets import SOLVENTS


class GeneralConfig(BasePartConfig):
    """Config class for general settings"""

    temperature: float = Field(gt=0, default=298.15)
    """Temperature evaluated for mRRHO calculations."""

    # multitemp: bool = True

    # trange: tuple[float, float, float] = (273.15, 373.15, 5)

    evaluate_rrho: bool = True
    """Whether to run mRRHO calculations for GmRRHO contributions."""

    sm_rrho: XtbSolvMod = XtbSolvMod.GBSA
    """Solvation model to use for mRRHO calculations."""

    # consider_sym: bool = True

    # bhess: bool = True
    # """Wether to use single-point Hessians for mRRHO calculations."""

    # rmsdbias: bool = False

    imagthr: float = Field(lt=0, default=-100.0)
    """Highest acceptable imaginary mode (used for imagthr in xcontrol files)."""

    sthr: float = Field(default=50.0)
    """Value used for RR cutoff for mRRHO calculations using xtb."""

    # scale: float = Field(gt=0, default=1.0)

    solvent: str = "h2o"
    """Solvent to use for all calculations involving solvation."""

    gas_phase: bool = False
    """Whether to run all calculations in the gas-phase."""

    copy_mo: bool = True
    """Wether to copy MO files (.gbw for ORCA, mos/alpha, beta for TM) (can improve runtime significantly)."""

    balance: bool = True
    """Wether to use static load balancing to optimize processor usage."""

    ignore_failed: bool = False
    """Whether to ignore failed conformers or raise an error."""

    # @field_validator("trange")
    # @classmethod
    # def trange_must_be_valid(cls, v: tuple[float, float, float]):
    #     if v[0] < 0:
    #         raise ValueError("Negative temperature values not allowed (trange).")
    #     if v[1] < v[0]:
    #         raise ValueError(
    #             f"Order of arguments incorrect for trange: {v[1]} must be larger than {v[0]}."
    #         )
    #     if not v[2] > 0:
    #         raise ValueError(f"Step size for trange must be positive.")
    #
    #     return v
    #
    # @field_validator("trange", mode="before")
    # @classmethod
    # def trange_cast(cls, v: Any):
    #     if isinstance(v, str):
    #         parsed = ast.literal_eval(v)
    #         try:
    #             return tuple(parsed)
    #         except TypeError:
    #             return parsed
    #     return v

    @model_validator(mode="after")
    def solvent_must_be_valid_for_sm(self):
        """
        Validate that the solvent is available for the chosen solvation model.

        :return: The validated instance.
        """
        available_solvents = [
            solvent
            for solvent, keywords in SOLVENTS.items()
            if keywords.get(self.sm_rrho, None) is not None
        ]
        if self.solvent not in available_solvents:
            raise ValueError(f"Solvent {self.solvent} not defined for {self.sm_rrho}.")
        return self
