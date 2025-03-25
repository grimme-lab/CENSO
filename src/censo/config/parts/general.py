from pydantic import field_validator, Field, model_validator


from .base import BasePartConfig
from ...params import XtbSolvMod
from ...assets import SOLVENTS


class GeneralConfig(BasePartConfig):
    """Config class for general settings"""

    temperature: float = Field(gt=0, default=298.15)
    multitemp: bool = True
    trange: tuple[float, float, float] = (273.15, 373.15, 5)
    evaluate_rrho: bool = True
    sm_rrho: XtbSolvMod = XtbSolvMod.GBSA
    consider_sym: bool = True
    bhess: bool = True
    rmsdbias: bool = False
    imagthr: float = Field(lt=0, default=-100.0)
    sthr: float = Field(default=0.0)
    scale: float = Field(gt=0, default=1.0)
    solvent: str = "h2o"
    gas_phase: bool = False
    copy_mo: bool = True
    balance: bool = True

    @field_validator("trange")
    @classmethod
    def trange_must_be_valid(cls, v: tuple[float, float, float]):
        if v[0] < 0:
            raise ValueError("Negative temperature values not allowed (trange).")
        if v[1] < v[0]:
            raise ValueError(
                f"Order of arguments incorrect for trange: {v[1]} must be larger than {v[0]}."
            )
        if not v[2] > 0:
            raise ValueError(f"Step size for trange must be positive.")

        return v

    @model_validator(mode="after")
    def solvent_must_be_valid_for_sm(self):
        available_solvents = [
            solvent
            for solvent, keywords in SOLVENTS.items()
            if keywords.get(self.sm_rrho, None) is not None
        ]
        if self.solvent not in available_solvents:
            raise ValueError(f"Solvent {self.solvent} not defined for {self.sm_rrho}.")
        return self
