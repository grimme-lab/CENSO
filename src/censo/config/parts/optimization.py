from typing import Literal
from pydantic import field_validator, ValidationInfo, Field

from .base import BasePartConfig
from ...params import GfnVersion, QmProg, TmSolvMod, OrcaSolvMod
from ...assets import FUNCTIONALS


class OptimizationConfig(BasePartConfig):
    """Config class for Opimization"""

    prog: QmProg = QmProg.TM
    func: str = "r2scan-3c"
    basis: str = "def2-mTZVPP"
    sm: OrcaSolvMod | Literal[TmSolvMod.DCOSMORS, TmSolvMod.COSMO] = TmSolvMod.DCOSMORS
    gfnv: GfnVersion = GfnVersion.GFN2
    optcycles: int = Field(gt=0, default=8)
    maxcyc: int = Field(gt=0, default=200)
    optlevel: Literal[
        "crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"
    ] = "normal"
    threshold: float = Field(gt=0, default=1.5)
    gradthr: float = Field(gt=0, default=0.01)
    hlow: float = Field(gt=0, default=0.01)
    macrocycles: bool = True
    constrain: bool = False
    xtb_opt: bool = True
    run: bool = True
    template: bool = False

    @field_validator("func")
    @classmethod
    def func_must_be_known_in_prog(cls, v: str, info: ValidationInfo):
        prog: str = info.data["prog"]
        try:
            assert FUNCTIONALS[v][prog] is not None
            assert FUNCTIONALS[v]["disp"] is not None
            assert FUNCTIONALS[v]["type"] is not None
        except (KeyError, AssertionError):
            raise ValueError(f"Functional {v} not (fully) defined for prog {prog}.")
        return v
