from pydantic import field_validator, Field, ValidationInfo

from .base import BasePartConfig
from ...params import QmProg, TmSolvMod, OrcaSolvMod, GfnVersion
from ...assets import FUNCTIONALS


class RefinementConfig(BasePartConfig):
    """Config class for Refinement"""

    prog: QmProg = QmProg.TM
    func: str = "wb97x-v"
    basis: str = "def2-TZVP"
    sm: TmSolvMod | OrcaSolvMod = TmSolvMod.COSMORS
    gfnv: GfnVersion = GfnVersion.GFN2
    threshold: float = Field(gt=0, le=1.0, default=0.95)
    gsolv_included: bool = False
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
