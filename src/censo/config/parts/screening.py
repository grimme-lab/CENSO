from pydantic import Field, field_validator, ValidationInfo

from .base import BasePartConfig
from ...params import OrcaSolvMod, QmProg, TmSolvMod, GfnVersion
from ...assets import FUNCTIONALS


class ScreeningConfig(BasePartConfig):
    """Config class for Screening"""

    prog: QmProg = QmProg.TM
    func: str = "r2scan-3c"
    basis: str = "def2-mTZVPP"
    sm: TmSolvMod | OrcaSolvMod = TmSolvMod.COSMORS
    gfnv: GfnVersion = GfnVersion.GFN2
    threshold: float = Field(gt=0, default=3.5)
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
