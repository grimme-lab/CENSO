from pydantic import Field, model_validator

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

    @model_validator(mode="after")
    def func_must_be_known_in_prog(self):
        prog: str = self.prog
        try:
            assert FUNCTIONALS[self.func][prog] is not None
            assert FUNCTIONALS[self.func]["disp"] is not None
            assert FUNCTIONALS[self.func]["type"] is not None
        except (KeyError, AssertionError):
            raise ValueError(
                f"Functional {self.func} not (fully) defined for prog {prog}."
            )
        return self
