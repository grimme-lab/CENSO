from pydantic import model_validator, Field

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
