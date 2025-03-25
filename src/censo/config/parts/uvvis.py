from typing import Literal
from pydantic import model_validator, Field

from .base import BasePartConfig
from ...assets import FUNCTIONALS
from ...params import GfnVersion, OrcaSolvMod, QmProg


class UVVisConfig(BasePartConfig):
    """Config class for UVVis"""

    prog: Literal[QmProg.ORCA] = QmProg.ORCA
    func: str = "wb97x-d4"
    basis: str = "def2-TZVP"
    sm: OrcaSolvMod = OrcaSolvMod.SMD
    gfnv: GfnVersion = GfnVersion.GFN2
    nroots: int = Field(gt=0, default=20)
    run: bool = False
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
