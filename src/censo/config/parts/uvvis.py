from typing import Literal
from pydantic import ValidationInfo, field_validator, Field

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
