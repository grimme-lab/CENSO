from pydantic import ValidationInfo, field_validator, Field

from .base import BasePartConfig
from ...params import GfnVersion, QmProg
from ...assets import FUNCTIONALS


class PrescreeningConfig(BasePartConfig):
    """Config class for Prescreening"""

    prog: QmProg = QmProg.TM
    func: str = "pbe-d4"
    basis: str = "def2-SV(P)"
    gfnv: GfnVersion = GfnVersion.GFN2
    threshold: float = Field(gt=0, default=4.0)
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
