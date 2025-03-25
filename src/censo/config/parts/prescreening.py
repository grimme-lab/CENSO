from pydantic import model_validator, Field

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
