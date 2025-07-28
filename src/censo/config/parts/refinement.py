from pydantic import model_validator, Field

from .base import BasePartConfig
from ...params import QmProg, TmSolvMod, OrcaSolvMod, GfnVersion
from ...assets import FUNCTIONALS


class RefinementConfig(BasePartConfig):
    """Config class for Refinement"""

    prog: QmProg = QmProg.TM
    """Program that should be used for calculations."""

    func: str = "wb97m-v"
    """Functional that should be used for calculations."""

    basis: str = "def2-TZVP"
    """Basis set that should be used for calculations."""

    sm: TmSolvMod | OrcaSolvMod = TmSolvMod.COSMORS
    """Solvation model that should be used for calculations."""

    gfnv: GfnVersion = GfnVersion.GFN2
    """GFN version to be used for xtb calculations."""

    threshold: float = Field(gt=0, le=1.0, default=0.95)
    """Boltzmann population threshold."""

    gsolv_included: bool = False
    """Whether to explicitly calculate Gsolv contributions."""

    template: bool = False
    """Whether to use template files."""

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
