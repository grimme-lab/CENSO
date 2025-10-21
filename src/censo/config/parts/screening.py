from typing import Self
from pydantic import Field, model_validator

from .base import BasePartConfig
from ...params import OrcaSolvMod, QmProg, TmSolvMod, GfnVersion
from ...assets import FUNCTIONALS


class ScreeningConfig(BasePartConfig):
    """Config class for Screening"""

    prog: QmProg = QmProg.TM
    """Program that should be used for calculations."""

    func: str = "r2scan-3c"
    """Functional that should be used for calculations."""

    basis: str = "def2-mTZVPP"
    """Basis set that should be used for calculations."""

    sm: TmSolvMod | OrcaSolvMod = TmSolvMod.COSMORS
    """Solvation model that should be used for calculations."""

    gfnv: GfnVersion = GfnVersion.GFN2
    """GFN version to be used for xtb calculations."""

    threshold: float = Field(gt=0, default=3.5)
    """ΔGtot threshold."""

    gsolv_included: bool = True
    """Whether to explicitly calculate Gsolv contributions (False) or not (True)."""

    template: bool = False
    """Whether to use template files."""

    @model_validator(mode="after")
    def func_must_be_known_in_prog(self) -> Self:
        """
        Validate that the functional is known for the chosen program.

        :return: The validated instance.
        """
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
