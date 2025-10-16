from pydantic import model_validator, Field
from typing import Self

from .base import BasePartConfig
from ...params import GfnVersion, QmProg
from ...assets import FUNCTIONALS


class PrescreeningConfig(BasePartConfig):
    """Config class for Prescreening"""

    prog: QmProg = QmProg.TM
    """Program that should be used for calculations."""

    func: str = "pbe-d3"
    """Functional that should be used for calculations."""

    basis: str = "def2-SV(P)"
    """Basis set that should be used for calculations."""

    gfnv: GfnVersion = GfnVersion.GFN2
    """GFN version used for xtb calculations."""

    threshold: float = Field(gt=0, default=4.0)
    """ΔGtot threshold."""

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
