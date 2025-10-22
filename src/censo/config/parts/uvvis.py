from typing import Literal, Self
from pydantic import model_validator, Field

from .base import BasePartConfig
from ...assets import FUNCTIONALS
from ...params import OrcaSolvMod, QmProg


class UVVisConfig(BasePartConfig):
    """Config class for UVVis"""

    prog: Literal[QmProg.ORCA] = QmProg.ORCA
    """Program that should be used for calculations."""

    func: str = "wb97x-d4"
    """Functional that should be used for calculations."""

    basis: str = "def2-TZVP"
    """Basis set that should be used for calculations."""

    sm: OrcaSolvMod = OrcaSolvMod.SMD
    """Solvation model that should be used for calculations."""

    nroots: int = Field(gt=0, default=20)
    """Number of roots to calculate."""

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
