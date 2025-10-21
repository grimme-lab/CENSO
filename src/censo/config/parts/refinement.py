from typing import Self
from pydantic import model_validator, Field

from .base import BasePartConfig
from ...params import QmProg, TmSolvMod, OrcaSolvMod, GfnVersion
from ...assets import FUNCTIONALS
from ...logging import setup_logger

logger = setup_logger(__name__)


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

    @model_validator(mode="after")
    def gsolv_included_must_be_false_for_cosmors(self) -> Self:
        """
        In case the user chose cosmors or cosmors-fine as solvation model, it might be necessary to
        switch the gsolv_included flag to False.

        :return: The validated instance.
        """
        sm: TmSolvMod | OrcaSolvMod = self.sm
        if sm in [TmSolvMod.COSMORS, TmSolvMod.COSMORS_FINE] and self.gsolv_included:
            logger.warning(
                f"Found {sm.value} as solvation model but gsolv_included is set to True. Setting to False automaticaly."
            )
            self.gsolv_included = False

        return self
