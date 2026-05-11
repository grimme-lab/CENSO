from typing import Literal, Self
from pydantic import model_validator, Field
import warnings

from .base import BasePartConfig
from ...params import GfnVersion, QmProg, TmSolvMod, OrcaSolvMod
from ...assets import FUNCTIONALS


class OptimizationConfig(BasePartConfig):
    """Config class for Opimization"""

    prog: QmProg = QmProg.TM
    """Program that should be used for calculations."""

    func: str = "r2scan-3c"
    """Functional that should be used for calculations."""

    basis: str = "def2-mTZVPP"
    """Basis set that should be used for calculations."""

    sm: OrcaSolvMod | TmSolvMod = TmSolvMod.DCOSMORS
    """Solvation model that should be used for calculations.

    When cosmors or cosmors-fine is chosen, geometry optimizations will use dcosmors
    and a final gsolv calculation with the chosen COSMO-RS model will be performed
    after all geometries are done optimizing, for final reranking.
    """

    gsolv_included: bool = False
    """Whether to explicitly calculate Gsolv contributions
    (False means gsolv will be calculated, True means it is already included in energy).
    """

    gfnv: GfnVersion = GfnVersion.GFN2
    """GFN version to be used for xtb calculations."""

    optcycles: int = Field(gt=0, default=8)
    """Number of microcycles per macrocycle."""

    maxcyc: int = Field(gt=0, default=200)
    """Maximum number of microcycles."""

    optlevel: Literal[
        "crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"
    ] = "normal"
    """Optlevel (ANCOPT settings, mapped to approximately equivalent settings for ORCA native optimizer)."""

    threshold: float = Field(gt=0, default=3.0)
    """ΔGtot threshold."""

    gradthr: float = Field(gt=0, default=0.01)
    """Gradnorm threshold below the normal ΔGtot threshold will be applied."""

    hlow: float = Field(gt=0, default=0.01)
    """Passed to hlow keyword in xcontrol files."""

    macrocycles: bool = True
    """Whether to use macrocycle optimization."""

    constrain: bool = False
    """Whether to use constraints for geometry optimizations."""

    xtb_opt: bool = True
    """Whether to use the ANCOPT optimizer as driver."""

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
    def constraints_only_with_ancopt(self):
        """
        Validate that constraints are only used with ANCOPT.

        :return: The validated instance.
        """
        if self.constrain and not self.xtb_opt:
            raise ValueError(
                "Constraints can currently only be used with ANCOPT. Enable xtb_opt to use constraints."
            )
        return self

    @model_validator(mode="after")
    def cosmors_uses_dcosmors_for_opt(self) -> Self:
        """
        When cosmors or cosmors-fine is chosen as solvent model, geometry optimizations
        will use dcosmors. Warn the user and set gsolv_included to False so that a
        final gsolv calculation with the chosen COSMO-RS model is performed for reranking.

        :return: The validated instance.
        """
        if self.sm in [TmSolvMod.COSMORS, TmSolvMod.COSMORS_FINE]:
            warnings.warn(
                f"Solvation model {self.sm.value} selected for optimization. "
                f"Geometry optimizations will use dcosmors, with a final reranking "
                f"after all geometries are done optimizing using {self.sm.value}."
            )
            self.gsolv_included = False
        return self
