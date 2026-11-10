import ast
from typing import Any, Self

from pydantic import Field, field_validator, model_validator


from .base import BasePartConfig
from ...params import XtbSolvMod, XTB_MAX_TEMPERATURES
from ...assets import SOLVENTS


class GeneralConfig(BasePartConfig):
    """Config class for general settings"""

    temperature: float = Field(gt=0, default=298.15)
    """Temperature evaluated for mRRHO calculations."""

    multitemp: bool = False
    """Whether to evaluate temperature-dependent contributions (mRRHO, COSMO-RS) for all temperatures in trange."""

    trange: tuple[float, float, float] = (273.15, 373.15, 5.0)
    """Temperature range (start, end, step) in K, end inclusive. Only used if multitemp is set."""

    evaluate_rrho: bool = True
    """Whether to run mRRHO calculations for GmRRHO contributions."""

    sm_rrho: XtbSolvMod = XtbSolvMod.GBSA
    """Solvation model to use for mRRHO calculations."""

    # consider_sym: bool = True

    # bhess: bool = True
    # """Wether to use single-point Hessians for mRRHO calculations."""

    # rmsdbias: bool = False

    imagthr: float = Field(lt=0, default=-100.0)
    """Highest acceptable imaginary mode (used for imagthr in xcontrol files)."""

    sthr: float = Field(default=50.0)
    """Value used for RR cutoff for mRRHO calculations using xtb."""

    # scale: float = Field(gt=0, default=1.0)

    solvent: str = "h2o"
    """Solvent to use for all calculations involving solvation."""

    gas_phase: bool = False
    """Whether to run all calculations in the gas-phase."""

    copy_mo: bool = True
    """Wether to copy MO files (.gbw for ORCA, mos/alpha, beta for TM) (can improve runtime significantly)."""

    balance: bool = True
    """Wether to use static load balancing to optimize processor usage."""

    ignore_failed: bool = False
    """Whether to ignore failed conformers or raise an error."""

    @field_validator("trange", mode="before")
    @classmethod
    def trange_cast(cls, v: Any):
        """
        Cast string input to a tuple of floats.

        :param v: The input value.
        :return: The casted tuple.
        """
        parsed = ast.literal_eval(v) if isinstance(v, str) else v
        if not isinstance(parsed, (list, tuple)) or not (
            len(parsed) == 3 and all(isinstance(x, (int, float)) for x in parsed)
        ):
            raise ValueError(
                f"No valid temperature range provided ({v}). Must be (start, end, step) with numeric values."
            )
        return tuple(parsed)

    @field_validator("trange")
    @classmethod
    def trange_must_be_valid(cls, v: tuple[float, float, float]):
        """
        Validate that the temperature range is positive and ordered with a positive step.

        :param v: The temperature range (start, end, step).
        :return: The validated temperature range.
        """
        if v[0] <= 0:
            raise ValueError("Temperatures in trange must be positive.")
        if v[1] < v[0]:
            raise ValueError(
                f"trange end ({v[1]}) must not be smaller than start ({v[0]})."
            )
        if v[2] <= 0:
            raise ValueError("Step size for trange must be positive.")
        return v

    def get_temperatures(self) -> list[float]:
        """
        Temperatures to evaluate (rounded to 0.01 K), ascending, main temperature always last.

        :return: List of temperatures in K.
        """
        main = round(self.temperature, 2)
        if not self.multitemp:
            return [main]
        start, end, step = self.trange
        nsteps = int((end - start) / step + 1e-9)
        grid = {round(start + i * step, 2) for i in range(nsteps + 1)}
        grid.discard(main)
        return sorted(grid) + [main]

    @model_validator(mode="after")
    def temperatures_must_fit_xtb(self) -> Self:
        """
        Validate that xtb can evaluate all temperatures of the grid in one run.

        :return: The validated instance.
        """
        n = len(self.get_temperatures())
        if n > XTB_MAX_TEMPERATURES:
            raise ValueError(
                f"multitemp evaluates {n} temperatures, but xtb accepts at most {XTB_MAX_TEMPERATURES}. Increase the trange step or narrow the range."
            )
        return self

    @model_validator(mode="after")
    def solvent_must_be_valid_for_sm(self):
        """
        Validate that the solvent is available for the chosen solvation model.

        :return: The validated instance.
        """
        available_solvents = [
            solvent
            for solvent, keywords in SOLVENTS.items()
            if keywords.get(self.sm_rrho, None) is not None
        ]
        if self.solvent not in available_solvents:
            raise ValueError(f"Solvent {self.solvent} not defined for {self.sm_rrho}.")
        return self
