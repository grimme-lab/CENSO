from typing import Literal, Any, Self
import ast
from pydantic import model_validator, Field, field_validator

from .base import BasePartConfig
from ...assets import FUNCTIONALS
from ...params import QmProg


class RotConfig(BasePartConfig):
    """Config class for optical rotation calculations (gas-phase only)."""

    prog: Literal[QmProg.TM] = QmProg.TM
    """Program that should be used for calculations."""

    func: str = "pbe-d4"
    """Functional that should be used for calculations."""

    basis: str = "def2-SVPD"
    """Basis set that should be used for calculations."""

    freq: list[float] = Field(default_factory=lambda: [589.0, 633.0])
    """List of frequencies for which optical rotation will be calculated."""

    template: bool = False
    """Whether to use template files."""

    @field_validator("freq")
    @classmethod
    def freq_must_be_valid(cls, v: list[float]):
        """
        Validate that frequencies are positive and non-empty.

        :param v: The list of frequencies.
        :return: The validated list.
        """
        if len(v) == 0:
            raise ValueError("No frequencies provided.")
        if any(x <= 0 for x in v):
            raise ValueError("Negative or 0 frequencies not allowed.")

        return v

    @field_validator("freq", mode="before")
    @classmethod
    def freq_cast(cls, v: Any):
        """
        Cast string input to list of floats.

        :param v: The input value.
        :return: The casted list.
        """
        if isinstance(v, str):
            parsed = ast.literal_eval(v)
            if not isinstance(parsed, list) or not all(
                isinstance(x, float) or isinstance(x, int) for x in parsed
            ):
                raise ValueError(
                    f"No valid list for optical rotation frequencies provided ({v}). Must have correct list syntax and contain only numeric values."
                )
            return parsed
        return v

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
