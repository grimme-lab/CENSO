from typing import override, Self
import enum
import warnings

from pydantic import model_validator


from ...assets import FUNCTIONALS
from ...params import PLENGTH, QmProg
from ..generic import GenericConfig
from ...utilities import h2


class BasePartConfig(GenericConfig):
    """Base class for configuration classes for all parts"""

    @override
    def __str__(self) -> str:
        lines: list[str] = []
        lines.append(h2(f"{self.__class__.__name__.split("Config")[0]}"))
        kv = [(k, v) for k, v in self]
        kv.sort(key=lambda x: type(x[1]).__name__)
        for name, value in kv:
            display_value = value.value if isinstance(value, enum.Enum) else value
            lines.append(f"{name:>{PLENGTH // 2 - 2}} : {display_value}")

        return str("\n".join(lines))

    @model_validator(mode="after")
    def tm_gcp_d4_bug_check(self) -> Self:
        """
        Check if the model has basis, func, and prog fields and get their values. Backcheck with GCP basis sets.

        :return: The validated instance.
        """
        if all(s in self.__class__.model_fields for s in ["func", "basis", "prog"]):
            func: str = getattr(self, "func")

            try:
                disp = FUNCTIONALS[func]["disp"]
            except KeyError:
                raise ValueError("Received invalid functional key.")

            gcp_keywords = {
                "minis": "MINIS",
                "sv": "SV",
                "6-31g(d)": "631GD",
                "def2-sv(p)": "SV(P)",
                "def2-svp": "SVP",
                "def2-tzvp": "TZ",
            }
            if (
                getattr(self, "basis") in gcp_keywords
                and disp == "d4"
                and getattr(self, "prog") == QmProg.TM.value
            ):
                warnings.warn(
                    "Small basis set detected: Current version of TM includes a bug when combining GCP with DFT-D4. Switching to D3."
                )
                setattr(self, "func", func.replace("d4", "d3"))
                # TODO: there should be a way to validate this afterwards

        return self
