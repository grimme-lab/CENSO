from typing import override, Any
import warnings

from pydantic import model_validator


from ...assets import FUNCTIONALS
from ...params import PLENGTH, GenericConfig, QmProg
from ...utilities import h2


class BasePartConfig(GenericConfig):
    """Base class for configuration classes for all parts"""

    @override
    def __str__(self):
        lines: list[str] = []
        lines.append(h2(f"{self.__class__.__name__.split("Config")[0]}"))
        kv = [(k, v) for k, v in self]
        kv.sort(key=lambda x: type(x[1]).__name__)
        for name, value in kv:
            lines.append(f"{name:>{PLENGTH // 2 - 2}} : {value}")

        return str("\n".join(lines))

    @model_validator(mode="before")
    @classmethod
    def convert_to_lower(cls, data: Any):
        """Make string settings case insensitive."""
        if isinstance(data, dict):
            for name, value in data.items():
                if (
                    isinstance(value, str)
                    and name in cls.model_fields
                    and cls.model_fields[name].annotation is str
                ):
                    data[name] = value.lower()

        return data

    @model_validator(mode="after")
    def tm_gcp_d4_bug_check(self):
        """Check if the model has basis, func, and prog fields and get their values. Backcheck with GCP basis sets."""
        if all(s in self.model_fields for s in ["func", "basis", "prog"]):
            func: str = getattr(self, "func")
            disp = FUNCTIONALS[func]["disp"]
            gcp_keywords = {
                "minis": "MINIS",
                "sv": "SV",
                "6-31g(d)": "631GD",
                "def2-sv(p)": "SV(P)",
                "def2-svp": "SVP",
                "def2-tzvp": "TZ",
            }
            if (
                getattr(self, "prog") == QmProg.TM
                and getattr(self, "basis") in gcp_keywords
                and disp == "d4"
            ):
                warnings.warn(
                    "Small basis set detected: Current version of TM includes a bug when combining GCP with DFT-D4. Switching to D3."
                )
                setattr(self, "func", func.replace("d4", "d3"))
                self.model_validate(self)

        return self
