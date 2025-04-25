from typing import override, Any

from pydantic import model_validator


from ...params import PLENGTH, GenericConfig
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
        for name, value in data:
            if isinstance(value, str) and cls.model_fields[name].annotation is str:
                data[name] = value.lower()

        return data
