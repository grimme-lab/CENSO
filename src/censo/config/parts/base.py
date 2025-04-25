from typing import override


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
            lines.append(f"{name:>{PLENGTH // 2}} : {value}")

        return str("\n".join(lines))
