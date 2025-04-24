from typing import override


from ...params import PLENGTH, GenericConfig
from ...utilities import h2


class BasePartConfig(GenericConfig):
    """Base class for configuration classes for all parts"""

    @override
    def __str__(self):
        lines: list[str] = []
        lines.append(h2(f"{self.__class__.__name__.split("Config")[0]}"))
        for name, value in self:
            lines.append(f"{name:>{PLENGTH / 2}} : {value}")

        return str("\n".join(lines))
