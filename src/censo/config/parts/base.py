from typing import override


from ...params import PLENGTH, GenericConfig
from ...utilities import h2


class BasePartConfig(GenericConfig):
    """Base class for configuration classes for all parts"""

    @override
    def __str__(self):
        lines: list[str] = []
        lines.append(h2(f"{self.__class__.__name__}"))
        for name, value in self:
            lines.append(f"{name} : {value}".center(PLENGTH, " "))

        return str("\n".join(lines))
