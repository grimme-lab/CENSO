from typing import override
from pydantic import BaseModel


from ...params import DIGILEN
from ...utilities import h2


class BasePartConfig(BaseModel):
    """Base class for configuration classes for all parts"""

    @override
    def __str__(self):
        lines: list[str] = []
        lines.append(h2(f"{self.__class__.__name__}"))
        for name, value in self:
            lines.append(f"{name:>{DIGILEN // 2}}: {value}")

        return str("\n".join(lines))
