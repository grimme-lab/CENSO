from typing import Optional, Literal
from .qm_config import QmConfig
from ...params import TmSolvMod, GridLevel
from ...utilities import Factory


class TmConfig(QmConfig):
    """Configuration settings for TURBOMOLE processor."""

    func: str
    basis: str
    grid: GridLevel
    sm: Optional[TmSolvMod]

    # Optimization
    macrocycles: Optional[bool]
    optlevel: Optional[
        Literal[
            "crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"
        ]
    ]

    # NMR
    couplings: Optional[bool]
    shieldings: Optional[bool]

    # UVVis
    # nroots: Optional[int]


Factory.register_builder("tm_config", TmConfig)
