from typing import Literal, Optional
from pydantic import BaseModel, field_validator

from .qm_config import QmConfig
from ...params import OrcaSolvMod, GridLevel
from ...utilities import Factory


class OrcaConfig(QmConfig):
    """Configuration settings for ORCA processor."""

    func: str
    basis: str
    grid: GridLevel
    sm: Optional[OrcaSolvMod]

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
    nroots: Optional[int]


Factory.register_builder("orca_config", OrcaConfig)
