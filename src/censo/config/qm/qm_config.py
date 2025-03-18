from typing import Optional, Literal
from pydantic import BaseModel


from ...params import GfnVersion, XtbSolvMod, GridLevel, OrcaSolvMod


class QmConfig(BaseModel):
    """
    Base class for QM processor configurations.
    Contains general/xtb-related settings.
    """

    gfnv: GfnVersion
    copy_mo: bool
    gas_phase: bool
    sm_rrho: Optional[XtbSolvMod]
    solvent: Optional[str]
    consider_sym: Optional[bool]
    multitemp: Optional[bool]
    temperature: Optional[float]
    trange: Optional[tuple[float, float, float]]
    sthr: Optional[float]
    imagthr: Optional[float]
    bhess: Optional[bool]
    rmsdbias: Optional[bool]

    # DFT related
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
    optcycles: Optional[int]
    hlow: Optional[float]

    # NMR
    couplings: Optional[bool]
    shieldings: Optional[bool]

    # UVVis
    nroots: Optional[int]
