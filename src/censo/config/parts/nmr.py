from typing import Literal
from pydantic import field_validator, Field, ValidationInfo

from .base import BasePartConfig
from ...params import QmProg, TmSolvMod, OrcaSolvMod, GfnVersion
from ...assets import FUNCTIONALS


class NMRConfig(BasePartConfig):
    """Config class for NMR"""

    prog: QmProg = QmProg.ORCA
    func: str = "pbe0-d4"
    basis: str = "def2-TZVP"
    sm: Literal[TmSolvMod.COSMO, TmSolvMod.DCOSMORS] | OrcaSolvMod = OrcaSolvMod.SMD
    gfnv: GfnVersion = GfnVersion.GFN2
    resonance_frequency: float = Field(gt=0, default=300.0)
    ss_cutoff: float = Field(gt=0, default=8.0)
    fc_only: bool = True
    shieldings: bool = True
    couplings: bool = True
    active_nuclei: str = "h,c"
    run: bool = False
    template: bool = False

    @field_validator("func")
    @classmethod
    def func_must_be_known_in_prog(cls, v: str, info: ValidationInfo):
        prog: str = info.data["prog"]
        try:
            assert FUNCTIONALS[v][prog] is not None
            assert FUNCTIONALS[v]["disp"] is not None
            assert FUNCTIONALS[v]["type"] is not None
        except (KeyError, AssertionError):
            raise ValueError(f"Functional {v} not (fully) defined for prog {prog}.")
        return v

    @field_validator("active_nuclei")
    @classmethod
    def active_nuclei_must_be_available(cls, v: str):
        nuclei = [s.lower() for s in v.split(",")]
        if not all(n in ["h", "c", "f", "si", "p"] for n in nuclei):
            raise ValueError(
                f"Selection of active nuclei ({v}) invalid. Allowed nuclei: h, c, f, si, p."
            )
        return v
