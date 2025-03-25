from typing import Literal
from pydantic import field_validator, Field, model_validator

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

    @model_validator(mode="after")
    def func_must_be_known_in_prog(self):
        prog: str = self.prog
        try:
            assert FUNCTIONALS[self.func][prog] is not None
            assert FUNCTIONALS[self.func]["disp"] is not None
            assert FUNCTIONALS[self.func]["type"] is not None
        except (KeyError, AssertionError):
            raise ValueError(
                f"Functional {self.func} not (fully) defined for prog {prog}."
            )
        return self

    @field_validator("active_nuclei")
    @classmethod
    def active_nuclei_must_be_available(cls, v: str):
        nuclei = [s.lower() for s in v.split(",")]
        if not all(n in ["h", "c", "f", "si", "p"] for n in nuclei):
            raise ValueError(
                f"Selection of active nuclei ({v}) invalid. Allowed nuclei: h, c, f, si, p."
            )
        return v
