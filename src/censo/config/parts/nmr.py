from typing import Literal
from pydantic import field_validator, Field, model_validator

from .base import BasePartConfig
from ...params import QmProg, TmSolvMod, OrcaSolvMod
from ...assets import FUNCTIONALS


class NMRConfig(BasePartConfig):
    """Config class for NMR"""

    prog: QmProg = QmProg.ORCA
    """Program that should be used for calculations."""

    func: str = "pbe0-d4"
    """Functional that should be used for calculations."""

    basis: str = "def2-tzvp"
    """Basis set that should be used for calculations."""

    sm: Literal[TmSolvMod.COSMO, TmSolvMod.DCOSMORS] | OrcaSolvMod = OrcaSolvMod.SMD
    """Solvation model that should be used for calculations."""

    resonance_frequency: float = Field(gt=0, default=300.0)
    """Resonance frequency assumed for calculation of coupling constants."""

    ss_cutoff: float = Field(gt=0, default=8.0)
    """Only for ORCA: passed to SpinSpinThresh keyword."""

    fc_only: bool = True
    """Whether to only calculate the Fermi-Contact term for spin-spin couplings."""

    shieldings: bool = True
    """Whether to calculate shieldings."""

    couplings: bool = True
    """Whether to calculate spin-spin couplings."""

    active_nuclei: str = "h,c"
    """Nuclei active for NMR calculations."""

    run: bool = False
    """Whether to run NMR calculations (only for CLI)."""

    template: bool = False
    """Whether to use template files."""

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
