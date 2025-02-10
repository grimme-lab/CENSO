from typing import override
from pydantic import BaseModel, ValidationInfo, field_validator

from censo.assets.assets import FUNCTIONALS

from ..params import OrcaSolvMod, QmProg, TmSolvMod

from ..assets import SOLVENTS

from .parts import (
    GeneralConfig,
    PrescreeningConfig,
    ScreeningConfig,
    OptimizationConfig,
    RefinementConfig,
    NMRConfig,
    UVVisConfig,
)


class PartsConfig(BaseModel):
    """
    Class to store all part-related settings for CENSO in one place.
    """

    general: GeneralConfig = GeneralConfig()
    prescreening: PrescreeningConfig = PrescreeningConfig()
    screening: ScreeningConfig = ScreeningConfig()
    optimization: OptimizationConfig = OptimizationConfig()
    refinement: RefinementConfig = RefinementConfig()
    nmr: NMRConfig = NMRConfig()
    uvvis: UVVisConfig = UVVisConfig()

    @override
    def __str__(self):
        """Create a formatted string for printing the settings"""
        return str("\n".join(f"{config}" for (_, config) in self))

    # TODO: validate solvent, func for all active parts if called from cli, otherwise check all parts
    # TODO: how to handle the case of setting settings, temporarily yielding invalid settings

    # SOLVENT/SM VALIDATION
    # NOTE: since solvent is a general settings this is validated here because we need access
    # to this setting

    @field_validator("screening")
    @classmethod
    def screening_sm_func_check(cls, v: ScreeningConfig, info: ValidationInfo):
        solvent: str = info.data["general"].solvent
        solvent_model: OrcaSolvMod | TmSolvMod = v.sm
        available_solvents = [
            s for s, keywords in SOLVENTS.items() if solvent_model in keywords
        ]
        if solvent not in available_solvents:
            raise ValueError(
                f"Solvent {solvent} not available with {solvent_model} in screening."
            )

    @field_validator("optimization")
    @classmethod
    def optimization_sm_func_check(cls, v: OptimizationConfig, info: ValidationInfo):
        solvent: str = info.data["general"].solvent
        solvent_model: OrcaSolvMod | TmSolvMod = v.sm
        available_solvents = [
            s for s, keywords in SOLVENTS.items() if solvent_model in keywords
        ]
        if solvent not in available_solvents:
            raise ValueError(
                f"Solvent {solvent} not available with {solvent_model} in optimization."
            )

    @field_validator("refinement")
    @classmethod
    def refinement_sm_func_check(cls, v: RefinementConfig, info: ValidationInfo):
        solvent: str = info.data["general"].solvent
        solvent_model: OrcaSolvMod | TmSolvMod = v.sm
        available_solvents = [
            s for s, keywords in SOLVENTS.items() if solvent_model in keywords
        ]
        if solvent not in available_solvents:
            raise ValueError(
                f"Solvent {solvent} not available with {solvent_model} in refinement."
            )

    @field_validator("nmr")
    @classmethod
    def nmr_sm_func_check(cls, v: NMRConfig, info: ValidationInfo):
        solvent: str = info.data["general"].solvent
        solvent_model: OrcaSolvMod | TmSolvMod = v.sm
        available_solvents = [
            s for s, keywords in SOLVENTS.items() if solvent_model in keywords
        ]
        if solvent not in available_solvents:
            raise ValueError(
                f"Solvent {solvent} not available with {solvent_model} in nmr."
            )

    @field_validator("uvvis")
    @classmethod
    def uvvis_sm_func_check(cls, v: UVVisConfig, info: ValidationInfo):
        solvent: str = info.data["general"].solvent
        solvent_model: OrcaSolvMod | TmSolvMod = v.sm
        available_solvents = [
            s for s, keywords in SOLVENTS.items() if solvent_model in keywords
        ]
        if solvent not in available_solvents:
            raise ValueError(
                f"Solvent {solvent} not available with {solvent_model} in uvvis."
            )
