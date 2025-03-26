from typing import override
from pydantic import model_validator


from ..params import OrcaSolvMod, TmSolvMod, GenericConfig

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


class PartsConfig(GenericConfig):
    """
    Class to store all part-related settings for CENSO in one place.
    """

    general: GeneralConfig = GeneralConfig()
    """General settings"""

    prescreening: PrescreeningConfig = PrescreeningConfig()
    """Prescreening settings"""

    screening: ScreeningConfig = ScreeningConfig()
    """Screening settings"""

    optimization: OptimizationConfig = OptimizationConfig()
    """Optimization settings"""

    refinement: RefinementConfig = RefinementConfig()
    """Refinement settings"""

    nmr: NMRConfig = NMRConfig()
    """NMR settings"""

    uvvis: UVVisConfig = UVVisConfig()
    """UV/Vis settings"""

    @override
    def __str__(self):
        """Create a formatted string for printing the settings"""
        return str("\n".join(f"{config}" for (_, config) in self))

    # TODO: validate only active parts (run=True) if called from cli, otherwise check all parts
    # TODO: how to handle the case of setting settings, temporarily yielding invalid settings

    # SOLVENT/SM VALIDATION
    # NOTE: since solvent is a general settings this is validated here because we need access
    # to this setting

    @model_validator(mode="after")
    def sm_check(self):
        solvent: str = self.general.solvent

        for name, part in [
            ("screening", self.screening),
            ("optimization", self.optimization),
            ("refinement", self.refinement),
            ("nmr", self.nmr),
            ("uvvis", self.uvvis),
        ]:
            solvent_model: OrcaSolvMod | TmSolvMod = part.sm
            available_solvents = [
                s for s, keywords in SOLVENTS.items() if solvent_model in keywords
            ]
            if solvent not in available_solvents:
                raise ValueError(
                    f"Solvent {solvent} not available with {solvent_model} in {name}."
                )

        return self
