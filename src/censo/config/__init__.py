"""
Contains configuration related classes and functions.

Important note for configuration in general:
Many of the models created with pydantic have field_validators that check other already
validated values to validate certain settings. If you set a setting in a way that (temporarily)
an invalid configuration would result, this will raise an error. If you want to avoid this,
e.g. for setting the functional to something not available for the program you want to use
you need to disable assignment validation manually: <model>.model_config.validate_assignment = False.
"""

from .generic import GenericConfig
from .parts_config import PartsConfig
from .parts import (
    PrescreeningConfig,
    ScreeningConfig,
    OptimizationConfig,
    RefinementConfig,
    NMRConfig,
    UVVisConfig,
    GeneralConfig,
    BasePartConfig,
)

__all__ = [
    "PartsConfig",
    "PrescreeningConfig",
    "ScreeningConfig",
    "OptimizationConfig",
    "RefinementConfig",
    "NMRConfig",
    "UVVisConfig",
    "GeneralConfig",
    "BasePartConfig",
    "GenericConfig",
]
