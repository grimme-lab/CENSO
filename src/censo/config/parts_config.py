from typing import override, Any
from pydantic import model_validator, Field, field_validator, ValidationInfo
import warnings


from ..params import OrcaSolvMod, QmProg, TmSolvMod
from .generic import GenericConfig
from ..assets import SOLVENTS
from .parts import (
    GeneralConfig,
    PrescreeningConfig,
    ScreeningConfig,
    OptimizationConfig,
    RefinementConfig,
    NMRConfig,
    UVVisConfig,
    RotConfig,
)
from .paths import PathsConfig


class PartsConfig(GenericConfig):
    """
    Class to store all part-related settings for CENSO in one place.
    """

    general: GeneralConfig = Field(default_factory=GeneralConfig)
    """General settings"""

    prescreening: PrescreeningConfig = Field(default_factory=PrescreeningConfig)
    """Prescreening settings"""

    screening: ScreeningConfig = Field(default_factory=ScreeningConfig)
    """Screening settings"""

    optimization: OptimizationConfig = Field(default_factory=OptimizationConfig)
    """Optimization settings"""

    refinement: RefinementConfig = Field(default_factory=RefinementConfig)
    """Refinement settings"""

    nmr: NMRConfig = Field(default_factory=NMRConfig)
    """NMR settings"""

    rot: RotConfig = Field(default_factory=RotConfig)
    """Optical rotation settings"""

    uvvis: UVVisConfig = Field(default_factory=UVVisConfig)
    """UV/Vis settings"""

    paths: PathsConfig = Field(default_factory=PathsConfig.model_construct)
    """Paths settings"""

    def _selected_parts(self, context: dict[str, Any] | None):
        """
        Helper to yield (name, part) tuples for the parts to check or skip depending on context.

        :param context: The context dictionary.
        :type context: dict[str, Any] | None
        :returns: Generator of (name, part) tuples.
        :rtype: Generator[tuple[str, Any], None, None]
        """
        available = [
            ("prescreening", self.prescreening),
            ("screening", self.screening),
            ("optimization", self.optimization),
            ("refinement", self.refinement),
            ("nmr", self.nmr),
            ("rot", self.rot),
            ("uvvis", self.uvvis),
        ]
        if context:
            check = context.get("check", [])
            check_all = context.get("check_all", False)
            if check:
                return [(name, part) for name, part in available if name in check]
            elif check_all:
                return available
        return list()

    @override
    def __str__(self):
        """Create a formatted string for printing the settings"""
        return str("\n".join(f"{config}" for (_, config) in self))

    @field_validator("paths", mode="before")
    def setup_paths_without_validation(cls, value: PathsConfig | dict[str, set[str]]):
        """
        Set up paths without validation if a dict is passed.

        :param value: The paths configuration or dict.
        :return: The paths configuration instance.
        """
        if isinstance(value, dict):
            return PathsConfig.model_construct(**value)
        return value

    @model_validator(mode="after")
    def validate_parts_sm_and_paths(self, info: ValidationInfo):
        """
        Validate part configs, including solvent models and paths for selected parts.

        :param info: Validation info.
        :return: The validated instance.
        """
        context = info.context
        if context:
            parts_to_check = self._selected_parts(context)
            self._parts_check(parts_to_check)

            check_paths = context.get("check_paths", True)
            check_sm = context.get("check_sm", True)
            if check_sm:
                self._sm_check(parts_to_check)
            if check_paths:
                self._paths_check(parts_to_check)
        else:
            warnings.warn(
                "No context found in config. Skipping solvent and program path validation."
            )

        return self

    def _parts_check(self, parts_to_check: list[tuple[str, Any]]):
        """
        Call attribute validators.
        """
        for name, part in parts_to_check:
            checked = part.model_validate(part)
            setattr(self, name, checked)

    # SOLVENT/SM VALIDATION
    # NOTE: since solvent is a general settings this is validated here because we need access
    # to this setting

    def _sm_check(self, parts_to_check: list[tuple[str, Any]]):
        solvent: str = self.general.solvent
        for name, part in parts_to_check:
            solvent_model: OrcaSolvMod | TmSolvMod | None = getattr(part, "sm", None)
            if solvent_model:
                available_solvents = [
                    s for s, keywords in SOLVENTS.items() if solvent_model in keywords
                ]
                if solvent not in available_solvents:
                    raise ValueError(
                        f"Solvent {solvent} not available with {solvent_model} in {name}."
                    )

        return self

    # PATHS VALIDATION
    # NOTE: we need to know which programs are going to be used before we check the paths

    def _paths_check(self, parts_to_check: list[tuple[str, Any]]):
        required_progs: set[str] = set()
        for name, part in parts_to_check:
            # Check for main program
            prog: QmProg | None = getattr(part, "prog", None)
            if prog is not None:
                required_progs.add(prog.value)

            # Special cases
            # NOTE: in principle you would also need to check if a part that actually runs xtb_rrho is checked
            if (
                (name == "general" and part.evaluate_rrho)
                or (name == "prescreening" and not self.general.gas_phase)
                or (name == "optimization" and part.xtb_opt)
            ):
                required_progs.add("xtb")

            # Check for solvent model specific programs
            sm: TmSolvMod | OrcaSolvMod | None = getattr(part, "sm", None)
            if (
                sm is not None
                and sm
                in [
                    TmSolvMod.COSMORS,
                    TmSolvMod.COSMORS_FINE,
                ]
                and not self.general.gas_phase
            ):
                required_progs.add("cosmotherm")
                required_progs.add("cosmorssetup")

        # Now check if the required paths are actually set and run the validators
        for p in required_progs:
            path = getattr(self.paths, p, None)
            if not path:
                raise ValueError(
                    f"Path for '{p}' is required but it is not set in the configuration."
                )
            # Re-assign to trigger validation since `validate_assignment` is True
            # on failure this should raise a ValueError
            setattr(self.paths, p, path)

        return self
