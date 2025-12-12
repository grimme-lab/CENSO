"""Tests for PartsConfig"""

import pytest
from unittest.mock import patch
from censo.config.parts_config import PartsConfig
from censo.config.parts import (
    GeneralConfig,
    PrescreeningConfig,
    ScreeningConfig,
    OptimizationConfig,
    RefinementConfig,
    NMRConfig,
    UVVisConfig,
)
from censo.params import OrcaSolvMod, TmSolvMod, QmProg


def test_parts_config_default_initialization():
    """Test that PartsConfig initializes with default values"""
    config = PartsConfig()

    # Test that all parts are initialized with their default configurations
    assert isinstance(config.general, GeneralConfig)
    assert isinstance(config.prescreening, PrescreeningConfig)
    assert isinstance(config.screening, ScreeningConfig)
    assert isinstance(config.optimization, OptimizationConfig)
    assert isinstance(config.refinement, RefinementConfig)
    assert isinstance(config.nmr, NMRConfig)
    assert isinstance(config.uvvis, UVVisConfig)


def test_parts_config_str_representation():
    """Test the string representation of PartsConfig"""
    config = PartsConfig()
    str_repr = str(config)

    # Verify that the string representation contains all parts
    assert str_repr.count("\n") > 0  # Should have multiple lines
    for part_name in [
        "general",
        "prescreening",
        "screening",
        "optimization",
        "refinement",
        "nmr",
        "uvvis",
    ]:
        assert part_name in str_repr.lower()


@pytest.mark.parametrize(
    "solvent,sm_model,should_pass",
    [
        ("water", OrcaSolvMod.CPCM, True),  # Valid combination
        ("chloroform", OrcaSolvMod.SMD, True),  # Valid combination
        ("invalid_solvent", OrcaSolvMod.CPCM, False),  # Invalid solvent
    ],
)
def test_solvent_model_validation(solvent, sm_model, should_pass):
    """Test solvent model validation with different combinations"""
    config = PartsConfig()
    config.general.solvent = solvent

    # Set the same solvent model for all parts that use it
    config.screening.sm = sm_model
    config.optimization.sm = sm_model
    config.refinement.sm = sm_model
    config.nmr.sm = sm_model
    config.uvvis.sm = sm_model

    if should_pass:
        # Should not raise any validation errors
        config.model_validate(config, context={"check_all": True, "check_paths": False})
    else:
        # Should raise ValueError for invalid combinations
        with pytest.raises(ValueError, match="not available with"):
            PartsConfig.model_validate(config, context={"check_all": True})


def test_custom_config_values():
    """Test PartsConfig with custom values"""
    custom_config = PartsConfig(
        general=GeneralConfig(solvent="water"),
        screening=ScreeningConfig(sm=OrcaSolvMod.CPCM),
        optimization=OptimizationConfig(sm=OrcaSolvMod.CPCM),
    )

    assert custom_config.general.solvent == "water"
    assert custom_config.screening.sm == OrcaSolvMod.CPCM
    assert custom_config.optimization.sm == OrcaSolvMod.CPCM


def test_invalid_solvent_combination():
    """Test that invalid solvent/model combinations raise appropriate errors"""
    config = PartsConfig()
    config.general.solvent = "dmf"
    config.screening.sm = OrcaSolvMod.CPCM
    config.optimization.sm = TmSolvMod.DCOSMORS  # Different solvent model

    # Should raise ValueError due to incompatible solvent/model combination
    with pytest.raises(ValueError):
        config.model_validate(config, context={"check_all": True})


def test_paths_model_validation():
    """Test that missing paths raise proper validation errors when check_paths=True"""
    config = PartsConfig()
    # Do not set required paths; assume default config is missing some required paths
    with pytest.raises(ValueError, match="is not set in the configuration"):
        PartsConfig.model_validate(config, context={"check_all": True})


def test_parts_validation():
    config = PartsConfig()
    config.screening.func = "invalid"
    with pytest.raises(ValueError):
        PartsConfig.model_validate(config, context={"check": ["screening"]})

    # Assert that PartsConfig._sm_check and PartsConfig._paths_check are called once each
    config = PartsConfig()
    config.screening.gsolv_included = True

    with (
        patch.object(config, "_sm_check", wraps=config._sm_check) as mock_sm_check,
        patch.object(
            config, "_paths_check", wraps=config._paths_check
        ) as mock_paths_check,
    ):
        config = PartsConfig.model_validate(
            config, context={"check": ["screening"], "check_paths": False}
        )
        mock_sm_check.assert_called_once()
        # _paths_check is not called when check_paths=False
        mock_paths_check.assert_not_called()

    assert not config.screening.gsolv_included

    # Now test that both are called when paths validation is enabled
    # Use model_construct to bypass path validation
    from censo.config.paths import PathsConfig

    config = PartsConfig()
    config.paths = PathsConfig.model_construct(
        tm="/fake/path/tm",
        cosmotherm="/fake/path/cosmotherm",
        cosmorssetup="fake_setup",
    )

    with (
        patch.object(config, "_sm_check", wraps=config._sm_check) as mock_sm_check,
        patch.object(
            config, "_paths_check", wraps=config._paths_check
        ) as mock_paths_check,
        patch("censo.config.paths.Path") as mock_path,
    ):
        # Mock Path to make validation pass
        mock_path.return_value.is_file.return_value = True
        config = PartsConfig.model_validate(config, context={"check": ["screening"]})
        mock_sm_check.assert_called_once()
        mock_paths_check.assert_called_once()


def test_screening_requires_cosmors_and_turbomole_paths():
    """Test that screening with COSMO-RS requires both cosmotherm, cosmorssetup, and turbomole paths"""
    config = PartsConfig()
    config.screening.sm = TmSolvMod.COSMORS
    config.screening.prog = QmProg.TM

    # Without the required paths, should raise ValueError
    with pytest.raises(
        ValueError, match="Path for '(cosmotherm|cosmorssetup|tm)'.*is not set"
    ):
        PartsConfig.model_validate(
            config,
            context={"check": ["screening"], "check_sm": False},
        )


def test_context_specific_path_validation_for_all_parts():
    """Test that path validation is context-specific for all parts.

    This test ensures that when validating with context={'check': '<part>'}, only the
    program paths required for that specific part are validated, not paths for other parts.

    We set up a config where different parts use different programs and validate each
    part separately to ensure only its required paths are checked.
    """
    from censo.config.paths import PathsConfig

    # Configure all parts with different programs
    config = PartsConfig()

    # Set parts to use ORCA
    config.prescreening.prog = QmProg.ORCA
    config.nmr.prog = QmProg.ORCA
    config.uvvis.prog = QmProg.ORCA

    # Set parts to use TM
    config.screening.prog = QmProg.TM
    config.screening.sm = TmSolvMod.COSMORS  # Requires cosmotherm/cosmorssetup
    config.optimization.prog = QmProg.TM
    config.refinement.prog = QmProg.TM
    config.rot.prog = QmProg.TM

    # Disable features that would require additional paths
    config.optimization.xtb_opt = False  # Avoid requiring xtb for optimization
    config.general.gas_phase = True  # Avoid requiring xtb for RRHO and prescreening

    # Test 1: Validate screening with only TM and COSMO-RS paths (no ORCA)
    config.paths = PathsConfig.model_construct(
        tm="/fake/path/tm",
        cosmotherm="/fake/path/cosmotherm",
        cosmorssetup="fake_setup",
    )

    with patch("censo.config.paths.Path") as mock_path:
        mock_path.return_value.is_file.return_value = True
        mock_path.return_value.resolve.return_value = mock_path.return_value
        mock_path.return_value.parent = mock_path.return_value
        mock_path.return_value.__truediv__.return_value = mock_path.return_value

        # Screening should pass (has TM and COSMO-RS paths)
        validated = PartsConfig.model_validate(
            config, context={"check": ["screening"], "check_sm": False}
        )
        assert validated.screening.prog == QmProg.TM

        # Optimization should pass (has TM path, no COSMO-RS needed)
        validated = PartsConfig.model_validate(
            config, context={"check": ["optimization"], "check_sm": False}
        )
        assert validated.optimization.prog == QmProg.TM

        # Refinement should pass (has TM path)
        validated = PartsConfig.model_validate(
            config, context={"check": ["refinement"], "check_sm": False}
        )
        assert validated.refinement.prog == QmProg.TM

        # Rot should pass (has TM path)
        validated = PartsConfig.model_validate(
            config, context={"check": ["rot"], "check_sm": False}
        )
        assert validated.rot.prog == QmProg.TM

        # NMR should fail (needs ORCA, which is not provided)
        with pytest.raises(ValueError, match="orca.*is not set"):
            PartsConfig.model_validate(
                config, context={"check": ["nmr"], "check_sm": False}
            )

        # UVVis should fail (needs ORCA)
        with pytest.raises(ValueError, match="orca.*is not set"):
            PartsConfig.model_validate(
                config, context={"check": ["uvvis"], "check_sm": False}
            )

    # Test 2: Validate prescreening, NMR, and UVVis with only ORCA path (no TM or COSMO-RS)
    config.paths = PathsConfig.model_construct(
        orca="/fake/path/orca",
    )

    with (
        patch("censo.config.paths.Path") as mock_path,
        patch("builtins.open", create=True) as mock_open,
    ):
        mock_path.return_value.is_file.return_value = True
        mock_file = mock_open.return_value.__enter__.return_value
        mock_file.read.return_value = b"Program Version 5.0.3"

        # Prescreening should pass (has ORCA, gas_phase=True so no xtb needed)
        validated = PartsConfig.model_validate(
            config, context={"check": ["prescreening"], "check_sm": False}
        )
        assert validated.prescreening.prog == QmProg.ORCA

        # NMR should pass (has ORCA)
        validated = PartsConfig.model_validate(
            config, context={"check": ["nmr"], "check_sm": False}
        )
        assert validated.nmr.prog == QmProg.ORCA

        # UVVis should pass (has ORCA)
        validated = PartsConfig.model_validate(
            config, context={"check": ["uvvis"], "check_sm": False}
        )
        assert validated.uvvis.prog == QmProg.ORCA

        # Screening should fail (needs TM and COSMO-RS paths)
        with pytest.raises(
            ValueError, match="(tm|cosmotherm|cosmorssetup).*is not set"
        ):
            PartsConfig.model_validate(
                config, context={"check": ["screening"], "check_sm": False}
            )

        # Optimization should fail (needs TM)
        with pytest.raises(ValueError, match="tm.*is not set"):
            PartsConfig.model_validate(
                config, context={"check": ["optimization"], "check_sm": False}
            )

        # Refinement should fail (needs TM)
        with pytest.raises(ValueError, match="tm.*is not set"):
            PartsConfig.model_validate(
                config, context={"check": ["refinement"], "check_sm": False}
            )

        # Rot should fail (needs TM)
        with pytest.raises(ValueError, match="tm.*is not set"):
            PartsConfig.model_validate(
                config, context={"check": ["rot"], "check_sm": False}
            )
