import pytest
from censo.config.parts.refinement import RefinementConfig
from censo.params import QmProg, GfnVersion, TmSolvMod, OrcaSolvMod


def test_refinement_config_default_values():
    """Test default values of RefinementConfig"""
    config = RefinementConfig()

    assert config.prog == QmProg.TM
    assert config.func == "wb97m-v"
    assert config.basis == "def2-tzvp"
    assert config.sm == TmSolvMod.COSMORS
    assert config.gfnv == GfnVersion.GFN2
    assert config.threshold == 0.95
    assert config.gsolv_included is False
    assert config.template is False


def test_threshold_validation():
    """Test threshold validation"""
    # Test valid thresholds
    config = RefinementConfig(threshold=0.5)
    assert config.threshold == 0.5

    # Test invalid thresholds
    with pytest.raises(ValueError):
        RefinementConfig(threshold=0)
    with pytest.raises(ValueError):
        RefinementConfig(threshold=-0.5)
    with pytest.raises(ValueError):
        RefinementConfig(threshold=1.1)


@pytest.mark.parametrize("sm", [TmSolvMod.COSMORS, OrcaSolvMod.CPCM, OrcaSolvMod.SMD])
def test_valid_solvent_models(sm):
    """Test valid solvent model combinations"""
    config = RefinementConfig(sm=sm)
    assert config.sm == sm


@pytest.mark.parametrize(
    "prog,func,should_pass",
    [
        (QmProg.TM, "wb97m-v", True),
        (QmProg.ORCA, "wb97x-v", True),
        (QmProg.TM, "invalid-func", False),
        (QmProg.ORCA, "invalid-func", False),
    ],
)
def test_functional_validation(prog, func, should_pass):
    """Test functional validation with different program combinations"""
    if should_pass:
        config = RefinementConfig(prog=prog, func=func)
        assert config.prog == prog
        assert config.func == func
    else:
        with pytest.raises(ValueError):
            RefinementConfig(prog=prog, func=func)


def test_gsolv_included_validation(caplog):
    """Test gsolv_included validation with COSMORS"""

    # Test that gsolv_included is set to False when COSMORS is used and gsolv_included is True
    with caplog.at_level("WARNING"):
        config = RefinementConfig(sm=TmSolvMod.COSMORS, gsolv_included=True)
        assert config.gsolv_included is False
        assert (
            "Found cosmors as solvation model but gsolv_included is set to True"
            in caplog.text
        )

    # Test COSMORS_FINE
    caplog.clear()
    with caplog.at_level("WARNING"):
        config2 = RefinementConfig(sm=TmSolvMod.COSMORS_FINE, gsolv_included=True)
        assert config2.gsolv_included is False
        assert (
            "Found cosmors-fine as solvation model but gsolv_included is set to True"
            in caplog.text
        )

    # Test that no warning when gsolv_included is already False
    caplog.clear()
    with caplog.at_level("WARNING"):
        config3 = RefinementConfig(sm=TmSolvMod.COSMORS, gsolv_included=False)
        assert config3.gsolv_included is False
        assert "gsolv_included is set to True" not in caplog.text

    # Test that no change for other solvation models
    config4 = RefinementConfig(sm=OrcaSolvMod.CPCM, gsolv_included=True)
    assert config4.gsolv_included is True
