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
