import pytest
from censo.config.parts.nmr import NMRConfig
from censo.params import QmProg, GfnVersion, TmSolvMod, OrcaSolvMod


def test_nmr_config_default_values():
    """Test default values of NMRConfig"""
    config = NMRConfig()

    assert config.prog == QmProg.ORCA
    assert config.func == "pbe0-d4"
    assert config.basis == "def2-tzvp"  # should be lower case
    assert config.sm == OrcaSolvMod.SMD
    assert config.gfnv == GfnVersion.GFN2
    assert config.resonance_frequency == 300.0
    assert config.ss_cutoff == 8.0
    assert config.fc_only is True
    assert config.shieldings is True
    assert config.couplings is True
    assert config.active_nuclei == "h,c"
    assert config.run is False
    assert config.template is False


def test_positive_float_validation():
    """Test validation of positive float parameters"""
    with pytest.raises(ValueError):
        NMRConfig(resonance_frequency=0)
    with pytest.raises(ValueError):
        NMRConfig(resonance_frequency=-300.0)
    with pytest.raises(ValueError):
        NMRConfig(ss_cutoff=0)
    with pytest.raises(ValueError):
        NMRConfig(ss_cutoff=-8.0)


@pytest.mark.parametrize(
    "nuclei,should_pass",
    [
        ("h,c", True),
        ("h,c,f", True),
        ("h,c,f,si,p", True),
        ("H,C", True),  # Should work with uppercase
        ("h", True),
        ("invalid", False),
        ("h,x", False),
        ("h,c,invalid", False),
    ],
)
def test_active_nuclei_validation(nuclei, should_pass):
    """Test validation of active nuclei parameter"""
    if should_pass:
        config = NMRConfig(active_nuclei=nuclei)
        assert config.active_nuclei == nuclei.lower()
    else:
        with pytest.raises(ValueError):
            NMRConfig(active_nuclei=nuclei)


@pytest.mark.parametrize(
    "sm", [TmSolvMod.COSMO, TmSolvMod.DCOSMORS, OrcaSolvMod.CPCM, OrcaSolvMod.SMD]
)
def test_valid_solvent_models(sm):
    """Test valid solvent model combinations"""
    config = NMRConfig(sm=sm)
    assert config.sm == sm


def test_invalid_solvent_model():
    """Test invalid solvent model"""
    with pytest.raises(ValueError):
        NMRConfig(sm=TmSolvMod.COSMORS)  # Not in Literal options
