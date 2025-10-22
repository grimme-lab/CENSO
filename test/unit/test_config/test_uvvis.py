import pytest
from censo.config.parts.uvvis import UVVisConfig
from censo.params import QmProg, OrcaSolvMod


def test_uvvis_config_default_values():
    """Test default values of UVVisConfig"""
    config = UVVisConfig()

    assert config.prog == QmProg.ORCA
    assert config.func == "wb97x-d4"
    assert config.basis == "def2-tzvp"
    assert config.sm == OrcaSolvMod.SMD
    assert config.nroots == 20
    assert config.template is False


def test_nroots_validation():
    """Test validation of nroots parameter"""
    # Test valid nroots
    config = UVVisConfig(nroots=10)
    assert config.nroots == 10

    # Test invalid nroots
    with pytest.raises(ValueError):
        UVVisConfig(nroots=0)
    with pytest.raises(ValueError):
        UVVisConfig(nroots=-1)


def test_program_restriction():
    """Test that only ORCA is allowed as program"""
    config = UVVisConfig(prog=QmProg.ORCA)
    assert config.prog == QmProg.ORCA

    with pytest.raises(ValueError):
        UVVisConfig(prog=QmProg.TM)  # type: ignore[arg-type]


@pytest.mark.parametrize("sm", [OrcaSolvMod.CPCM, OrcaSolvMod.SMD])
def test_valid_solvent_models(sm):
    """Test valid solvent model combinations"""
    config = UVVisConfig(sm=sm)
    assert config.sm == sm


@pytest.mark.parametrize(
    "func,should_pass",
    [
        ("wb97x-v", True),
        ("pbe0-d4", True),
        ("invalid-func", False),
    ],
)
def test_functional_validation(func, should_pass):
    """Test functional validation for ORCA program"""
    if should_pass:
        config = UVVisConfig(func=func)
        assert config.func == func
    else:
        with pytest.raises(ValueError):
            UVVisConfig(func=func)
