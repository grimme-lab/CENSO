"""Tests for GeneralConfig"""

import pytest
from censo.config.parts.general import GeneralConfig
from censo.params import XtbSolvMod


def test_general_config_default_values():
    """Test default values of GeneralConfig"""
    config = GeneralConfig()

    assert config.temperature == 298.15
    # assert config.multitemp is True
    # assert config.trange == (273.15, 373.15, 5)
    assert config.evaluate_rrho is True
    assert config.sm_rrho == XtbSolvMod.GBSA
    # assert config.consider_sym is True
    # assert config.bhess is True
    # assert config.rmsdbias is False
    assert config.imagthr == -100.0
    assert config.sthr == 50.0
    assert config.solvent == "h2o"
    assert config.gas_phase is False
    assert config.copy_mo is True
    assert config.balance is True


def test_temperature_validation():
    """Test temperature validation"""
    with pytest.raises(ValueError):
        GeneralConfig(temperature=-1.0)


# def test_trange_validation():
#     """Test trange validation"""
#     # Test negative temperature
#     with pytest.raises(ValueError, match="Negative temperature"):
#         GeneralConfig(trange=(-1.0, 373.15, 5))
#
#     # Test wrong order
#     with pytest.raises(ValueError, match="must be larger than"):
#         GeneralConfig(trange=(373.15, 273.15, 5))
#
#     # Test invalid step size
#     with pytest.raises(ValueError, match="must be positive"):
#         GeneralConfig(trange=(273.15, 373.15, 0))
#     with pytest.raises(ValueError, match="must be positive"):
#         GeneralConfig(trange=(273.15, 373.15, -1))


def test_imagthr_validation():
    """Test imagthr validation"""
    with pytest.raises(ValueError):
        GeneralConfig(imagthr=0.0)


@pytest.mark.parametrize(
    "solvent,sm_rrho,should_pass",
    [
        ("h2o", XtbSolvMod.GBSA, True),
        ("invalid_solvent", XtbSolvMod.GBSA, False),
    ],
)
def test_solvent_validation(solvent, sm_rrho, should_pass):
    """Test solvent validation with different combinations"""
    if should_pass:
        config = GeneralConfig(solvent=solvent, sm_rrho=sm_rrho)
        assert config.solvent == solvent
        assert config.sm_rrho == sm_rrho
    else:
        with pytest.raises(ValueError):
            GeneralConfig(solvent=solvent, sm_rrho=sm_rrho)
