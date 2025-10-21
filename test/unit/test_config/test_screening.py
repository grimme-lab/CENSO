import pytest
from censo.config.parts.screening import ScreeningConfig
from censo.params import OrcaSolvMod, QmProg, TmSolvMod, GfnVersion


def test_screening_config_default_values():
    """Test default values of ScreeningConfig"""
    config = ScreeningConfig()

    assert config.prog == QmProg.TM
    assert config.func == "r2scan-3c"
    assert config.basis == "def2-mtzvpp"
    assert config.sm == TmSolvMod.COSMORS
    assert config.gfnv == GfnVersion.GFN2
    assert config.threshold == 3.5
    assert config.gsolv_included is False
    assert config.template is False


def test_threshold_validation():
    """Test threshold validation"""
    with pytest.raises(ValueError):
        ScreeningConfig(threshold=0.0)
    with pytest.raises(ValueError):
        ScreeningConfig(threshold=-1.0)


@pytest.mark.parametrize(
    "prog,func,should_pass",
    [
        (QmProg.TM, "r2scan-3c", True),
        (QmProg.ORCA, "r2scan-3c", True),
        (QmProg.TM, "invalid-func", False),
    ],
)
def test_functional_validation(prog, func, should_pass):
    """Test functional validation with different combinations"""
    if should_pass:
        config = ScreeningConfig(prog=prog, func=func)
        assert config.prog == prog
        assert config.func == func
    else:
        with pytest.raises(ValueError):
            ScreeningConfig(prog=prog, func=func)


@pytest.mark.parametrize(
    "prog,sm,should_pass",
    [
        (QmProg.TM, TmSolvMod.COSMORS, True),
        (QmProg.ORCA, OrcaSolvMod.CPCM, True),
        (QmProg.TM, OrcaSolvMod.CPCM, True),
    ],
)
def test_solvent_model_validation(prog, sm, should_pass):
    """Test solvent model validation with different combinations"""
    if should_pass:
        config = ScreeningConfig(prog=prog, sm=sm)
        assert config.prog == prog
        assert config.sm == sm
    else:
        with pytest.raises(ValueError):
            ScreeningConfig(prog=prog, sm=sm)


def test_gsolv_included_validation(caplog):
    """Test gsolv_included validation with COSMORS"""

    # Test that gsolv_included is set to False when COSMORS is used and gsolv_included is True
    with caplog.at_level("WARNING"):
        config = ScreeningConfig(sm=TmSolvMod.COSMORS, gsolv_included=True)
        assert config.gsolv_included is False
        assert (
            "Found cosmors as solvation model but gsolv_included is set to True"
            in caplog.text
        )

    # Test COSMORS_FINE
    caplog.clear()
    with caplog.at_level("WARNING"):
        config2 = ScreeningConfig(sm=TmSolvMod.COSMORS_FINE, gsolv_included=True)
        assert config2.gsolv_included is False
        assert (
            "Found cosmors-fine as solvation model but gsolv_included is set to True"
            in caplog.text
        )

    # Test that no warning when gsolv_included is already False
    caplog.clear()
    with caplog.at_level("WARNING"):
        config3 = ScreeningConfig(sm=TmSolvMod.COSMORS, gsolv_included=False)
        assert config3.gsolv_included is False
        assert "gsolv_included is set to True" not in caplog.text

    # Test that no change for other solvation models
    config4 = ScreeningConfig(sm=OrcaSolvMod.CPCM, gsolv_included=True)
    assert config4.gsolv_included is True
