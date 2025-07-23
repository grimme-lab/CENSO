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
        (QmProg.ORCA, "r2scan-3c", True),  # Assuming this is valid
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
        (QmProg.TM, OrcaSolvMod.CPCM, True),  # Assuming this cross-usage is valid
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
