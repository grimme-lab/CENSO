import pytest
from censo.config.parts.prescreening import PrescreeningConfig
from censo.params import GfnVersion, QmProg


def test_prescreening_config_default_values():
    """Test default values of PrescreeningConfig"""
    config = PrescreeningConfig()

    assert config.prog == QmProg.TM
    assert config.func == "pbe-d3"  # should be changed to d3 from d4 because of GCP bug
    assert config.basis == "def2-sv(p)"
    assert config.gfnv == GfnVersion.GFN2
    assert config.threshold == 4.0
    assert config.template is False


def test_threshold_validation():
    """Test threshold validation"""
    with pytest.raises(ValueError):
        PrescreeningConfig(threshold=0.0)
    with pytest.raises(ValueError):
        PrescreeningConfig(threshold=-1.0)


@pytest.mark.parametrize(
    "prog,func,should_pass",
    [
        (QmProg.TM, "pbe-d3", True),
        (QmProg.ORCA, "pbe-d4", True),  # Assuming this is valid
        (QmProg.TM, "invalid-func", False),
    ],
)
def test_functional_validation(prog, func, should_pass):
    """Test functional validation with different combinations"""
    if should_pass:
        config = PrescreeningConfig(prog=prog, func=func)
        assert config.prog == prog
        assert config.func == func
    else:
        with pytest.raises(ValueError):
            PrescreeningConfig(prog=prog, func=func)
