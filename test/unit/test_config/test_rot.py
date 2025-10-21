import pytest
from censo.config.parts.rot import RotConfig
from censo.params import QmProg


def test_rot_config_default_values():
    """Test default values of RotConfig"""
    config = RotConfig()

    assert config.prog == QmProg.TM
    assert config.func == "pbe-d4"
    assert config.basis == "def2-svpd"
    assert config.freq == [589.0, 633.0]
    assert config.template is False


@pytest.mark.parametrize(
    "freq,should_pass",
    [
        ([589.0], True),
        ([589.0, 633.0, 1000.0], True),
        ([], False),
        ([0.0], False),
        ([-1.0], False),
        ([589.0, -1.0], False),
    ],
)
def test_freq_validation(freq, should_pass):
    """Test frequency validation"""
    if should_pass:
        config = RotConfig(freq=freq)
        assert config.freq == freq
    else:
        with pytest.raises(ValueError):
            RotConfig(freq=freq)


@pytest.mark.parametrize(
    "freq_str,expected",
    [
        ("[589.0]", [589.0]),
        ("[589.0, 633.0]", [589.0, 633.0]),
    ],
)
def test_freq_casting(freq_str, expected):
    """Test frequency string casting"""
    config = RotConfig(freq=freq_str)  # type: ignore[arg-type]
    assert config.freq == expected


def test_invalid_freq_casting():
    """Test invalid frequency string casting"""
    with pytest.raises(ValueError):
        RotConfig(freq="invalid")  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "func,should_pass",
    [
        ("pbe-d4", True),
        ("r2scan-3c", True),
        ("invalid-func", False),
    ],
)
def test_functional_validation(func, should_pass):
    """Test functional validation for TM program"""
    if should_pass:
        config = RotConfig(func=func)
        assert config.func == func
    else:
        with pytest.raises(ValueError):
            RotConfig(func=func)
