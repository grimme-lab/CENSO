"""Tests for BasePartConfig"""

import enum
import warnings

from censo.config.parts.base import BasePartConfig
from censo.params import PLENGTH, QmProg


def test_base_part_config_str():
    """Test the string representation of BasePartConfig"""

    class TestEnum(enum.Enum):
        """Test enum"""

        A = "a"
        B = "b"

    class TestBaseConfig(BasePartConfig):
        """Test implementation of BasePartConfig"""

        test_value: str = "test"
        test_number: int = 42
        test_enum: TestEnum = TestEnum.A

    config = TestBaseConfig()
    result = str(config)

    # Check header formatting
    assert "TestBase" in result

    # Check value formatting
    assert "test_value" in result
    assert "test" in result
    assert "test_number" in result
    assert "42" in result
    assert "TestEnum" not in result

    # Check alignment
    lines = result.split("\n")
    for line in lines[1:]:  # Skip header
        if ":" in line:
            name = line.split(":")[0]
            assert len(name) < PLENGTH // 2


def test_tm_gcp_d4_bug_check():
    """Test the tm_gcp_d4_bug_check validator"""

    class TestConfigWithFields(BasePartConfig):
        """Test config with func, basis, prog fields"""

        func: str = "pbe-d4"
        basis: str = "def2-sv(p)"
        prog: QmProg = QmProg.TM

    # Test case that triggers the warning and correction
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        config = TestConfigWithFields(basis="def2-sv(p)", func="pbe-d4", prog=QmProg.TM)
        assert config.func == "pbe-d3"  # Should be corrected
        assert len(w) == 1
        assert "Small basis set detected" in str(w[0].message)

    # Test case that does not trigger (different basis)
    config2 = TestConfigWithFields(basis="def2-qzvp", func="pbe-d4", prog=QmProg.TM)
    assert config2.func == "pbe-d4"  # Should not be changed

    # Test case that does not trigger (different disp)
    config3 = TestConfigWithFields(basis="def2-sv(p)", func="pbe-d3", prog=QmProg.TM)
    assert config3.func == "pbe-d3"  # Should not be changed

    # Test case that does not trigger (different prog)
    config4 = TestConfigWithFields(basis="def2-sv(p)", func="pbe-d4", prog=QmProg.ORCA)
    assert config4.func == "pbe-d4"  # Should not be changed
