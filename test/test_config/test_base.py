"""Tests for BasePartConfig"""

from censo.config.parts.base import BasePartConfig
from censo.params import PLENGTH


def test_base_part_config_str():
    """Test the string representation of BasePartConfig"""

    class TestBaseConfig(BasePartConfig):
        """Test implementation of BasePartConfig"""

        test_value: str = "test"
        test_number: int = 42

    config = TestBaseConfig()
    result = str(config)

    # Check header formatting
    assert "TestBaseConfig" in result

    # Check value formatting
    assert "test_value" in result
    assert "test" in result
    assert "test_number" in result
    assert "42" in result

    # Check alignment
    lines = result.split("\n")
    for line in lines[1:]:  # Skip header
        if ":" in line:
            name = line.split(":")[0]
            assert len(name) >= PLENGTH // 2
