"""Tests for BasePartConfig"""

import pytest
from censo.config.parts.base import BasePartConfig
from censo.params import DIGILEN


class TestBaseConfig(BasePartConfig):
    """Test implementation of BasePartConfig"""

    test_value: str = "test"
    test_number: int = 42


def test_base_part_config_str():
    """Test the string representation of BasePartConfig"""
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
            assert len(name) >= DIGILEN // 2
