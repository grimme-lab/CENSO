"""
Configuration file for pytest.
"""

import pytest
from pathlib import Path


@pytest.fixture
def example_ensemble_file():
    return Path(__file__).parent / "fixtures" / "crest_conformers.xyz"


@pytest.fixture(autouse=True)
def tmp_wd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    yield
