import pytest
from pathlib import Path

from censo.params import QmProg, XtbSolvMod, TmSolvMod, OrcaSolvMod


@pytest.fixture
def fixtures_path():
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def example_ensemble_file():
    return Path(__file__).parent / "fixtures" / "crest_conformers.xyz"


@pytest.fixture(autouse=True)
def tmp_wd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    yield


@pytest.fixture
def mock_functionals():
    """Mock functional definitions for testing"""
    return {
        "pbe-d4": {
            QmProg.TM: "actual_value",
            QmProg.ORCA: "actual_value",
            "disp": "d4",
            "type": "gga",
        },
        "r2scan-3c": {
            QmProg.TM: "actual_value",
            QmProg.ORCA: "actual_value",
            "disp": "3c",
            "type": "mgga",
        },
    }


@pytest.fixture
def mock_solvents():
    """Mock solvent definitions for testing"""
    return {
        "h2o": {
            XtbSolvMod.GBSA: "actual_value",
            TmSolvMod.COSMORS: "actual_value",
            OrcaSolvMod.CPCM: "actual_value",
        },
        "chcl3": {
            XtbSolvMod.GBSA: "actual_value",
            TmSolvMod.COSMORS: "actual_value",
            OrcaSolvMod.CPCM: "actual_value",
        },
    }
