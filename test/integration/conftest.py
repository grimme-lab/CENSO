import pytest
import os
from pathlib import Path
import shutil
from censo.ensembledata import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig
from censo.params import QmProg


@pytest.fixture(autouse=True)
def tmp_wd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    orig = os.getcwd()
    monkeypatch.chdir(tmp_path)
    yield
    os.chdir(orig)


@pytest.fixture
def ensemble_from_xyz(tmp_path: Path) -> EnsembleData:
    def load_xyz(filepath: str) -> EnsembleData:
        ensemble = EnsembleData()
        ensemble.read_input(filepath)
        return ensemble

    src = Path(__file__).parent / "fixtures" / "crest_conformers_small.xyz"
    dst = tmp_path / src.name
    shutil.copy(src, dst)
    return load_xyz(str(dst))


@pytest.fixture
def mock_parallel_config():
    return ParallelConfig(ncores=4, omp=2, ompmin=1, ompmax=4)


@pytest.fixture
def mock_parts_config_orca():
    config = PartsConfig()
    config.prescreening.prog = QmProg.ORCA
    config.general.gas_phase = True
    return config


@pytest.fixture
def mock_parts_config_tm():
    config = PartsConfig()
    config.prescreening.prog = QmProg.TM
    config.general.gas_phase = False
    return config
