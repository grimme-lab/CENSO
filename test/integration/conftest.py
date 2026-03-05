import pytest
import os
from pathlib import Path
import shutil
from censo.config.paths import PathsConfig
from censo.ensemble import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.setup import find_program_paths
from censo.parallel import get_cluster


@pytest.fixture
def ensemble_from_xyz(tmp_path: Path) -> EnsembleData:
    def load_xyz(filepath: str) -> EnsembleData:
        ensemble = EnsembleData()
        ensemble.read_input(filepath, nconf=10)  # read only 10 conformers
        return ensemble

    src = Path(__file__).parent / "fixtures" / "crest_conformers_small.xyz"
    dst = tmp_path / src.name
    shutil.copy(src, dst)
    return load_xyz(str(dst))


@pytest.fixture(scope="session")
def client():
    cluster = get_cluster()
    client = cluster.get_client()
    yield client


@pytest.fixture
def config():
    config = PartsConfig()
    config.paths = PathsConfig.model_construct(None, **find_program_paths())
    return config


@pytest.fixture(autouse=True)
def sync_worker_cwd(tmp_path, client):
    client.run(os.chdir, str(tmp_path))
