import pytest
from pathlib import Path
import shutil
from censo.ensembledata import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig
from censo.params import QmProg
from censo.config.setup import find_program_paths


@pytest.fixture(autouse=True)
def set_program_paths():
    program_paths = find_program_paths()
    for prog in ["xtb", "orca"]:
        if program_paths[prog] == "":
            raise RuntimeError(f"{prog} is not present in your path.")

    from censo.processing import GenericProc

    GenericProc.paths.update(program_paths)


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
