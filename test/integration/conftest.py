import pytest
import os
from pathlib import Path
import shutil
from censo.ensembledata import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig
from censo.params import QmProg
from censo.config.setup import find_program_paths


def pytest_runtest_setup(item):
    if "requires_xtb" in item.keywords and not has_xtb():
        pytest.skip("xtb is not present in your path.")
    if "requires_orca" in item.keywords and not has_orca():
        pytest.skip("ORCA is not present in your path.")
    if "requires_turbomole" in item.keywords and not has_turbomole():
        pytest.skip("Turbomole (ridft) is not present in your path.")


# Utility functions for program availability checks


def has_xtb():
    program_paths = find_program_paths()
    return program_paths.get("xtb", "") != ""


def has_orca():
    program_paths = find_program_paths()
    return program_paths.get("orca", "") != ""


def has_turbomole():
    return shutil.which("ridft") is not None


@pytest.fixture(autouse=True)
def set_program_paths():
    """
    Always set GenericProc.paths, but do not skip tests globally
    """
    program_paths = find_program_paths()
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
def parallel_config():
    ncores = os.cpu_count() or 4
    return ParallelConfig(ncores=ncores, omp=1)


@pytest.fixture
def parts_config_orca():
    config = PartsConfig()
    config.prescreening.prog = QmProg.ORCA
    config.general.gas_phase = True
    return config


@pytest.fixture
def parts_config_tm():
    config = PartsConfig()
    config.prescreening.prog = QmProg.TM
    config.general.gas_phase = False
    return config
