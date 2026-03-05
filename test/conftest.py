import pytest
import os
from pathlib import Path
import shutil

from censo.logging import set_filehandler, set_loglevel


def _has_xtb() -> bool:
    return shutil.which("xtb") is not None


def _has_orca() -> bool:
    return shutil.which("orca") is not None


def _has_turbomole() -> bool:
    return shutil.which("ridft") is not None


def _has_cosmotherm() -> bool:
    return shutil.which("cosmotherm") is not None


def pytest_runtest_setup(item):
    if "requires_xtb" in item.keywords and not _has_xtb():
        pytest.skip("xtb is not present in your path.")
    if "requires_orca" in item.keywords and not _has_orca():
        pytest.skip("ORCA is not present in your path.")
    if "requires_turbomole" in item.keywords and not _has_turbomole():
        pytest.skip("Turbomole (ridft) is not present in your path.")
    if "requires_cosmotherm" in item.keywords and not _has_cosmotherm():
        pytest.skip("CosmoTherm is not present in your path.")


@pytest.fixture(autouse=True)
def tmp_wd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, request):
    orig = os.getcwd()
    if request.config.getoption("--keep-log"):
        set_filehandler(Path(orig) / "censo.log")
    set_loglevel("DEBUG")
    monkeypatch.chdir(tmp_path)
    yield
    os.chdir(orig)
    # Clean up temporary directory
    for item in tmp_path.iterdir():
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink()


def pytest_addoption(parser):
    parser.addoption(
        "--run-optional", action="store_true", default=False, help="run optional tests"
    )
    parser.addoption(
        "--keep-log",
        action="store_true",
        default=False,
        help="keep the log file during test execution",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-optional"):
        # --run-optional given in cli: do not skip optional tests
        return
    skip_optional = pytest.mark.skip(reason="need --run-optional option to run")
    for item in items:
        if "optional" in item.keywords:
            item.add_marker(skip_optional)
