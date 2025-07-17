import pytest
import os
from pathlib import Path
import shutil

from censo.logging import set_filehandler, set_loglevel


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
