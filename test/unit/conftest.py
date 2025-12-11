import pytest
from pathlib import Path
import shutil

from censo.params import XtbSolvMod, TmSolvMod, OrcaSolvMod, Prog
from censo.config.parts_config import PartsConfig
from censo.parallel import get_cluster


@pytest.fixture
def fixtures_path():
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def example_ensemble_file(tmp_path: Path):
    """Copies the fixture ensemble file to the test dir and yield the path to that copy."""
    src = Path(__file__).parent / "fixtures" / "crest_conformers.xyz"
    dst = tmp_path / src.name
    shutil.copy(src, dst)
    return dst


# @pytest.fixture(autouse=True)
# def tmp_wd(tmp_path, monkeypatch):
#     orig = os.getcwd()
#     monkeypatch.chdir(tmp_path)
#     yield
#     os.chdir(orig)


@pytest.fixture
def mock_functionals():
    """Mock functional definitions for testing"""
    return {
        "pbe-d4": {
            Prog.TM.value: "actual_value",
            Prog.ORCA.value: "actual_value",
            "disp": "d4",
            "type": "gga",
        },
        "r2scan-3c": {
            Prog.TM.value: "actual_value",
            Prog.ORCA.value: "actual_value",
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


# Patch PartsConfig.model_validate to default context['check_paths'] = False unless explicitly set
# This is an opt-in fixture - use it in tests where you want to skip path validation by default
@pytest.fixture
def skip_paths_validation(monkeypatch):
    original_model_validate = PartsConfig.model_validate

    def patched_model_validate(self, *args, **kwargs):
        # Support both positional and keyword context
        context = kwargs.get("context")
        if context is None and len(args) > 0:
            context = args[0]
        if context is None:
            context = {}
        # If check_paths is not set, default to False
        if "check_paths" not in context:
            context["check_paths"] = False
        # Ensure context is passed as a kwarg for consistency
        kwargs["context"] = context
        return original_model_validate(self, **kwargs)

    monkeypatch.setattr(PartsConfig, "model_validate", patched_model_validate)


@pytest.fixture(scope="session")
def parallel_setup():
    """Provide real parallel setup for tests that need it."""

    cluster = get_cluster()
    client = cluster.get_client()
    return client, cluster
