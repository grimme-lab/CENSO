from pathlib import Path
from unittest.mock import patch
from censo.config.setup import configure


def test_missing_paths_filled_in(fixtures_path: Path, example_ensemble_file: Path):
    """Test that find_program_paths is called and its results are in PartsConfig.paths"""
    rcpath = fixtures_path / "example.censo2rc"

    # Mock find_program_paths to return specific paths
    mock_paths = {
        "xtb": "/path/to/xtb",
        "orca": "",
        "tm": "",
        "cosmotherm": "",
        "cosmorssetup": "",
    }

    with patch("censo.config.setup.find_program_paths", return_value=mock_paths):
        config = configure(rcpath=str(rcpath))

        # Verify that xtb path is set correctly in PartsConfig.paths
        assert config.paths.xtb == "/path/to/xtb"
