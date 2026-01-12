import pytest
from pathlib import Path
import json

from censo.ensemble import EnsembleData


class TestEnsembleDataReadOutput:
    """Tests for EnsembleData.read_output method"""

    def test_read_output_success(self, fixtures_path: Path):
        """Test successful reading of output file with matching conformers"""
        # Create ensemble from crest_conformers.xyz
        ensemble = EnsembleData()
        ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix())

        # Load the reload fixture
        reload_file = fixtures_path / "1_SCREENING.json.reload"

        # Read the reload file to see which conformers it has
        reload_data = json.loads(reload_file.read_text())
        conformer_names_in_reload = set(reload_data["data"].keys())

        # Remove conformers not in reload file
        ensemble.remove_conformers(
            lambda conf: conf.name not in conformer_names_in_reload
        )

        # Store original conformer names
        original_conformers = [conf.name for conf in ensemble.conformers]

        # Read the output file
        ensemble.read_output(reload_file)

        # Verify all conformers are still there
        assert len(ensemble.conformers) == len(original_conformers)

        # Verify that contributions were loaded correctly for a few conformers
        # Check CONF1
        conf1 = next(c for c in ensemble.conformers if c.name == "CONF1")
        assert conf1.energy == pytest.approx(-344.73078844674)
        assert conf1.gsolv == pytest.approx(-0.0060278603308606534)
        assert conf1.grrho == pytest.approx(0.082858561455)

        # Check CONF2
        conf2 = next(c for c in ensemble.conformers if c.name == "CONF2")
        assert conf2.energy == pytest.approx(-344.72393564073)
        assert conf2.gsolv == pytest.approx(-0.011852745882101001)
        assert conf2.grrho == pytest.approx(0.082637163592)

    def test_read_output_missing_conformer(self, fixtures_path: Path):
        """Test that RuntimeError is raised when conformer in ensemble is missing from output"""
        # Create ensemble from crest_conformers.xyz
        ensemble = EnsembleData()
        ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix())

        # Load the reload fixture
        reload_file = fixtures_path / "1_SCREENING.json.reload"

        reload_data = json.loads(reload_file.read_text())
        conformer_names_in_reload = set(reload_data["data"].keys())

        # Keep only conformers NOT in reload file (if any exist)
        all_conformer_names = {conf.name for conf in ensemble.conformers}
        missing_conformers = all_conformer_names - conformer_names_in_reload

        if missing_conformers:
            # Keep only the missing ones
            ensemble.remove_conformers(
                lambda conf: conf.name in conformer_names_in_reload
            )

            # This should raise RuntimeError
            with pytest.raises(
                RuntimeError,
                match="Not all conformers from the current ensemble are found in the output data",
            ):
                ensemble.read_output(reload_file)
        else:
            # If all conformers are in the reload file, we need to add a fake one
            # that won't be in the reload file to trigger the error
            pytest.skip(
                "All conformers from fixture are in reload file, skipping this test case"
            )

    def test_read_output_extra_conformers_in_file(self, fixtures_path: Path):
        """Test that extra conformers in output file are ignored (not an error)"""
        # Create ensemble from crest_conformers.xyz
        ensemble = EnsembleData()
        ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix())

        # Load the reload fixture
        reload_file = fixtures_path / "1_SCREENING.json.reload"

        # Keep only a subset of conformers (e.g., just CONF1 and CONF2)
        ensemble.remove_conformers(lambda conf: conf.name not in ["CONF1", "CONF2"])

        # Should succeed - extra conformers in file are ignored
        ensemble.read_output(reload_file)

        # Verify we still have only 2 conformers
        assert len(ensemble.conformers) == 2

        # Verify the data was loaded correctly
        conf1 = next(c for c in ensemble.conformers if c.name == "CONF1")
        assert conf1.energy == pytest.approx(-344.73078844674)

    def test_read_output_flat_format(self, fixtures_path: Path, tmp_path: Path):
        """Test backward compatibility with flat JSON format (without 'data' key)"""
        # Create ensemble from crest_conformers.xyz
        ensemble = EnsembleData()
        ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix())

        # Keep only a few conformers
        ensemble.remove_conformers(lambda conf: conf.name not in ["CONF1", "CONF2"])

        # Create a flat format JSON file (as produced by dump_json)
        flat_data = {
            "CONF1": {
                "energy": -344.73078844674,
                "gsolv": -0.0060278603308606534,
                "grrho": 0.082858561455,
            },
            "CONF2": {
                "energy": -344.72393564073,
                "gsolv": -0.011852745882101001,
                "grrho": 0.082637163592,
            },
        }

        flat_file = tmp_path / "flat_output.json"
        flat_file.write_text(json.dumps(flat_data, indent=4))

        # Read the flat format file
        ensemble.read_output(flat_file)

        # Verify the data was loaded correctly
        conf1 = next(c for c in ensemble.conformers if c.name == "CONF1")
        assert conf1.energy == pytest.approx(-344.73078844674)
        assert conf1.gsolv == pytest.approx(-0.0060278603308606534)
        assert conf1.grrho == pytest.approx(0.082858561455)

        conf2 = next(c for c in ensemble.conformers if c.name == "CONF2")
        assert conf2.energy == pytest.approx(-344.72393564073)
