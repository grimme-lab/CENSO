import pytest
from pathlib import Path
import json

from censo.config.parts.general import GeneralConfig
from censo.ensemble import EnsembleData
from censo.molecules import Contributions


@pytest.fixture
def two_confs(fixtures_path: Path) -> EnsembleData:
    """Ensemble of the first two conformers from the crest fixture."""
    ensemble = EnsembleData()
    ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix(), nconf=2)
    return ensemble


def _set(ensemble: EnsembleData, **per_conf: Contributions) -> None:
    ensemble.update_contributions(per_conf)


class TestTemperatureContributions:
    """Tests for per-temperature contributions and populations"""

    def test_gtot_at_uses_map_and_falls_back(self, two_confs):
        conf = two_confs.conformers[0]
        conf.update(
            Contributions(
                energy=-1.0, gsolv=-0.01, grrho=0.1, grrho_t={250.0: 0.2, 298.15: 0.1}
            )
        )
        assert conf.gtot_at(250.0) == pytest.approx(-1.0 - 0.01 + 0.2)
        assert conf.gtot_at(298.15) == pytest.approx(conf.gtot)
        # 300 K not in the map: scalar fallback
        assert conf.gtot_at(300.0) == pytest.approx(conf.gtot)
        assert conf.gsolv_at(250.0) == -0.01

    def test_populations_flip_with_temperature(self, two_confs):
        c1, c2 = two_confs.conformers
        _set(
            two_confs,
            **{
                c1.name: Contributions(
                    energy=-1.0, grrho=0.0, grrho_t={250.0: 0.0, 400.0: 0.0}
                ),
                c2.name: Contributions(
                    energy=-1.002, grrho=0.0, grrho_t={250.0: 0.0, 400.0: 0.005}
                ),
            },
        )
        assert (
            max(two_confs.get_populations(250.0).items(), key=lambda x: x[1])[0]
            == c2.name
        )
        assert (
            max(two_confs.get_populations(400.0).items(), key=lambda x: x[1])[0]
            == c1.name
        )

    def test_dump_json_read_output_roundtrip(self, two_confs, tmp_path):
        c1, c2 = two_confs.conformers
        maps = {
            c.name: Contributions(
                energy=-1.0 - i,
                gsolv=-0.01,
                grrho=0.1,
                gsolv_t={273.15: -0.02, 298.15: -0.01},
                grrho_t={273.15: 0.11, 298.15: 0.1},
            )
            for i, c in enumerate(two_confs.conformers)
        }
        _set(two_confs, **maps)
        file = tmp_path / "dump.json"
        two_confs.dump_json(file)

        reloaded = EnsembleData()
        reloaded.conformers = [c1, c2]
        for conf in reloaded.conformers:
            conf.update(Contributions())
        reloaded.read_output(file)

        for conf in reloaded.conformers:
            assert conf.gsolv_t == maps[conf.name].gsolv_t
            assert conf.grrho_t == maps[conf.name].grrho_t
            assert conf.energy == maps[conf.name].energy

    def test_property_json_goes_to_property_results(self, two_confs, tmp_path):
        for conf in two_confs.conformers:
            conf.update(Contributions(energy=-5.0))
        data = {
            conf.name: {"shieldings": [[0, 30.0]], "couplings": []}
            for conf in two_confs.conformers
        }
        file = tmp_path / "4_NMR.json"
        file.write_text(json.dumps({"part_name": "nmr", "data": data, "settings": {}}))

        two_confs.read_output(file)

        assert two_confs.property_results["nmr"] == data
        assert all(conf.energy == -5.0 for conf in two_confs.conformers)

    def test_temperature_range(self, two_confs):
        general = GeneralConfig(multitemp=True, trange=(250, 300, 25))
        grid = [250.0, 275.0, 298.15, 300.0]

        # multitemp off
        assert two_confs.temperature_range(GeneralConfig(), "screening") == []

        # no temperature-dependent data at all
        assert two_confs.temperature_range(general, "screening") == []

        # rrho only at a single temperature
        for conf in two_confs.conformers:
            conf.update(Contributions(energy=-1.0, grrho=0.1))
        with pytest.raises(RuntimeError, match="single temperature"):
            two_confs.temperature_range(general, "nmr")

        # temperature missing for one conformer
        c1, c2 = two_confs.conformers
        full = {t: 0.1 for t in grid}
        c1.update(Contributions(energy=-1.0, grrho=0.1, grrho_t=full))
        c2.update(
            Contributions(energy=-1.0, grrho=0.1, grrho_t={250.0: 0.1, 298.15: 0.1})
        )
        with pytest.raises(RuntimeError, match=r"missing for \[275.0, 300.0\]"):
            two_confs.temperature_range(general, "nmr")

        c2.update(Contributions(energy=-1.0, grrho=0.1, grrho_t=full))
        assert two_confs.temperature_range(general, "nmr") == grid


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
