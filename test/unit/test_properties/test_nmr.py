import pytest
from unittest.mock import patch

from censo.config.parts_config import PartsConfig
from censo.ensemble import EnsembleData
from censo.molecules import Contributions
from censo.properties.nmr import nmr, read_chemeq


@patch("censo.properties.nmr.execute", side_effect=AssertionError("no calculations"))
def test_nmr_retemp_reuses_results(
    mock_execute, fixtures_path, tmp_path, monkeypatch, skip_paths_validation
):
    """--retemp averages the reloaded NMR results over the temperature range."""
    monkeypatch.chdir(tmp_path)
    ensemble = EnsembleData()
    ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix(), nconf=2)

    config = PartsConfig()
    config.general.multitemp = True
    config.general.trange = (273.15, 323.15, 25.0)
    grid = sorted(config.general.get_temperatures())
    for i, conf in enumerate(ensemble.conformers):
        conf.update(
            Contributions(
                energy=-100.0 - 0.001 * i,
                grrho=0.05,
                grrho_t={t: 0.05 + 2e-5 * i * (t - 298.15) for t in grid},
            )
        )
    ensemble.property_results["nmr"] = {
        conf.name: {
            "shieldings": [[0, 30.0 + i], [1, 25.0]],
            "couplings": [[[0, 1], 7.0 + i]],
        }
        for i, conf in enumerate(ensemble.conformers)
    }

    results_dict = nmr(ensemble, config, client=None, retemp="xtb")  # type: ignore[arg-type]

    out = (tmp_path / "4_NMR_RETEMP.out").read_text()
    assert "T = 273.15 K" in out
    assert (tmp_path / "4_NMR_RETEMP.json").is_file()
    assert list(results_dict["averages_t"]) == [str(t) for t in grid]
    sigma_0 = {
        t: dict(av["shieldings"])[0] for t, av in results_dict["averages_t"].items()
    }
    assert sigma_0["273.15"] != pytest.approx(sigma_0["323.15"])


def test_read_chemeq_fixture(fixtures_path, monkeypatch, tmp_path):
    """Test read_chemeq with data based on the fixture file."""
    monkeypatch.chdir(tmp_path)

    # The fixture has two sections (multiple nuclei types)
    # Extract only the first section (first 29 lines) for this test
    fixture_file = fixtures_path / "anmr_nucinfo"
    fixture_lines = fixture_file.read_text().split("\n")
    # Take first 29 lines (header + 14 atom pairs)
    first_section = "\n".join(fixture_lines[:29])

    anmr_file = tmp_path / "anmr_nucinfo"
    anmr_file.write_text(first_section)

    result = read_chemeq()

    # Verify structure - should have 14 atoms based on fixture
    assert len(result) == 14

    # Check specific equivalences from the fixture (converted to 0-indexed)
    assert result[0] == [0]
    assert result[1] == [1, 2]
    assert result[2] == [2, 1]
    assert result[3] == [3]
    assert result[4] == [4, 5]
    assert result[5] == [5, 4]
    assert result[6] == [6, 8]
    assert result[7] == [7]
    assert result[8] == [8, 6]
    assert result[9] == [9]
    assert result[10] == [10, 12]
    assert result[11] == [11, 13]
    assert result[12] == [12, 10]
    assert result[13] == [13, 11]


def test_read_chemeq_file_not_found(tmp_path, monkeypatch):
    """Test read_chemeq raises FileNotFoundError when file doesn't exist."""
    monkeypatch.chdir(tmp_path)

    with pytest.raises(FileNotFoundError):
        read_chemeq()


def test_read_chemeq_single_atom(tmp_path, monkeypatch):
    """Test read_chemeq with a single atom."""
    monkeypatch.chdir(tmp_path)

    content = "          1\n   1   1\n 1"
    anmr_file = tmp_path / "anmr_nucinfo"
    anmr_file.write_text(content)

    result = read_chemeq()

    assert result == {0: [0]}


def test_read_chemeq_all_equivalent(tmp_path, monkeypatch):
    """Test read_chemeq where all atoms are equivalent to each other."""
    monkeypatch.chdir(tmp_path)

    content = "          3\n   1   3\n 1 2 3\n   2   3\n 2 1 3\n   3   3\n 3 1 2"
    anmr_file = tmp_path / "anmr_nucinfo"
    anmr_file.write_text(content)

    result = read_chemeq()

    assert result == {
        0: [0, 1, 2],
        1: [1, 0, 2],
        2: [2, 0, 1],
    }


def test_read_chemeq_returns_dict_with_int_keys(tmp_path, monkeypatch):
    """Test that read_chemeq returns dict with int keys and list[int] values."""
    monkeypatch.chdir(tmp_path)

    content = "          2\n   1   1\n 1\n   2   1\n 2"
    anmr_file = tmp_path / "anmr_nucinfo"
    anmr_file.write_text(content)

    result = read_chemeq()

    # Type checks
    assert isinstance(result, dict)
    for key, value in result.items():
        assert isinstance(key, int)
        assert isinstance(value, list)
        assert all(isinstance(v, int) for v in value)
