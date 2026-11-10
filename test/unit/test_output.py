"""Tests for the shared part output in censo.output."""

import json
from pathlib import Path

from censo.config.parts_config import PartsConfig
from censo.ensemble import EnsembleData
from censo.molecules import Contributions
from censo.output import write_output, write_ranking


def _ensemble(fixtures_path: Path) -> EnsembleData:
    ensemble = EnsembleData()
    ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix(), nconf=3)
    ensemble.update_contributions(
        {
            conf.name: Contributions(-100.0 - 0.001 * i, -0.01, 0.002 * i)
            for i, conf in enumerate(ensemble)
        }
    )
    return ensemble


def test_write_output_json_holds_part_data_and_reloads(fixtures_path: Path):
    """The json dump holds energies, part fields and the part's own settings and can be reloaded."""
    ensemble = _ensemble(fixtures_path)
    config = PartsConfig()

    dump = write_output(
        "4_NMR", ensemble, config, "REPORT", lambda conf: {"nat": conf.geom.nat}
    )

    assert json.loads(Path("4_NMR.json").read_text()) == dump
    assert dump["part_name"] == "nmr"
    assert dump["settings"] == config.nmr.model_dump(mode="json")
    for conf in ensemble:
        assert dump["data"][conf.name] == {
            "energy": conf.energy,
            "gsolv": conf.gsolv,
            "grrho": conf.grrho,
            "gtot": conf.gtot,
            "nat": conf.geom.nat,
        }
    assert Path("4_NMR.out").read_text().startswith("REPORT\n")

    reloaded = EnsembleData()
    reloaded.read_input((fixtures_path / "crest_conformers.xyz").as_posix(), nconf=3)
    reloaded.read_output("4_NMR.json")
    assert [conf.gtot for conf in reloaded] == [conf.gtot for conf in ensemble]


def test_write_ranking_reports_boltzmann_averages(fixtures_path: Path):
    """The ranking report states Boltzmann averages at the configured temperature for the part."""
    ensemble = _ensemble(fixtures_path)
    config = PartsConfig()
    populations = ensemble.get_populations(config.general.temperature)
    av_e = sum(populations[conf.name] * conf.energy for conf in ensemble)
    av_g = sum(populations[conf.name] * conf.gtot for conf in ensemble)

    write_ranking(
        "2_OPTIMIZATION",
        ensemble,
        config,
        "OPTIMIZATION RESULTS",
        geometries="DFT optimized",
    )

    out = Path("2_OPTIMIZATION.out").read_text()
    assert "ensemble on input geometries (DFT optimized):" in out
    assert f"{av_e:.7f}" in out
    assert f"{av_g:.7f}" in out
    assert "<<==part2==" in out
