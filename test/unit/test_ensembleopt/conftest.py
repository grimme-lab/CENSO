from dataclasses import fields
import json
from collections import defaultdict
from collections.abc import Mapping
import pytest
from pathlib import Path

from censo.ensemble import EnsembleData
from censo.config.job_config import (
    SPResult,
    OptResult,
    GsolvResult,
    RRHOResult,
    NMRResult,
    UVVisResult,
)


def tree():
    return defaultdict(tree)


def to_plain_dict(d):
    if isinstance(d, Mapping):
        return {k: to_plain_dict(v) for k, v in d.items()}
    return d


# ============= Fixtures and Mock Classes =============


@pytest.fixture
def mock_ensemble(fixtures_path: Path):
    """Create a glycerol ensemble"""
    ensemble = EnsembleData()
    ensemble.read_input((fixtures_path / "crest_conformers.xyz").as_posix())

    return ensemble


@pytest.fixture
def mock_execute_results(fixtures_path: Path):
    """Load real results for the glycerol ensemble from fixtures"""
    # Load data
    prescreening = json.loads((fixtures_path / "0_PRESCREENING.json").read_text())[
        "results"
    ]
    screening = json.loads((fixtures_path / "1_SCREENING.json").read_text())["results"]
    optimization = json.loads((fixtures_path / "2_OPTIMIZATION.json").read_text())[
        "results"
    ]
    refinement = json.loads((fixtures_path / "3_REFINEMENT.json").read_text())[
        "results"
    ]
    nmr = json.loads((fixtures_path / "4_NMR.json").read_text())["results"]
    uvvis = json.loads((fixtures_path / "6_UVVIS.json").read_text())["results"]

    # Collect
    results = tree()
    for conf in prescreening:
        res = SPResult(**prescreening[conf]["sp"])
        results["prescreening"]["sp"][conf] = res
        res = GsolvResult(**prescreening[conf]["xtb_gsolv"])
        results["prescreening"]["xtb_gsolv"][conf] = res

    for conf in screening:
        res = GsolvResult(**screening[conf]["gsolv"])
        results["screening"]["gsolv"][conf] = res
        res = SPResult("", res.energy_gas)
        results["screening"]["sp"][conf] = res
        res = RRHOResult(**screening[conf]["xtb_rrho"])
        results["screening"]["xtb_rrho"][conf] = res

    for conf in optimization:
        res = OptResult(**optimization[conf]["opt"])
        results["optimization"]["opt"][conf] = res
        res = RRHOResult(**optimization[conf]["xtb_rrho"])
        results["optimization"]["xtb_rrho"][conf] = res

    for conf in refinement:
        res = SPResult(**refinement[conf]["sp"])
        results["refinement"]["sp"][conf] = res
        res = RRHOResult(**refinement[conf]["xtb_rrho"])
        results["refinement"]["xtb_rrho"][conf] = res

    for conf in nmr:
        res = NMRResult(
            **{k: v for k, v in nmr[conf]["nmr"].items() if k in fields(NMRResult)}
        )
        results["nmr"][conf] = res

    for conf in uvvis:
        res = UVVisResult(
            **{
                k: v
                for k, v in uvvis[conf]["uvvis"].items()
                if k in fields(UVVisResult)
            }
        )
        results["uvvis"][conf] = res

    return to_plain_dict(results)
