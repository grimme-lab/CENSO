from dataclasses import fields
import json
from collections import defaultdict
from collections.abc import Mapping
import pytest
from pathlib import Path

from censo.ensemble import EnsembleData
from censo.processing.results import (
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
        res_sp = SPResult(**prescreening[conf]["sp"])
        results["prescreening"]["sp"][conf] = res_sp
        res_gsolv = GsolvResult(**prescreening[conf]["xtb_gsolv"])
        results["prescreening"]["xtb_gsolv"][conf] = res_gsolv

    for conf in screening:
        res_gsolv = GsolvResult(**screening[conf]["gsolv"])
        results["screening"]["gsolv"][conf] = res_gsolv
        res_sp = SPResult("", res_gsolv.energy_gas)
        results["screening"]["sp"][conf] = res_sp
        res_rrho = RRHOResult(**screening[conf]["xtb_rrho"])
        results["screening"]["xtb_rrho"][conf] = res_rrho

    for conf in optimization:
        res_opt = OptResult(**optimization[conf]["opt"])
        results["optimization"]["opt"][conf] = res_opt
        res_rrho = RRHOResult(**optimization[conf]["xtb_rrho"])
        results["optimization"]["xtb_rrho"][conf] = res_rrho

    for conf in refinement:
        res_sp = SPResult(**refinement[conf]["sp"])
        results["refinement"]["sp"][conf] = res_sp
        res_rrho = RRHOResult(**refinement[conf]["xtb_rrho"])
        results["refinement"]["xtb_rrho"][conf] = res_rrho

    for conf in nmr:
        res_nmr = NMRResult(
            **{k: v for k, v in nmr[conf]["nmr"].items() if k in fields(NMRResult)}
        )
        results["nmr"][conf] = res_nmr

    for conf in uvvis:
        res_uvvis = UVVisResult(
            **{
                k: v
                for k, v in uvvis[conf]["uvvis"].items()
                if k in fields(UVVisResult)
            }
        )
        results["uvvis"][conf] = res_uvvis

    return to_plain_dict(results)
