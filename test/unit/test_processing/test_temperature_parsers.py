"""Tests for the per-temperature output parsers of xtb and COSMOtherm."""

import json
import math

from censo.params import AU2KCAL, R
from censo.processing.tm_processor import _parse_cosmotherm_tab
from censo.processing.xtb_processor import _parse_thermo_table


def test_parse_thermo_table_last_entry_matches_json(fixtures_path):
    lines = (fixtures_path / "xtb_rrho_multitemp.out").read_text().splitlines()
    enso = json.loads((fixtures_path / "xtb_enso_multitemp.json").read_text())

    gibbs = _parse_thermo_table(lines)

    assert set(gibbs) == {273.15, 298.15, 373.15}
    # xtb_enso.json reports G(T) of the last temp= entry (298.15)
    assert abs(gibbs[298.15] - enso["G(T)"]) < 1e-6
    assert gibbs[373.15] < gibbs[273.15]


def test_parse_thermo_table_without_table():
    assert _parse_thermo_table(["no thermo here\n"]) == {}


def test_parse_cosmotherm_tab():
    lines = [
        "COSMOtherm tab\n",
        " Settings  job   1 : T=   298.150 K ; x(1)= 1.0000 x(2)= 0.0000\n",
        "  Nr Compound      H        ln(gamma)    pv     Gsolv\n",
        "   1 h2o           0.0      0.0          0.0    -6.000\n",
        "   2 out           0.0      0.0          0.0    -5.000\n",
        " Settings  job   2 : T=   350.004 K ; x(1)= 1.0000 x(2)= 0.0000\n",
        "   1 h2o           0.0      0.0          0.0    -5.500\n",
        "   2 out           0.0      0.0          0.0    -4.000\n",
    ]
    videal = 24.789561955 / 298.15

    gsolv = _parse_cosmotherm_tab(lines)

    assert set(gsolv) == {298.15, 350.0}
    for t, value in [(298.15, -5.0), (350.004, -4.0)]:
        vwork = R * t * math.log(videal * t)
        assert math.isclose(gsolv[round(t, 2)], (value + vwork) / AU2KCAL)
