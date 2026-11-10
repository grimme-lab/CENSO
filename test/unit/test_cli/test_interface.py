import json
from pathlib import Path
from unittest.mock import patch

import pytest

from censo.cli.cml_parser import parse
from censo.cli.interface import startup, entry_point
from censo.params import Returncode


@pytest.fixture(autouse=True)
def mock_find_program_paths():
    """Mock find_program_paths to return empty strings for all program paths."""
    mock_paths = {
        "orca": "",
        "tm": "",
        "xtb": "",
        "cosmotherm": "",
        "cosmorssetup": "",
    }
    with patch("censo.config.setup.find_program_paths", return_value=mock_paths):
        yield


def test_blank_startup(capsys):
    ret = entry_point([])
    captured = capsys.readouterr()
    assert ret == Returncode.OK
    assert "No tasks enabled" in captured.out


def test_help_startup():
    argv = ["-h"]
    entry_point(argv)


def test_general_startup(example_ensemble_file, skip_paths_validation):
    argv = str(
        f"-i {example_ensemble_file} --prescreening --solvent water -c 0 -u 0"
    ).split()
    ensemble, config = startup(parse(argv), {"check": ["prescreening"]})


def test_writeconfig():
    argv = "--new-config".split()
    entry_point(argv)

    assert Path("censo2rc_NEW").is_file()


def test_writereadconfig(example_ensemble_file, skip_paths_validation):
    argv = "--new-config".split()
    entry_point(argv)

    argv = str(
        f"-i {example_ensemble_file} --prescreening --solvent water -c 0 -u 0 --inprc censo2rc_NEW"
    ).split()
    ensemble, config = startup(parse(argv), {"check": ["prescreening"]})


def test_rc_override(example_ensemble_file, skip_paths_validation):
    argv = "--new-config".split()
    entry_point(argv)

    argv = str(
        f"--inprc censo2rc_NEW -i {example_ensemble_file} --prescreening --solvent water -c 1 -u 1 --gas-phase"
    ).split()
    args = parse(argv)
    ensemble, config = startup(args, {"check": ["prescreening"]})

    assert config.general.gas_phase is True
    assert config.general.solvent == "water"
    assert ensemble.conformers[0].charge == 1
    assert ensemble.conformers[0].unpaired == 1


def test_trange_implies_multitemp(example_ensemble_file, skip_paths_validation):
    argv = f"-i {example_ensemble_file} --prescreening --gas-phase --trange 250 300 10".split()
    _, config = startup(parse(argv), {"check": ["prescreening"]})

    assert config.general.multitemp is True
    assert config.general.trange == (250.0, 300.0, 10.0)


def test_multitemp_flag(example_ensemble_file, skip_paths_validation):
    argv = f"-i {example_ensemble_file} --prescreening --gas-phase --multitemp".split()
    _, config = startup(parse(argv), {"check": ["prescreening"]})

    assert config.general.multitemp is True


def test_retemp_parsing():
    assert parse(["--retemp"]).retemp == "xtb"
    assert parse(["--retemp", "cosmors"]).retemp == "cosmors"
    assert parse([]).retemp is None


@pytest.mark.parametrize(
    "flags",
    ["--screening --retemp", "--screening --refinement --retemp --reload x.json"],
)
def test_retemp_argument_errors(example_ensemble_file, flags):
    args = parse(f"-i {example_ensemble_file} {flags}".split())
    context = {
        "check": [p for p in ("screening", "refinement") if getattr(args, p, False)]
    }
    with pytest.raises(SystemExit) as exc:
        startup(args, context)
    assert exc.value.code == Returncode.ARGUMENT_ERROR


def test_property_reload_requires_retemp(
    example_ensemble_file, skip_paths_validation, tmp_path
):
    from censo.ensemble import EnsembleData

    ensemble = EnsembleData()
    ensemble.read_input(example_ensemble_file)
    data: dict[str, dict[str, list]] = {
        conf.name: {"shieldings": [], "couplings": []} for conf in ensemble.conformers
    }
    nmr_json = tmp_path / "4_NMR.json"
    nmr_json.write_text(json.dumps({"part_name": "nmr", "data": data, "settings": {}}))

    args = parse(
        f"-i {example_ensemble_file} --nmr --gas-phase --reload {nmr_json}".split()
    )
    with pytest.raises(SystemExit) as exc:
        startup(args, {"check": ["nmr"]})
    assert exc.value.code == Returncode.ARGUMENT_ERROR
