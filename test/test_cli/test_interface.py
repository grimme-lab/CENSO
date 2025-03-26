from pathlib import Path
import pytest

from censo.cli.cml_parser import parse
from censo.cli.interface import startup, entry_point


def test_blank_startup():
    with pytest.raises(FileNotFoundError, match="No such file or directory"):
        entry_point([])


def test_help_startup():
    argv = ["-h"]
    entry_point(argv)


def test_general_startup(example_ensemble_file):
    argv = str(f"-i {example_ensemble_file} --solvent water -c 0 -u 0").split()
    ensemble, config = startup(parse(argv))


def test_writeconfig():
    argv = str("--new-config").split()
    entry_point(argv)

    assert Path("censo2rc_NEW").is_file()


def test_writereadconfig(example_ensemble_file):
    argv = str("--new-config").split()
    entry_point(argv)

    argv = str(
        f"-i {example_ensemble_file} --solvent water -c 0 -u 0 --inprc censo2rc_NEW"
    ).split()
    ensemble, config = startup(parse(argv))


def test_rc_override(example_ensemble_file):
    argv = str("--new-config").split()
    entry_point(argv)

    argv = str(
        f"--inprc censo2rc_NEW -i {example_ensemble_file} --solvent water -c 1 -u 1 --gas-phase"
    ).split()
    args = parse(argv)
    ensemble, config = startup(args)

    assert config.general.gas_phase is True
    assert config.general.solvent == "water"
    assert ensemble.conformers[0].charge == 1
    assert ensemble.conformers[0].unpaired == 1
