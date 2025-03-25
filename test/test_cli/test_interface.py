import os
import shutil
import pytest
from pathlib import Path

from censo.cli.cml_parser import parse
from censo.cli.interface import startup, entry_point
from censo.params import DESCR


def test_blank_startup():
    entry_point([])


def test_help_startup():
    argv = ["-h"]
    entry_point(argv)


def test_general_startup():
    argv = str("-i testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0").split()
    ensemble, config = startup(parse(argv))


def test_writeconfig():
    argv = str("-newconfig").split()
    entry_point(argv)

    assert Path("censo2rc_NEW").is_file()


def test_writereadconfig():
    argv = str("-newconfig").split()
    entry_point(argv)

    argv = str(
        "-inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0 -inprc censo2rc_NEW"
    ).split()
    ensemble, config = startup(parse(argv))


def test_rc_override():
    argv = str("-newconfig").split()
    entry_point(argv)

    argv = str(
        "-inprc censo2rc_NEW -inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0 -gp"
    ).split()
    args = parse(argv)
    ensemble, config = startup(args)

    assert config.general.gas_phase is True
