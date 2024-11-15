import os
import pytest
from censo.cli.cml_parser import parse
from censo.params import DESCR
from censo.ensembledata import EnsembleData


def test_read_input(conf_name_list: tuple[str]):
    """
    Test reading an input
    """
    fixtures_dir = os.path.join(os.path.split(__file__)[0], "fixtures")
    ensemble_dir = os.path.join(fixtures_dir, "crest_conformers.xyz")

    # Test reading
    ensemble = EnsembleData()
    ensemble.read_input(ensemble_dir, charge=2, unpaired=2)
    assert len(ensemble.conformers) == 7
    assert ensemble.runinfo["charge"] == 2
    assert ensemble.runinfo["unpaired"] == 2

    # Test reading only some number of confs
    ensemble = EnsembleData()
    ensemble.read_input(ensemble_dir, nconf=5)
    assert len(ensemble.conformers) == 5
    assert ensemble.runinfo["charge"] == 0
    assert ensemble.runinfo["unpaired"] == 0

    # Test appending confs
    ensemble = EnsembleData()
    ensemble.read_input(ensemble_dir, nconf=5)
    ensemble.read_input(ensemble_dir, append=True, nconf=2)
    assert len(ensemble.conformers) == 7
    assert ensemble.runinfo["charge"] == 0
    assert ensemble.runinfo["unpaired"] == 0

    # Test reading ensemble from censo output
    ensemble_dir = os.path.join(fixtures_dir, "2_OPTIMIZATION.xyz")
    ensemble = EnsembleData()
    ensemble.read_input(ensemble_dir)

    for i, conf in enumerate(ensemble.conformers):
        assert conf.name == conf_name_list[i]


def test_remove_conformers():
    """
    Test removing conformers from the ensemble.
    """
    # Read ensemble from test fixture
    fixtures_dir = os.path.join(os.path.split(__file__)[0], "fixtures")
    ensemble_dir = os.path.join(fixtures_dir, "crest_conformers.xyz")
    ensemble = EnsembleData()
    ensemble.read_input(ensemble_dir)

    # Test removing conformers
    remove_test = ["CONF3", "CONF30"]
    ensemble.remove_conformers(remove_test)
    assert all(
        confname not in [conf.name for conf in ensemble.conformers]
        for confname in remove_test
    )

    # Test removing conformer that is not in ensemble
    ensemble.remove_conformers(remove_test)

    # Test invalid remove
    with pytest.raises(TypeError):
        ensemble.remove_conformers(123)


@pytest.fixture
def conf_name_list() -> tuple[str]:
    """
    Return the list of conformer names in order to be found in 2_OPTIMIZATION.xyz.
    """
    names = (
        "CONF3",
        "CONF1",
        "CONF59",
        "CONF5",
        "CONF10",
        "CONF11",
        "CONF42",
        "CONF31",
        "CONF12",
        "CONF15",
        "CONF21",
        "CONF7",
        "CONF9",
        "CONF6",
        "CONF30",
        "CONF25",
        "CONF18",
        "CONF17",
    )
    return names
