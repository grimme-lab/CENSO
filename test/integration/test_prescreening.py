import pytest
from censo.ensembleopt.prescreening import prescreening
from censo.ensemble import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.params import QmProg
from dask.distributed import Client


@pytest.mark.optional
@pytest.mark.parametrize(
    "prog, gas_phase",
    [
        pytest.param(QmProg.ORCA, True, marks=pytest.mark.requires_orca),
        pytest.param(QmProg.ORCA, False, marks=pytest.mark.requires_orca),
        pytest.param(QmProg.TM, True, marks=pytest.mark.requires_turbomole),
        pytest.param(QmProg.TM, False, marks=pytest.mark.requires_turbomole),
    ],
)
def test_prescreening_integration(
    config: PartsConfig,
    ensemble_from_xyz: EnsembleData,
    client: Client,
    prog: QmProg,
    gas_phase: bool,
):
    config.general.gas_phase = gas_phase
    config.prescreening.prog = prog
    result = prescreening(ensemble_from_xyz, config, client, cut=True)

    assert result is not None
    # Add more assertions to verify the output
