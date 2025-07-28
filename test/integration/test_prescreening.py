import pytest
from censo.ensembleopt.prescreening import prescreening
from censo.ensemble import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig
from censo.params import QmProg


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
    parallel_config: ParallelConfig,
    prog: QmProg,
    gas_phase: bool,
):
    config.general.gas_phase = gas_phase
    config.prescreening.prog = prog
    config = PartsConfig.model_validate(config, context={"check": "prescreening"})
    timing = prescreening(ensemble_from_xyz, config, parallel_config, cut=True)

    assert timing is not None
    # Add more assertions to verify the output
