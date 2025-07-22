import pytest
from censo.ensembleopt.prescreening import prescreening
from censo.ensembledata import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig


@pytest.mark.optional
@pytest.mark.parametrize(
    "prog_name, gas_phase",
    [
        pytest.param("orca", True, marks=pytest.mark.requires_orca),
        pytest.param("orca", False, marks=pytest.mark.requires_orca),
        pytest.param("tm", True, marks=pytest.mark.requires_turbomole),
        pytest.param("tm", False, marks=pytest.mark.requires_turbomole),
    ],
)
def test_prescreening_parameterized(
    ensemble_from_xyz: EnsembleData,
    parallel_config: ParallelConfig,
    parts_config_orca: PartsConfig,  # For orca configuration
    parts_config_tm: PartsConfig,  # For tm configuration
    prog_name: str,
    gas_phase: bool,
):
    if prog_name == "orca":
        config = parts_config_orca
    elif prog_name == "tm":
        config = parts_config_tm
    else:
        raise ValueError("Invalid program name")

    config.general.gas_phase = gas_phase
    result = prescreening(ensemble_from_xyz, config, parallel_config, cut=True)

    assert result is not None
    # Add more assertions to verify the output
