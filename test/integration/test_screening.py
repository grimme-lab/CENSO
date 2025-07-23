import pytest
from censo.ensembleopt.screening import screening
from censo.ensembledata import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig
from censo.params import SOLV_MODS, OrcaSolvMod, QmProg, TmSolvMod

# Only test orca and tm (not xtb)
PROGS = [QmProg.ORCA.value, QmProg.TM.value]


# Parameterized test cases for screening
@pytest.mark.optional
@pytest.mark.parametrize(
    "prog, solvation_model, gas_phase",
    [
        pytest.param(QmProg.ORCA, None, True, marks=pytest.mark.requires_orca),
        pytest.param(QmProg.TM, None, True, marks=pytest.mark.requires_turbomole),
        *[
            pytest.param(QmProg.ORCA, sm, False, marks=pytest.mark.requires_orca)
            for sm in SOLV_MODS[QmProg.ORCA.value]
        ],
        *[
            pytest.param(
                QmProg.TM, sm, False,
                marks=(pytest.mark.requires_turbomole, pytest.mark.requires_cosmotherm)
                if sm in [TmSolvMod.COSMORS, TmSolvMod.COSMORS_FINE]
                else pytest.mark.requires_turbomole,
            )
            for sm in SOLV_MODS[QmProg.TM.value]
        ],
    ],
)
def test_screening_integration(
    ensemble_from_xyz: EnsembleData,
    parallel_config: ParallelConfig,
    prog: QmProg,
    solvation_model: TmSolvMod | OrcaSolvMod | None,
    gas_phase: bool,
):
    config = PartsConfig()
    config.general.gas_phase = gas_phase
    config.screening.prog = prog
    if not gas_phase and solvation_model is not None:
        config.screening.sm = solvation_model
        config.general.solvent = "h2o"
    timing = screening(ensemble_from_xyz, config, parallel_config, cut=True)
    assert timing is not None
    # Optionally add more output checks
