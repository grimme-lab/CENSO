import pytest
from censo.ensembleopt.screening import screening
from censo.ensembledata import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.config.parallel_config import ParallelConfig
from censo.params import SOLV_MODS, OrcaSolvMod, QmProg, TmSolvMod

# Only test orca and tm (not xtb)
PROGS = [QmProg.ORCA.value, QmProg.TM.value]


# Generate test cases (prog_name, solvation_model, gas_phase)
def generate_screening_test_cases() -> list[tuple[str, str | None, bool]]:
    cases: list[tuple[str, str | None, bool]] = []
    # For gas_phase True: no solvation model
    for prog in PROGS:
        cases.append((prog, None, True))
    # For gas_phase False: all solvation models for orca and tm
    for prog in PROGS:
        for solv_mod in SOLV_MODS[prog]:
            cases.append((prog, solv_mod.value, False))
    return cases


@pytest.mark.optional
@pytest.mark.parametrize(
    "prog_name, solvation_model, gas_phase",
    generate_screening_test_cases(),
)
def test_screening_parameterized(
    ensemble_from_xyz: EnsembleData,
    parallel_config: ParallelConfig,
    parts_config_orca: PartsConfig,
    parts_config_tm: PartsConfig,
    prog_name: str,
    solvation_model: TmSolvMod | OrcaSolvMod,
    gas_phase: bool,
):
    if prog_name == QmProg.ORCA.value:
        config = parts_config_orca
    elif prog_name == QmProg.TM.value:
        config = parts_config_tm
    else:
        pytest.skip("Screening only tested for ORCA and TM.")

    config.general.gas_phase = gas_phase
    if not gas_phase and solvation_model is not None:
        config.screening.sm = solvation_model
        config.general.solvent = "h2o"

    result = screening(ensemble_from_xyz, config, parallel_config, cut=True)
    assert result is not None
    # Optionally add more output checks
