import pytest
from censo.ensembleopt.screening import screening
from censo.ensemble import EnsembleData
from censo.config.parts_config import PartsConfig
from censo.params import SOLV_MODS, OrcaSolvMod, QmProg, TmSolvMod
from dask.distributed import Client


# Parameterized test cases for screening
@pytest.mark.optional
@pytest.mark.parametrize(
    "prog, solvation_model, evaluate_rrho",
    [
        pytest.param(QmProg.ORCA, None, True, marks=pytest.mark.requires_orca),
        pytest.param(QmProg.TM, None, True, marks=pytest.mark.requires_turbomole),
        pytest.param(
            QmProg.ORCA, OrcaSolvMod.SMD, True, marks=pytest.mark.requires_orca
        ),
        pytest.param(
            QmProg.TM, TmSolvMod.DCOSMORS, True, marks=pytest.mark.requires_turbomole
        ),
        *[  # type: ignore[var-annotated]
            pytest.param(QmProg.ORCA, sm, False, marks=pytest.mark.requires_orca)
            for sm in SOLV_MODS[QmProg.ORCA.value]
        ],
        *[  # type: ignore[var-annotated]
            pytest.param(
                QmProg.TM,
                sm,
                False,
                marks=(
                    (pytest.mark.requires_turbomole, pytest.mark.requires_cosmotherm)
                    if sm in [TmSolvMod.COSMORS, TmSolvMod.COSMORS_FINE]
                    else pytest.mark.requires_turbomole
                ),
            )
            for sm in SOLV_MODS[QmProg.TM.value]
        ],
    ],
)
def test_screening_integration(
    config: PartsConfig,
    ensemble_from_xyz: EnsembleData,
    client: Client,
    prog: QmProg,
    solvation_model: TmSolvMod | OrcaSolvMod | None,
    evaluate_rrho: bool,
):
    config.general.evaluate_rrho = evaluate_rrho
    config.screening.prog = prog
    if solvation_model is not None:
        config.general.gas_phase = False
        config.screening.sm = solvation_model
        config.general.solvent = "h2o"
        # For COSMORS also test the case where the user might have forgotten to set gsolv_included to False
        if solvation_model == TmSolvMod.COSMORS:
            config.screening.gsolv_included = True
    else:
        config.general.gas_phase = True
    result = screening(ensemble_from_xyz, config, client, cut=True)
    assert result is not None
    # Optionally add more output checks
