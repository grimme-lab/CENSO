import pytest
from censo.config.parts.optimization import OptimizationConfig
from censo.params import QmProg, GfnVersion, TmSolvMod, OrcaSolvMod


def test_optimization_config_default_values():
    """Test default values of OptimizationConfig"""
    config = OptimizationConfig()

    assert config.prog == QmProg.TM
    assert config.func == "r2scan-3c"
    assert config.basis == "def2-mtzvpp"
    assert config.sm == TmSolvMod.DCOSMORS
    assert config.gfnv == GfnVersion.GFN2
    assert config.optcycles == 8
    assert config.maxcyc == 200
    assert config.optlevel == "normal"
    assert config.threshold == 3.0
    assert config.gradthr == 0.01
    assert config.hlow == 0.01
    assert config.macrocycles is True
    assert config.constrain is False
    assert config.xtb_opt is True
    assert config.template is False


def test_positive_integer_validation():
    """Test validation of positive integer parameters"""
    with pytest.raises(ValueError):
        OptimizationConfig(optcycles=0)
    with pytest.raises(ValueError):
        OptimizationConfig(optcycles=-1)
    with pytest.raises(ValueError):
        OptimizationConfig(maxcyc=0)
    with pytest.raises(ValueError):
        OptimizationConfig(maxcyc=-1)


def test_positive_float_validation():
    """Test validation of positive float parameters"""
    with pytest.raises(ValueError):
        OptimizationConfig(threshold=0)
    with pytest.raises(ValueError):
        OptimizationConfig(threshold=-1.0)
    with pytest.raises(ValueError):
        OptimizationConfig(gradthr=0)
    with pytest.raises(ValueError):
        OptimizationConfig(gradthr=-0.01)
    with pytest.raises(ValueError):
        OptimizationConfig(hlow=0)
    with pytest.raises(ValueError):
        OptimizationConfig(hlow=-0.01)


@pytest.mark.parametrize(
    "optlevel",
    ["crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"],
)
def test_valid_optlevels(optlevel):
    """Test all valid optimization levels"""
    config = OptimizationConfig(optlevel=optlevel)
    assert config.optlevel == optlevel


def test_invalid_optlevel():
    """Test invalid optimization level"""
    with pytest.raises(ValueError):
        OptimizationConfig(optlevel="invalid")  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "sm", [TmSolvMod.DCOSMORS, TmSolvMod.COSMO, OrcaSolvMod.CPCM, OrcaSolvMod.SMD]
)
def test_valid_solvent_models(sm):
    """Test valid solvent model combinations"""
    config = OptimizationConfig(sm=sm)
    assert config.sm == sm


def test_invalid_solvent_model():
    """Test invalid solvent model"""
    with pytest.raises(ValueError):
        OptimizationConfig(sm=TmSolvMod.COSMORS)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "prog,func,should_pass",
    [
        (QmProg.TM, "r2scan-3c", True),
        (QmProg.ORCA, "r2scan-3c", True),
        (QmProg.TM, "invalid-func", False),
        (QmProg.ORCA, "invalid-func", False),
    ],
)
def test_functional_validation(prog, func, should_pass):
    """Test functional validation with different program combinations"""
    if should_pass:
        config = OptimizationConfig(prog=prog, func=func)
        assert config.prog == prog
        assert config.func == func
    else:
        with pytest.raises(ValueError):
            OptimizationConfig(prog=prog, func=func)


def test_constraints_validation():
    """Test constraints validation"""
    # Valid: constraints with xtb_opt enabled
    config = OptimizationConfig(constrain=True, xtb_opt=True)
    assert config.constrain is True
    assert config.xtb_opt is True

    # Invalid: constraints without xtb_opt
    with pytest.raises(
        ValueError, match="Constraints can currently only be used with ANCOPT"
    ):
        OptimizationConfig(constrain=True, xtb_opt=False)
