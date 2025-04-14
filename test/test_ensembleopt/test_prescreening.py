import pytest
from pathlib import Path
from unittest.mock import patch

from censo.ensembleopt.prescreening import prescreening, _write_results, jsonify
from censo.params import QmProg, GfnVersion, XtbSolvMod


# ============= Fixtures and Mock Classes =============


class MockMoleculeData:
    """Mock class for MoleculeData that uses real attributes instead of Mock"""

    def __init__(
        self, name, energy=-100.0, gsolv=-0.01, grrho=0.1, gtot=-99.91, bmw=0.5
    ):
        self.name = name
        self.energy = energy
        self.gsolv = gsolv
        self.grrho = grrho
        self.gtot = gtot
        self.bmw = bmw


class MockEnsembleData:
    """Mock class for EnsembleData that supports iteration"""

    def __init__(self, conformers):
        self.conformers = conformers

    def __iter__(self):
        return iter(self.conformers)

    def remove_conformers(self, cond):
        self.conformers = [conf for conf in self.conformers if not cond(conf)]

    def set_populations(self, temperature):
        for conf in self.conformers:
            conf.bmw = 1.0 / len(self.conformers)

    def dump_xyz(self, path):
        pass


class MockResult:
    """Mock class for QM calculation results"""

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


@pytest.fixture
def mock_ensemble():
    """Create a mock ensemble with two conformers using real-like objects"""
    conf1 = MockMoleculeData(
        name="CONF1", energy=-100.0, gsolv=-0.01, grrho=0.1, gtot=-99.91, bmw=0.7
    )

    conf2 = MockMoleculeData(
        name="CONF2", energy=-99.5, gsolv=-0.02, grrho=0.12, gtot=-99.4, bmw=0.3
    )

    return MockEnsembleData([conf1, conf2])


@pytest.fixture
def mock_config():
    """Create a mock configuration using real-like objects"""

    class MockPrescreeningConfig:
        def __init__(self):
            self.prog = QmProg.TM
            self.func = "pbe-d4"
            self.basis = "def2-SV(P)"
            self.gfnv = GfnVersion.GFN2
            self.threshold = 4.0
            self.template = False

        def model_dump(self):
            return {
                "prog": self.prog.value,
                "func": self.func,
                "basis": self.basis,
                "gfnv": self.gfnv.value,
                "threshold": self.threshold,
                "template": self.template,
            }

    class MockGeneralConfig:
        def __init__(self):
            self.gas_phase = False
            self.solvent = "h2o"
            self.sm_rrho = XtbSolvMod.GBSA
            self.temperature = 298.15
            self.balance = True
            self.copy_mo = True

    class MockPartsConfig:
        def __init__(self):
            self.prescreening = MockPrescreeningConfig()
            self.general = MockGeneralConfig()

    return MockPartsConfig()


@pytest.fixture
def mock_execute_results():
    """Create mock results for execute function using real-like objects"""
    return {
        "xtb_gsolv": (
            {"CONF1": MockResult(gsolv=-0.01), "CONF2": MockResult(gsolv=-0.02)},
            None,
        ),
        "sp": (
            {"CONF1": MockResult(energy=-100.0), "CONF2": MockResult(energy=-99.5)},
            None,
        ),
    }


# ============= Tests for Core Prescreening Functionality =============


class TestPrescreening:
    """Tests for core prescreening functionality"""

    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_gas_phase(
        self,
        mock_execute,
        mock_factory,
        mock_ensemble,
        mock_config,
        mock_execute_results,
    ):
        """Test prescreening function in gas phase"""
        # Set up gas phase
        mock_config.general.gas_phase = True

        # Mock execute results for sp only (gas phase)
        mock_execute.return_value = mock_execute_results["sp"]

        # Run prescreening
        prescreening(mock_ensemble, mock_config, ncores=1, omp=1)

        # Verify calls
        assert mock_execute.call_count == 1  # Only sp calculation
        mock_factory[QmProg].create.assert_called_once()

    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_solution(
        self,
        mock_execute,
        mock_factory,
        mock_ensemble,
        mock_config,
        mock_execute_results,
    ):
        """Test prescreening function with solvation"""
        # Mock execute results for both xtb_gsolv and sp
        mock_execute.side_effect = [
            mock_execute_results["xtb_gsolv"],
            mock_execute_results["sp"],
        ]

        # Run prescreening
        prescreening(mock_ensemble, mock_config, ncores=1, omp=1)

        # Verify calls
        assert mock_execute.call_count == 2  # Both xtb_gsolv and sp calculations
        mock_factory[QmProg].create.assert_called_once()

    @pytest.mark.parametrize(
        "threshold,expected_count",
        [
            (1.0, 1),  # Small threshold should remove high energy conformer
            (1000.0, 2),  # Large threshold should keep all conformers
        ],
    )
    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_threshold(
        self,
        mock_execute,
        mock_factory,
        mock_ensemble,
        mock_config,
        mock_execute_results,
        threshold,
        expected_count,
    ):
        """Test energy threshold-based conformer removal"""
        mock_config.prescreening.threshold = threshold

        # Mock execute results
        mock_execute.side_effect = [
            mock_execute_results["xtb_gsolv"],
            mock_execute_results["sp"],
        ]

        # Run prescreening
        prescreening(mock_ensemble, mock_config, ncores=1, omp=1)

        # Verify number of remaining conformers
        assert len(mock_ensemble.conformers) == expected_count

    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_empty_ensemble(self, mock_execute, mock_factory, mock_config):
        """Test prescreening with empty ensemble"""
        ensemble = MockEnsembleData([])
        mock_execute.return_value = ({}, None)  # Return empty results

        # Run prescreening
        with pytest.raises(ValueError, match="empty"):
            prescreening(ensemble, mock_config, ncores=1, omp=1)

    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_single_conformer(
        self, mock_execute, mock_factory, mock_config
    ):
        """Test prescreening with single conformer"""
        conf = MockMoleculeData(
            name="CONF1", energy=-100.0, gsolv=-0.01, gtot=-100.01, bmw=1.0
        )
        ensemble = MockEnsembleData([conf])

        # Mock execute results
        mock_execute.side_effect = [
            ({"CONF1": MockResult(gsolv=-0.01)}, None),
            ({"CONF1": MockResult(energy=-100.0)}, None),
        ]

        # Run prescreening
        prescreening(ensemble, mock_config, ncores=1, omp=1)

        # Verify single conformer was processed
        assert len(ensemble.conformers) == 1
        assert mock_execute.call_count == 2


# ============= Tests for Result Handling =============


class TestResultHandling:
    """Tests for result handling and output functions"""

    def test_jsonify(self):
        """Test JSON conversion function"""
        conf = MockMoleculeData(
            name="CONF1", energy=-100.0, gsolv=-0.01, grrho=0.1, gtot=-99.91
        )
        ensemble = MockEnsembleData([conf])

        class MockPrescreeningConfig:
            def model_dump(self):
                return {"prog": "tm", "func": "pbe-d4"}

        config = MockPrescreeningConfig()

        # Run jsonify
        result = jsonify(ensemble, config)

        # Verify result
        assert isinstance(result, dict)
        assert "data" in result
        assert "settings" in result
        assert "CONF1" in result["data"]
