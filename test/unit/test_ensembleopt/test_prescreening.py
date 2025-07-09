import pytest
from unittest.mock import patch

from censo.config.parts_config import PartsConfig
from censo.ensembledata import EnsembleData
from censo.ensembleopt.prescreening import prescreening, jsonify
from censo.params import QmProg

# ============= Tests for Core Prescreening Functionality =============


class TestPrescreening:
    """Tests for core prescreening functionality"""

    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_gas_phase(
        self,
        mock_execute,
        mock_factory,
        mock_ensemble: EnsembleData,
        mock_execute_results,
    ):
        """Test prescreening function in gas phase"""
        # Set up gas phase
        config = PartsConfig()
        config.general.gas_phase = True

        # Mock execute results for sp only (gas phase)
        mock_execute.return_value = mock_execute_results["prescreening"]["sp"]

        # Run prescreening
        prescreening(mock_ensemble, config, None)

        # Verify calls
        assert mock_execute.call_count == 1  # Only sp calculation
        assert mock_factory[QmProg].create.call_count == 1

    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_solution(
        self,
        mock_execute,
        mock_factory,
        mock_ensemble: EnsembleData,
        mock_execute_results,
    ):
        """Test prescreening function with solvation"""
        # Mock execute results for both xtb_gsolv and sp
        mock_execute.side_effect = [
            mock_execute_results["prescreening"]["xtb_gsolv"],
            mock_execute_results["prescreening"]["sp"],
        ]

        config = PartsConfig()

        # Run prescreening
        prescreening(mock_ensemble, config, None)

        # Verify calls
        assert mock_execute.call_count == 2  # Both xtb_gsolv and sp calculations
        assert mock_factory[QmProg].create.call_count == 2

    @pytest.mark.parametrize(
        "threshold,expected_count",
        [
            (1.0, 5),  # Small threshold should remove high energy conformers
            (1000.0, 74),  # Large threshold should keep all conformers
        ],
    )
    @patch("censo.ensembleopt.prescreening.Factory")
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_threshold(
        self,
        mock_execute,
        mock_factory,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        threshold: float,
        expected_count: int,
    ):
        """Test energy threshold-based conformer removal"""
        config = PartsConfig()
        config.prescreening.threshold = threshold

        # Mock execute results
        mock_execute.side_effect = [
            mock_execute_results["prescreening"]["xtb_gsolv"],
            mock_execute_results["prescreening"]["sp"],
        ]

        # Run prescreening
        prescreening(mock_ensemble, config, None)

        # Verify number of remaining conformers
        assert len(mock_ensemble.conformers) == expected_count

    # @patch("censo.ensembleopt.prescreening.Factory")
    # @patch("censo.ensembleopt.prescreening.execute")
    # def test_prescreening_empty_ensemble(self, mock_execute, mock_factory, mock_config):
    #     """Test prescreening with empty ensemble"""
    #     ensemble = MockEnsembleData([])
    #     mock_execute.return_value = ({}, None)  # Return empty results
    #
    #     # Run prescreening
    #     with pytest.raises(ValueError, match="empty"):
    #         prescreening(ensemble, mock_config, ncores=1, omp=1)

    # @patch("censo.ensembleopt.prescreening.Factory")
    # @patch("censo.ensembleopt.prescreening.execute")
    # def test_prescreening_single_conformer(
    #     self, mock_execute, mock_factory, mock_config
    # ):
    #     """Test prescreening with single conformer"""
    #     conf = MockMoleculeData(
    #         name="CONF1", energy=-100.0, gsolv=-0.01, gtot=-100.01, bmw=1.0
    #     )
    #     ensemble = MockEnsembleData([conf])
    #
    #     # Mock execute results
    #     mock_execute.side_effect = [
    #         ({"CONF1": MockResult(gsolv=-0.01)}, None),
    #         ({"CONF1": MockResult(energy=-100.0)}, None),
    #     ]
    #
    #     # Run prescreening
    #     prescreening(ensemble, mock_config, ncores=1, omp=1)
    #
    #     # Verify single conformer was processed
    #     assert len(ensemble.conformers) == 1
    #     assert mock_execute.call_count == 2
    #


# ============= Tests for Result Handling =============


class TestResultHandling:
    """Tests for result handling and output functions"""

    def test_jsonify(self, mock_ensemble: EnsembleData):
        """Test JSON conversion function"""

        config = PartsConfig()

        # Run jsonify
        result = jsonify(mock_ensemble, config.prescreening)

        # Verify result
        assert isinstance(result, dict)
        assert "data" in result
        assert "settings" in result
        assert "CONF1" in result["data"]
