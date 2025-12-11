import pytest
from unittest.mock import patch

from censo.config.parts_config import PartsConfig
from censo.ensemble import EnsembleData
from censo.ensembleopt.prescreening import prescreening, jsonify

# ============= Tests for Core Prescreening Functionality =============


class TestPrescreening:
    """Tests for core prescreening functionality"""

    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_gas_phase(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        parallel_setup,
        skip_paths_validation,
    ):
        """Test prescreening function in gas phase and rrho"""
        # Set up gas phase
        config = PartsConfig()
        config.general.gas_phase = True

        # Mock execute results for sp only (gas phase)
        mock_execute.return_value = mock_execute_results["prescreening"]["sp"]

        # Run prescreening
        client, cluster = parallel_setup
        prescreening(mock_ensemble, config, client)

        # Verify calls
        assert mock_execute.call_count == 1

    @pytest.mark.parametrize(
        "threshold,expected_count",
        [
            (1.0, 5),  # Small threshold should remove high energy conformers
            (1000.0, 74),  # Large threshold should keep all conformers
        ],
    )
    @patch("censo.ensembleopt.prescreening.execute")
    def test_prescreening_threshold(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        threshold: float,
        expected_count: int,
        parallel_setup,
        skip_paths_validation,
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
        client, cluster = parallel_setup
        prescreening(mock_ensemble, config, client)

        assert mock_execute.call_count == 2

        # Verify number of remaining conformers
        assert len(mock_ensemble.conformers) == expected_count


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
