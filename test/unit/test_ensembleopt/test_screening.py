import pytest
from unittest.mock import patch

from censo.config.parts_config import PartsConfig
from censo.ensemble import EnsembleData
from censo.ensembleopt.screening import screening, jsonify

# ============= Tests for Core Screening Functionality =============


class TestScreening:
    """Tests for core screening functionality"""

    @patch("censo.ensembleopt.screening.execute")
    def test_screening_gas_phase(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        parallel_setup,
        skip_paths_validation,
    ):
        """Test screening function in gas phase"""
        # Set up gas phase
        config = PartsConfig()
        config.general.gas_phase = True

        # Mock execute results for sp only (gas phase)
        mock_execute.side_effect = [
            mock_execute_results["screening"]["sp"],
            mock_execute_results["screening"]["xtb_rrho"],
        ]

        # Prepare ensemble (remove surplus confs)
        mock_ensemble.remove_conformers(
            lambda conf: conf.name not in mock_execute_results["screening"]["sp"]
        )

        # Run screening
        client, cluster = parallel_setup
        screening(mock_ensemble, config, client)

        # Verify calls
        assert mock_execute.call_count == 2  # sp and xtb_rrho

    @patch("censo.ensembleopt.screening.execute")
    def test_screening_gas_phase_norrho(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        parallel_setup,
        skip_paths_validation,
    ):
        """Test screening function in gas phase"""
        # Set up gas phase
        config = PartsConfig()
        config.general.gas_phase = True
        config.general.evaluate_rrho = False

        # Mock execute results for sp only (gas phase)
        mock_execute.side_effect = [
            mock_execute_results["screening"]["sp"],
        ]

        # Prepare ensemble (remove surplus confs)
        mock_ensemble.remove_conformers(
            lambda conf: conf.name not in mock_execute_results["screening"]["sp"]
        )

        # Run screening
        client, cluster = parallel_setup
        screening(mock_ensemble, config, client)

        # Verify calls
        assert mock_execute.call_count == 1

    @patch("censo.ensembleopt.screening.execute")
    def test_screening_solution(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        parallel_setup,
        skip_paths_validation,
    ):
        """Test screening function with solvation"""
        # Mock execute results for gsolv (not included)
        mock_execute.side_effect = [
            mock_execute_results["screening"]["gsolv"],
            mock_execute_results["screening"]["xtb_rrho"],
        ]

        config = PartsConfig()

        # Prepare ensemble (remove surplus confs)
        mock_ensemble.remove_conformers(
            lambda conf: conf.name not in mock_execute_results["screening"]["gsolv"]
        )

        # Run screening
        client, cluster = parallel_setup
        screening(mock_ensemble, config, client)

        # Verify calls
        assert mock_execute.call_count == 2  # Both xtb_gsolv and sp calculations

    @pytest.mark.parametrize(
        "threshold,expected_count",
        [
            (1.0, 12),  # Small threshold should remove high energy conformer
            (1000.0, 45),  # Large threshold should keep all conformers
        ],
    )
    @patch("censo.ensembleopt.screening.execute")
    def test_screening_threshold(
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
        config.screening.threshold = threshold

        # Mock execute results
        mock_execute.side_effect = [
            mock_execute_results["screening"]["gsolv"],
            mock_execute_results["screening"]["xtb_rrho"],
        ]

        # Prepare ensemble (remove surplus confs)
        mock_ensemble.remove_conformers(
            lambda conf: conf.name not in mock_execute_results["screening"]["sp"]
        )

        # Run screening
        client, cluster = parallel_setup
        screening(mock_ensemble, config, client)

        # Verify number of remaining conformers
        assert len(mock_ensemble.conformers) == expected_count


# ============= Tests for Result Handling =============


class TestResultHandling:
    """Tests for result handling and output functions"""

    def test_jsonify(self, mock_ensemble):
        """Test JSON conversion function"""
        config = PartsConfig()

        # Run jsonify
        result = jsonify(mock_ensemble, config.screening)

        # Verify result
        assert isinstance(result, dict)
        assert "data" in result
        assert "settings" in result
        assert "CONF1" in result["data"]
