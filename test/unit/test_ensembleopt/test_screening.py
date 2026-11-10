import json
from pathlib import Path

import pytest
from unittest.mock import patch

from censo.config.parts_config import PartsConfig
from censo.ensemble import EnsembleData
from censo.ensembleopt.screening import screening, jsonify
from censo.molecules import Contributions
from censo.processing.results import RRHOResult


def _multitemp_rrho(ensemble: EnsembleData, grid: list[float]) -> dict:
    """Fake mRRHO results with a conformer-dependent temperature slope."""
    return {
        conf.name: RRHOResult(
            energy=0.08,
            gibbs={t: 0.08 + 1e-5 * (i % 3) * (t - grid[-1]) for t in grid},
        )
        for i, conf in enumerate(ensemble.conformers)
    }


class TestScreeningMultitemp:
    """Tests for multitemp evaluation and --retemp in screening"""

    @patch("censo.ensembleopt.screening.execute")
    def test_multitemp_gas_phase(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        mock_execute_results,
        parallel_setup,
        skip_paths_validation,
        tmp_path,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        config = PartsConfig()
        config.general.gas_phase = True
        config.general.multitemp = True
        grid = config.general.get_temperatures()

        mock_ensemble.remove_conformers(
            lambda conf: conf.name not in mock_execute_results["screening"]["sp"]
        )
        mock_execute.side_effect = [
            mock_execute_results["screening"]["sp"],
            _multitemp_rrho(mock_ensemble, grid),
        ]

        client, cluster = parallel_setup
        results_dict = screening(mock_ensemble, config, client)

        assert all(len(d["grrho_t"]) == 21 for d in results_dict["data"].values())
        assert len(results_dict["averages_t"]) == 21
        assert all(
            sum(pops.values()) == pytest.approx(1.0)
            for pops in results_dict["populations_t"].values()
        )
        marked = [
            line
            for line in Path("1_SCREENING.out").read_text().splitlines()
            if "<<==part1==" in line
        ]
        assert len(marked) == 1 and marked[0].split()[0] == "298.15"
        saved = json.loads(Path("1_SCREENING.json").read_text())
        assert "averages_t" in saved and "populations_t" in saved

    @patch("censo.ensembleopt.screening.execute")
    def test_retemp_xtb_reuses_energies(
        self,
        mock_execute,
        mock_ensemble: EnsembleData,
        parallel_setup,
        skip_paths_validation,
        tmp_path,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        mock_ensemble.conformers = mock_ensemble.conformers[:5]
        energies = {
            conf.name: -300.0 - 0.001 * i for i, conf in enumerate(mock_ensemble)
        }
        for conf in mock_ensemble:
            conf.update(Contributions(energy=energies[conf.name], grrho=0.07))

        config = PartsConfig()
        config.general.gas_phase = True
        config.general.multitemp = True
        config.general.trange = (250.0, 350.0, 25.0)
        mock_execute.return_value = _multitemp_rrho(
            mock_ensemble, config.general.get_temperatures()
        )

        client, cluster = parallel_setup
        screening(mock_ensemble, config, client, retemp="xtb")

        assert mock_execute.call_count == 1
        assert Path("1_SCREENING_RETEMP.json").is_file()
        assert not Path("1_SCREENING.json").exists()
        assert all(conf.energy == energies[conf.name] for conf in mock_ensemble)
        assert all(len(conf.grrho_t) == 6 for conf in mock_ensemble)


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
