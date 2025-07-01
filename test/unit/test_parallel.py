import time
import pytest
from unittest.mock import Mock
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.managers import SyncManager
from typing import Any

from censo.parallel import (
    setup_managers,
    ResourceMonitor,
    SPResult,
    GsolvResult,
    RRHOResult,
    OptResult,
    NMRResult,
    UVVisResult,
    MetaData,
    ParallelJob,
    set_omp_tmproc,
    set_omp,
    prepare_jobs,
    execute,
)
from censo.molecules import GeometryData, MoleculeData
from censo.params import GridLevel, QmProg, OMPMIN
from censo.config.job_config import SPJobConfig


# ============= Fixtures and Mock Classes =============


class MockGeometryData:
    """Mock class for GeometryData"""

    def __init__(self, name: str):
        self.name = name
        self.xyz = ["C 0.0 0.0 0.0\n", "H 1.0 0.0 0.0\n"]


@pytest.fixture
def mock_molecule_data():
    """Create a mock MoleculeData instance"""
    mol = Mock(spec=MoleculeData)
    mol.name = "CONF1"
    mol.charge = 0
    mol.unpaired = 0
    mol.mo_paths = {
        "tm": ["path/to/mos"],
        "orca": ["path/to/file.gbw"],
        QmProg.TM: ["path/to/mos"],
        QmProg.ORCA: ["path/to/file.gbw"],
    }
    mol.geom = GeometryData(name="CONF1", xyz=["C 0.0 0.0 0.0\n", "H 1.0 0.0 0.0\n"])
    return mol


@pytest.fixture
def mock_parallel_job(mock_molecule_data):
    """Create a mock ParallelJob instance"""
    return ParallelJob(
        conf=mock_molecule_data.geom,
        charge=mock_molecule_data.charge,
        unpaired=mock_molecule_data.unpaired,
    )


@pytest.fixture
def mock_job_config():
    """Create a mock job configuration"""
    return SPJobConfig(
        copy_mo=False,
        func="pbe-d4",
        basis="def2-sv(p)",
        grid=GridLevel.LOW,
        template=False,
        gas_phase=True,
    )


@pytest.fixture
def create_mock_conformers():
    """Factory fixture to create specified number of mock conformers"""

    def _create_conformers(count: int):
        conformers = []
        for i in range(1, count + 1):
            mock = Mock(spec=MoleculeData)
            mock.name = f"CONF{i}"
            mock.charge = 0
            mock.unpaired = 0
            mock.mo_paths = {
                QmProg.TM: ["path/to/mos"],
                QmProg.ORCA: ["path/to/file.gbw"],
            }
            mock.geom = GeometryData(
                name=f"CONF{i}", xyz=["C 0.0 0.0 0.0\n", "H 1.0 0.0 0.0\n"]
            )
            conformers.append(mock)
        return conformers

    return _create_conformers


# ============= Mock Task Functions =============


# Define task functions at module level for pickling
def task_success(
    job: ParallelJob, config: Any, resources: ResourceMonitor
) -> tuple[SPResult, MetaData]:
    """Task that always succeeds"""
    with resources.occupy_cores(job.omp):
        time.sleep(0.1)  # Shorter sleep for faster tests
    return SPResult(mo_path=f"path/to/mo_{job.conf.name}", energy=-100.0), MetaData(
        conf_name=job.conf.name, success=True
    )


def task_fail(
    job: ParallelJob, config: Any, resources: ResourceMonitor
) -> tuple[None, MetaData]:
    """Task that always fails"""
    with resources.occupy_cores(job.omp):
        time.sleep(0.1)
    return None, MetaData(conf_name=job.conf.name, success=False, error="Test error")


def task_conditional(
    job: ParallelJob, config: Any, resources: ResourceMonitor
) -> tuple[Any, MetaData]:
    """Task that succeeds for even-numbered conformers and fails for odd-numbered ones"""
    conf_num = int(job.conf.name.replace("CONF", ""))
    if conf_num % 2 == 0:
        return task_success(job, config, resources)
    return task_fail(job, config, resources)


# ============= Tests for Core Functionality =============


class TestCoreParallelComponents:
    """Tests for core parallel processing components"""

    def test_setup_managers(self):
        """Test setup_managers context manager creates correct manager instances"""
        max_workers = 2
        ncores = 4

        with setup_managers(max_workers, ncores) as (
            executor,
            manager,
            resource_manager,
        ):
            assert isinstance(executor, ProcessPoolExecutor)
            assert isinstance(manager, SyncManager)
            assert isinstance(resource_manager, ResourceMonitor)

    def test_resource_monitor(self):
        """Test ResourceMonitor core management functionality"""
        with setup_managers(2, 4) as (_, manager, resource_monitor):
            # Test initial state
            assert resource_monitor._ResourceMonitor__free_cores.value == 4

            # Test core occupation
            with resource_monitor.occupy_cores(2):
                assert resource_monitor._ResourceMonitor__free_cores.value == 2

            # Verify cores are released after context exit
            assert resource_monitor._ResourceMonitor__free_cores.value == 4


class TestResultDataclasses:
    """Tests for QM result dataclasses"""

    def test_sp_result(self):
        """Test SPResult dataclass"""
        sp_result = SPResult(mo_path="path/to/mo", energy=-100.0)
        assert sp_result.mo_path == "path/to/mo"
        assert sp_result.energy == -100.0

    def test_gsolv_result(self):
        """Test GsolvResult dataclass"""
        gsolv_result = GsolvResult(
            mo_path="path/to/mo", gsolv=-0.1, energy_gas=-100.0, energy_solv=-100.1
        )
        assert gsolv_result.mo_path == "path/to/mo"
        assert gsolv_result.gsolv == -0.1
        assert gsolv_result.energy_gas == -100.0
        assert gsolv_result.energy_solv == -100.1

    def test_rrho_result(self):
        """Test RRHOResult dataclass"""
        rrho_result = RRHOResult(
            mo_path="path/to/mo",
            energy=-100.0,
            rmsd=0.001,
            gibbs={298.15: -99.9},
            symmetry="c2v",
        )
        assert rrho_result.mo_path == "path/to/mo"
        assert rrho_result.energy == -100.0
        assert rrho_result.rmsd == 0.001
        assert rrho_result.gibbs[298.15] == -99.9
        assert rrho_result.symmetry == "c2v"

    def test_opt_result(self):
        """Test OptResult dataclass"""
        opt_result = OptResult(
            mo_path="path/to/mo",
            energy=-100.0,
            ecyc=[-99.0, -99.5, -100.0],
            converged=True,
        )
        assert opt_result.mo_path == "path/to/mo"
        assert opt_result.energy == -100.0
        assert len(opt_result.ecyc) == 3
        assert opt_result.converged

    def test_nmr_result(self):
        """Test NMRResult dataclass"""
        nmr_result = NMRResult(
            mo_path="path/to/mo", shieldings=[(1, 30.5)], couplings=[((1, 2), 7.2)]
        )
        assert nmr_result.mo_path == "path/to/mo"
        assert len(nmr_result.shieldings) == 1
        assert len(nmr_result.couplings) == 1
        assert nmr_result.shieldings[0] == (1, 30.5)
        assert nmr_result.couplings[0] == ((1, 2), 7.2)

    def test_uvvis_result(self):
        """Test UVVisResult dataclass"""
        uvvis_result = UVVisResult(
            mo_path="path/to/mo", excitations=[{"wavelength": 350.0, "f": 0.5}]
        )
        assert uvvis_result.mo_path == "path/to/mo"
        assert len(uvvis_result.excitations) == 1
        assert uvvis_result.excitations[0]["wavelength"] == 350.0
        assert uvvis_result.excitations[0]["f"] == 0.5


class TestMetaData:
    """Tests for MetaData class"""

    def test_metadata_success(self):
        """Test MetaData with successful execution"""
        meta = MetaData(conf_name="CONF1", success=True)
        assert meta.conf_name == "CONF1"
        assert meta.success
        assert meta.error == ""

    def test_metadata_failure(self):
        """Test MetaData with failed execution"""
        meta_failed = MetaData(conf_name="CONF2", success=False, error="Test error")
        assert meta_failed.conf_name == "CONF2"
        assert not meta_failed.success
        assert meta_failed.error == "Test error"


class TestParallelJob:
    """Tests for ParallelJob class"""

    def test_parallel_job_initialization(self, mock_molecule_data):
        """Test ParallelJob initialization and properties"""
        job = ParallelJob(mock_molecule_data.geom, 0, 0)
        assert job.conf == mock_molecule_data.geom
        assert job.omp == OMPMIN
        assert job.charge == 0
        assert job.unpaired == 0
        assert isinstance(job.flags, dict)


# ============= Tests for OMP Threading Configuration =============


class TestOMPConfiguration:
    """Tests for OMP threading configuration functions"""

    @pytest.mark.parametrize(
        "balance,omp,ncores,expected_omp",
        [
            (
                True,
                2,
                4,
                4,
            ),  # Test with balancing enabled, len(jobs) < ncores // OMPMIN (4 // 1 = 4), expect ncores // len(jobs) (4)
            (
                False,
                0,
                4,
                1,
            ),  # Test with balancing disabled, omp < OMPMIN, expect OMPMIN
            (
                False,
                5,
                4,
                4,
            ),  # Test with balancing disabled, omp > ncores, expect ncores
            (False, 5, 8, 5),  # Test with ncores > omp > OMPMIN, expect omp
        ],
    )
    def test_set_omp_tmproc(
        self, balance, omp, ncores, expected_omp, mock_parallel_job
    ):
        """Test set_omp_tmproc function handles thread allocation correctly"""
        jobs = [mock_parallel_job]
        set_omp_tmproc(jobs, balance, omp, ncores)
        assert jobs[0].omp == expected_omp

    @pytest.mark.parametrize(
        "balance,omp,ncores,njobs,expected_omp",
        [
            (True, 2, 4, 2, 2),  # Test with balancing enabled
            (False, 3, 4, 2, 3),  # Test with balancing disabled
            (True, 2, 8, 2, 4),  # Test optimal distribution
        ],
    )
    def test_set_omp(
        self, balance, omp, ncores, njobs, expected_omp, mock_parallel_job
    ):
        """Test set_omp function distributes threads properly"""
        jobs = [mock_parallel_job for _ in range(njobs)]
        set_omp(jobs, balance, omp, ncores)
        assert all(job.omp == expected_omp for job in jobs)


# ============= Tests for Job Preparation =============


class TestJobPreparation:
    """Tests for job preparation functions"""

    @pytest.mark.parametrize(
        "prog,copy_mo,expected_mo_guess",
        [
            (QmProg.ORCA, True, "path/to/file.gbw"),
            (QmProg.TM, True, "path/to/mos"),
            (QmProg.ORCA, False, None),
        ],
    )
    def test_prepare_jobs(self, prog, copy_mo, expected_mo_guess, mock_molecule_data):
        """Test prepare_jobs function creates jobs with correct parameters"""
        conformers = [mock_molecule_data]
        jobs = prepare_jobs(
            conformers, prog, 2, 4, "test", balance=True, copy_mo=copy_mo
        )

        assert len(jobs) == 1
        assert isinstance(jobs[0], ParallelJob)
        if copy_mo:
            assert str(jobs[0].mo_guess) == expected_mo_guess
        else:
            assert jobs[0].mo_guess is None


# ============= Tests for Job Execution =============


class TestJobExecution:
    """Tests for job execution functions"""

    def test_execute_success(self, mock_job_config, create_mock_conformers):
        """Test successful parallel execution of all jobs"""
        # Create two mock conformers
        conformers = create_mock_conformers(2)

        # Run execute
        results = execute(
            conformers=conformers,
            task=task_success,
            job_config=mock_job_config,
            prog=QmProg.ORCA,
            ncores=4,
            omp=2,
            from_part="test",
        )

        # Verify results
        assert len(results) == 2
        assert all(f"CONF{i}" in results for i in range(1, 3))

        # Check result contents
        for conf_name, result in results.items():
            assert isinstance(result, SPResult)
            assert result.mo_path == f"path/to/mo_{conf_name}"
            assert result.energy == -100.0

    def test_execute_all_failed(self, mock_job_config, create_mock_conformers):
        """Test execution when all jobs fail"""
        # Create two mock conformers
        conformers = create_mock_conformers(2)

        # Check that RuntimeError is raised when all jobs fail
        with pytest.raises(RuntimeError, match="All jobs failed to execute"):
            execute(
                conformers=conformers,
                task=task_fail,
                job_config=mock_job_config,
                prog=QmProg.ORCA,
                ncores=4,
                omp=2,
                from_part="test",
            )

    def test_execute_partial_failure(self, mock_job_config, create_mock_conformers):
        """Test execution with mixed success/failure"""
        # Create four mock conformers
        conformers = create_mock_conformers(4)

        # Run execute
        results = execute(
            conformers=conformers,
            task=task_conditional,
            job_config=mock_job_config,
            prog=QmProg.ORCA,
            ncores=4,
            omp=2,
            from_part="test",
            ignore_failed=True,
        )

        # Verify results (even-numbered conformers succeed, odd ones fail)
        assert len(results) == 2
        assert all(f"CONF{i}" in results for i in [2, 4])

        # Test ignore_failed
        with pytest.raises(RuntimeError, match="failed jobs"):
            results = execute(
                conformers=conformers,
                task=task_conditional,
                job_config=mock_job_config,
                prog=QmProg.ORCA,
                ncores=4,
                omp=2,
                from_part="test",
                ignore_failed=False,
            )

    def test_execute_resource_management(self, mock_job_config, create_mock_conformers):
        """Test resource management with multiple jobs"""
        # Create eight mock conformers
        n_jobs = 8
        conformers = create_mock_conformers(n_jobs)

        # Run execute with limited cores
        start_time = time.time()
        results = execute(
            conformers=conformers,
            task=task_success,
            job_config=mock_job_config,
            prog=QmProg.ORCA,
            ncores=4,  # Limit to 4 cores
            omp=2,  # 2 cores per job
            from_part="test",
            ignore_failed=False,
            balance=False,
            copy_mo=False,
        )
        total_time = time.time() - start_time

        # Verify results
        assert len(results) == n_jobs

        # With 4 cores total and 2 cores per job, we should be able to run 2 jobs in parallel
        # Each job takes 0.1s, so 8 jobs should take ~0.4s (plus overhead)
        # Allow for some timing variability in CI environments
        assert (
            0.3 < total_time < 2.0
        ), f"Expected ~0.4s execution time, got {total_time}s"
