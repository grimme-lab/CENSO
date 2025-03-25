# import pytest
# from unittest.mock import patch, MagicMock
# from src.censo.parallel import ParallelExecutor
# from src.censo.datastructure import MoleculeData, ParallelJob
# from src.censo.qm_processor import QmProc
# from concurrent.futures import ProcessPoolExecutor
#
#
# @pytest.fixture
# def mock_conformers():
#     return [MoleculeData(geom="geom1"), MoleculeData(geom="geom2")]
#
#
# @pytest.fixture
# def mock_prepinfo():
#     return {"info1": {"param1": "value1"}, "info2": {"param2": "value2"}}
#
#
# @pytest.fixture
# def mock_jobtype():
#     return ["jobtype1", "jobtype2"]
#
#
# @pytest.fixture
# def mock_processor():
#     processor = MagicMock(spec=QmProc)
#     processor.run = MagicMock(return_value="result")
#     return processor
#
#
# def test_prepare_jobs(mock_conformers, mock_prepinfo, mock_jobtype):
#     executor = ParallelExecutor(
#         mock_conformers, "workdir", "prog", mock_prepinfo, mock_jobtype
#     )
#     jobs = executor.prepare_jobs()
#     assert len(jobs) == len(mock_conformers)
#     for job in jobs:
#         assert job.prepinfo == mock_prepinfo
#
#
# @patch("src.censo.parallel.Factory.create")
# def test_execute(
#     mock_factory, mock_conformers, mock_prepinfo, mock_jobtype, mock_processor
# ):
#     mock_factory.return_value = mock_processor
#     executor = ParallelExecutor(
#         mock_conformers, "workdir", "prog", mock_prepinfo, mock_jobtype
#     )
#     with patch.object(
#         ProcessPoolExecutor, "submit", return_value=MagicMock(result=lambda: "result")
#     ) as mock_submit:
#         results = executor.execute()
#     assert len(results) == len(mock_conformers)
#
#
# def test_set_omp_chunking(mock_conformers, mock_prepinfo, mock_jobtype):
#     executor = ParallelExecutor(
#         mock_conformers, "workdir", "prog", mock_prepinfo, mock_jobtype
#     )
#     jobs = executor.prepare_jobs()
#     executor.__set_omp_chunking()
#     for job in jobs:
#         assert job.omp >= 1
#
#
# def test_retry_failed_jobs(
#     mock_conformers, mock_prepinfo, mock_jobtype, mock_processor
# ):
#     executor = ParallelExecutor(
#         mock_conformers, "workdir", "prog", mock_prepinfo, mock_jobtype
#     )
#     jobs = executor.prepare_jobs()
#     retried, failed_confs = executor.retry_failed_jobs(
#         jobs, mock_processor, balance=True
#     )
#     assert isinstance(retried, list)
#     assert isinstance(failed_confs, list)
#
#
# @patch("src.censo.parallel.multiprocessing.Manager")
# def test_dynamic_queue_processing(
#     mock_manager, mock_conformers, mock_prepinfo, mock_jobtype, mock_processor
# ):
#     executor = ParallelExecutor(
#         mock_conformers, "workdir", "prog", mock_prepinfo, mock_jobtype
#     )
#     jobs = executor.prepare_jobs()
#     results = executor.dqp(jobs, mock_processor)
#     assert len(results) == len(jobs)
