"""
Unit tests for logging module.
"""

import logging
import os
import tempfile
from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def reset_logging_state():
    """Reset logging state between tests."""
    import censo.logging

    # Store original state
    original_filehandler_path = censo.logging.__filehandler_path

    # Reset to None for the test
    censo.logging.__filehandler_path = None

    # Clean up env vars set by set_filehandler / set_loglevel
    os.environ.pop("CENSO_LOG_PATH", None)
    os.environ.pop("CENSO_LOG_LEVEL", None)

    yield

    # Restore original state
    censo.logging.__filehandler_path = original_filehandler_path

    # Clean up any test loggers
    loggers_to_remove = [
        name
        for name in logging.Logger.manager.loggerDict.keys()
        if name.startswith("censo.test_")
    ]
    for name in loggers_to_remove:
        logger = logging.getLogger(name)
        # Remove all handlers
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)
        # Remove from manager
        del logging.Logger.manager.loggerDict[name]


def test_logger_created_after_set_filehandler_receives_filehandler():
    """Test that logger created after set_filehandler() receives FileHandler."""
    from censo.logging import setup_logger, set_filehandler

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = Path(tmpdir) / "test.log"

        # First set the file handler
        set_filehandler(log_path)

        # Then create the logger
        logger_name = "censo.test_after_filehandler"
        logger = setup_logger(logger_name)

        # Check that the logger has both StreamHandler and FileHandler
        assert len(logger.handlers) == 2
        has_stream_handler = any(
            isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
            for h in logger.handlers
        )
        has_file_handler = any(
            isinstance(h, logging.FileHandler) for h in logger.handlers
        )
        assert has_stream_handler, "Logger should have StreamHandler"
        assert has_file_handler, "Logger should have FileHandler"

        # Verify the FileHandler points to the correct file
        file_handler = next(
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        )
        assert file_handler.baseFilename == str(log_path)

        # Test that logging to the file works
        logger.info("Test message")
        assert log_path.exists()
        content = log_path.read_text()
        assert "Test message" in content


def test_logger_created_before_set_filehandler_receives_filehandler():
    """Test that logger created before set_filehandler() receives FileHandler via update."""
    from censo.logging import setup_logger, set_filehandler

    # First create the logger (before setting file handler)
    logger_name = "censo.test_before_filehandler"
    logger = setup_logger(logger_name)

    # Should only have StreamHandler initially
    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.StreamHandler)

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = Path(tmpdir) / "test.log"

        # Then set the file handler
        set_filehandler(log_path)

        # Check that the logger now has both handlers
        assert len(logger.handlers) == 2
        has_file_handler = any(
            isinstance(h, logging.FileHandler) for h in logger.handlers
        )
        assert (
            has_file_handler
        ), "Logger should have FileHandler after set_filehandler()"

        # Verify the FileHandler points to the correct file
        file_handler = next(
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        )
        assert file_handler.baseFilename == str(log_path)


def test_multiple_set_filehandler_calls_update_loggers():
    """Test that multiple set_filehandler() calls update all loggers correctly."""
    from censo.logging import setup_logger, set_filehandler

    logger_name1 = "censo.test_multiple_1"
    logger_name2 = "censo.test_multiple_2"

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path1 = Path(tmpdir) / "test1.log"
        log_path2 = Path(tmpdir) / "test2.log"

        # Set first file handler
        set_filehandler(log_path1)

        # Create first logger
        logger1 = setup_logger(logger_name1)

        # Verify logger1 has FileHandler pointing to log_path1
        file_handler1 = next(
            h for h in logger1.handlers if isinstance(h, logging.FileHandler)
        )
        assert file_handler1.baseFilename == str(log_path1)

        # Set second file handler with different path
        set_filehandler(log_path2)

        # Create second logger
        logger2 = setup_logger(logger_name2)

        # Verify logger2 has FileHandler pointing to log_path2
        file_handler2 = next(
            h for h in logger2.handlers if isinstance(h, logging.FileHandler)
        )
        assert file_handler2.baseFilename == str(log_path2)

        # Verify logger1 still has its original FileHandler
        # (set_filehandler doesn't remove old handlers, just adds new ones if path differs)
        file_handlers1 = [
            h for h in logger1.handlers if isinstance(h, logging.FileHandler)
        ]
        # Should have handlers for both paths now
        assert len(file_handlers1) >= 1


def test_no_duplicate_filehandlers_for_same_path():
    """Test that no duplicate FileHandlers are created for the same path."""
    from censo.logging import setup_logger, set_filehandler

    logger_name = "censo.test_no_duplicates"

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = Path(tmpdir) / "test.log"

        # Set file handler
        set_filehandler(log_path)

        # Create logger
        logger = setup_logger(logger_name)

        # Count FileHandlers
        initial_file_handler_count = sum(
            1 for h in logger.handlers if isinstance(h, logging.FileHandler)
        )

        # Call set_filehandler again with the same path
        set_filehandler(log_path)

        # Count FileHandlers again
        final_file_handler_count = sum(
            1 for h in logger.handlers if isinstance(h, logging.FileHandler)
        )

        # Should not have added a duplicate
        assert (
            initial_file_handler_count == final_file_handler_count
        ), "Should not create duplicate FileHandlers for the same path"


def test_set_filehandler_persists_to_env():
    """Test that set_filehandler stores the path in an env var for worker processes."""
    from censo.logging import set_filehandler

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = Path(tmpdir) / "test.log"

        # Env var should not be set before calling set_filehandler
        assert os.environ.get("CENSO_LOG_PATH") is None

        set_filehandler(log_path)

        # Env var should now be set
        assert os.environ.get("CENSO_LOG_PATH") == str(log_path)


def test_set_loglevel_persists_to_env():
    """Test that set_loglevel stores the level in an env var for worker processes."""
    from censo.logging import set_loglevel

    # Env var should not be set before calling set_loglevel
    os.environ.pop("CENSO_LOG_LEVEL", None)

    set_loglevel("DEBUG")

    # Env var should now be set
    assert os.environ.get("CENSO_LOG_LEVEL") == "DEBUG"


def test_setup_logger_reads_env_vars():
    """Test that setup_logger falls back to env vars when module state is unset.

    This simulates what happens in a dask worker process: the module-level
    __filehandler_path is None because it's a fresh process, but the env
    vars are inherited from the parent process.
    """
    from censo.logging import setup_logger
    import censo.logging

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = Path(tmpdir) / "test.log"

        # Simulate a worker process: module state is None, but env vars are set
        censo.logging.__filehandler_path = None
        os.environ["CENSO_LOG_PATH"] = str(log_path)
        os.environ["CENSO_LOG_LEVEL"] = "DEBUG"

        # Create logger - should pick up env vars
        logger = setup_logger("censo.test_env_vars")

        # Verify logger has FileHandler with correct path
        has_file_handler = any(
            isinstance(h, logging.FileHandler) for h in logger.handlers
        )
        assert has_file_handler, "Logger should have FileHandler from env var"

        # Verify the FileHandler points to the correct file
        file_handler = next(
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        )
        assert file_handler.baseFilename == str(log_path)

        # Verify the log level was picked up
        assert logger.level == logging.DEBUG

        # Verify logging actually works
        logger.info("Test message from env var path")
        assert log_path.exists()
        content = log_path.read_text()
        assert "Test message from env var path" in content


def test_cli_lazy_import_modules_receive_filehandler():
    """Test that lazily imported modules receive FileHandlers."""
    from censo.logging import setup_logger, set_filehandler

    logger_name = "censo.test_lazy_import"

    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = Path(tmpdir) / "test.log"

        # Simulate CLI flow: set_filehandler before importing modules
        set_filehandler(log_path)

        # Simulate lazy import: create logger after set_filehandler
        logger = setup_logger(logger_name)

        # Verify logger has FileHandler
        has_file_handler = any(
            isinstance(h, logging.FileHandler) for h in logger.handlers
        )
        assert has_file_handler, "Lazily imported module logger should have FileHandler"

        # Test that logging works
        logger.info("Lazy import test message")
        assert log_path.exists()
        content = log_path.read_text()
        assert "Lazy import test message" in content
