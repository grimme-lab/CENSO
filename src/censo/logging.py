import logging
import sys
from pathlib import Path

__loglevel = logging.INFO

# _loglevel = logging.DEBUG


def setup_logger(name: str) -> logging.Logger:
    """
    Initializes and configures a logger with the specified name.

    :param name: The name of the logger.
    :return: The configured logger instance.
    """
    # Create a logger instance with the specified name
    logger = logging.getLogger(name)
    logger.setLevel(__loglevel)

    # Only add handlers if the logger doesn't have any
    if not logger.handlers:
        # Create a StreamHandler to log messages to stdout (only WARNING or above)
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(logging.WARNING)

        # Define the log message format
        stream_formatter = logging.Formatter("{levelname:<10s}- {message}", style="{")
        stream_handler.setFormatter(stream_formatter)

        # Add the StreamHandler to the logger
        logger.addHandler(stream_handler)

    return logger


def set_filehandler(path: str | Path):
    """
    Set filehandler for all censo loggers, avoiding duplicates.

    :param path: Path to the log file.
    :return: None
    """
    filehandler_path = str(path)
    formatter = logging.Formatter(
        "{asctime:24s}-{name:^24s}-{levelname:^10s}- {message}", style="{"
    )
    # Loop over all loggers registered in logging
    for logger_name, logger in logging.Logger.manager.loggerDict.items():
        if isinstance(logger, logging.Logger) and logger_name.startswith("censo"):
            # Check for pre-existing FileHandlers with same path
            filehandler_exists = any(
                isinstance(h, logging.FileHandler)
                and getattr(h, "baseFilename", None) == filehandler_path
                for h in logger.handlers
            )
            if not filehandler_exists:
                handler = logging.FileHandler(path)
                handler.setFormatter(formatter)
                logger.addHandler(handler)


def set_loglevel(loglevel: str) -> None:
    """
    Set the log level for all censo loggers.

    :param loglevel: The log level to set.
    :return: None
    """
    global __loglevel
    __loglevel = getattr(logging, loglevel)
    for logger_name, logger in logging.Logger.manager.loggerDict.items():
        if isinstance(logger, logging.Logger) and logger_name.startswith("censo"):
            logger.setLevel(__loglevel)
