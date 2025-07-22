import logging
import sys
from pathlib import Path

__loglevel = logging.INFO

__loggers = []

# _loglevel = logging.DEBUG


def setup_logger(name: str) -> logging.Logger:
    """
    Initializes and configures a logger with the specified name.

    Args:
        name (str): The name of the logger.

    Returns:
        logging.Logger: The configured logger instance.
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

    if logger not in __loggers:
        __loggers.append(logger)

    return logger


def set_filehandler(path: str | Path):
    """Set filehandler for all loggers."""
    handler = logging.FileHandler(path)
    formatter = logging.Formatter(
        "{asctime:24s}-{name:^24s}-{levelname:^10s}- {message}", style="{"
    )
    handler.setFormatter(formatter)
    for logger in __loggers:
        logger.addHandler(handler)


def set_loglevel(loglevel: str) -> None:
    """
    Set the log level for the logger.

    Args:
        loglevel (str): The log level to set.

    Returns:
        None
    """
    global __loglevel
    __loglevel = getattr(logging, loglevel)
    for logger in __loggers:
        logger.setLevel(__loglevel)
