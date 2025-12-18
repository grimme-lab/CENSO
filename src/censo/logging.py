import logging
import sys
from pathlib import Path

__loglevel = logging.INFO
__filehandler_path: str | Path | None = None

# _loglevel = logging.DEBUG


def setup_logger(name: str) -> logging.Logger:
    """
    Initializes and configures a logger with the specified name.

    If a file handler path has been configured via set_filehandler(),
    the logger will automatically receive a FileHandler in addition to
    the StreamHandler.

    :param name: The name of the logger.
    :return: The configured logger instance.
    """
    global __loglevel, __filehandler_path

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

        # Add FileHandler if one has been configured
        if __filehandler_path is not None:
            file_formatter = logging.Formatter(
                "{asctime:24s}-{name:^24s}-{levelname:^10s}- {message}", style="{"
            )
            file_handler = logging.FileHandler(__filehandler_path)
            file_handler.setLevel(__loglevel)
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)

    return logger


def set_filehandler(path: str | Path):
    """
    Set filehandler for all censo loggers, avoiding duplicates.

    This function stores the file handler path globally so that any loggers
    created after this call will automatically receive a FileHandler. It also
    updates all existing censo loggers to add the FileHandler if they don't
    already have one for this path.

    Multiple calls with different paths are supported - this will update all
    loggers to use the new path.

    :param path: Path to the log file.
    :return: None
    """
    global __loglevel, __filehandler_path
    __filehandler_path = path
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
                handler.setLevel(__loglevel)
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
