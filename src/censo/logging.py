import os
import logging
import sys

__logpath: str = os.path.join(os.getcwd(), "censo.log")
__loglevel = logging.INFO

__loggers = []

# _loglevel = logging.DEBUG


def setup_logger(name: str, silent: bool = True) -> logging.Logger:
    """
    Initializes and configures a logger with the specified name.

    Args:
        name (str): The name of the logger.
        silent (bool, optional): Whether to print logpath or not. Defaults to True.

    Returns:
        logging.Logger: The configured logger instance.
    """
    if not silent:
        print(f"LOGFILE CAN BE FOUND AT: {__logpath}")

    # Create a logger instance with the specified name
    logger = logging.getLogger(name)
    logger.setLevel(__loglevel)

    # Create a FileHandler to log messages to the logpath file
    handler = logging.FileHandler(__logpath)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.WARNING)

    # Define the log message format
    formatter = logging.Formatter(
        "{asctime:24s}-{name:^24s}-{levelname:^10s}- {message}", style="{"
    )
    stream_formatter = logging.Formatter("{levelname:^10s}- {message}", style="{")
    handler.setFormatter(formatter)
    stream_handler.setFormatter(stream_formatter)

    # Add the FileHandler and StreamHandler to the logger
    logger.addHandler(handler)
    logger.addHandler(stream_handler)

    __loggers.append(logger)

    return logger


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
