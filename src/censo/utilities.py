"""
Utility functions which are used in the CENSO modules. From creating folders to
printout routines.
"""

import functools
import os
import time
from collections.abc import Callable
from pydantic import BaseModel
from typing import Any
from datetime import timedelta

from .params import BOHR2ANG, PLENGTH
from .logging import setup_logger
import json
from pydantic import ValidationError

logger = setup_logger(__name__)


def print_validation_errors(e: ValidationError) -> None:
    """
    Print Pydantic validation errors in a human-readable format.

    :param e: The ValidationError instance containing the errors.
    :return: None
    """
    printf(f"Found {e.error_count()} validation error(s):\n")
    printf(f"Found {e.error_count()} validation error(s):\n")
    for error in e.errors():
        field = " -> ".join(map(str, error["loc"]))
        message = error["msg"]
        user_input = error["input"]
        # Handle model-level validator errors differently
        if not error["loc"] or (
            len(error["loc"]) == 1 and error["loc"][0] == "__root__"
        ):
            printf("  - Model-level error:")
            printf(f"    Message: {message}")
        else:
            try:
                user_input_str = json.dumps(user_input)
            except TypeError:
                user_input_str = str(user_input)
            printf(f"  - Field: '{field}'")
            printf(f"    Message: {message}")
            printf(f"    Your input: {user_input_str}")
        printf("-" * 20)


class Factory[T]:
    """
    Abstract object factory class.
    """

    _registry: dict[str, Callable[..., T]] = {}

    @classmethod
    def register_builder(cls, name: str, builder: Callable[..., T]) -> None:
        """
        Register a builder.

        :param name: Name of the builder.
        :param builder: The builder callable.
        :return: None
        """
        cls._registry[name] = builder
        cls._registry[name] = builder

    @classmethod
    def create(cls, name: str, *args, **kwargs) -> T:
        """
        Generic factory method.

        :param name: Name of the builder to create.
        :return: The created instance.
        """
        if name not in cls._registry:
            raise TypeError(f"No type was found for '{name}' in {list(cls._registry)}.")
        builder = cls._registry[name]
        return builder(*args, **kwargs)


class DataDump(BaseModel):
    """
    Model for dumping data with part name and settings.
    """

    part_name: str
    data: dict[
        str,
        dict[
            str,
            Any,
        ],
    ] = {}
    settings: dict[str, Any] = {}


def printf(*args, **kwargs):
    """
    patch print to always flush
    """
    print(*args, flush=True, **kwargs)


def frange(start: float, end: float, step: float = 1) -> list[float]:
    """
    Create a range of floats, adding 'step' to 'start' while it's less or equal than 'end'.

    :param start: The start of the range.
    :param end: The end of the range.
    :param step: The step size. Defaults to 1.
    :return: The list of floats.
    """
    result = []
    current = start
    while current <= end:
        result.append(current)
        current += step
    return result


def t2x(
    path: str, writexyz: bool = False, outfile: str = "original.xyz"
) -> tuple[list[str], int, str]:
    """
    Convert TURBOMOLE coord file to xyz data and/or write *.xyz output.

    :param path: Absolute path either to dir or file directly.
    :param writexyz: If True, directly write to outfile. Defaults to False.
    :param outfile: Filename of xyz file which is written into the same directory as path.
    :return: Tuple of (coordxyz list of strings, number of atoms, output path).
    """
    # read lines from coord file
    with open(path, newline=None) as f:
        coord = f.readlines()

    # read coordinates with atom labels directly into a string
    # and append the string to a list to be written/returned later
    xyzatom = []
    for line in coord:
        if "$end" in line:  # stop at $end ...
            break
        xyzatom.append(
            functools.reduce(
                lambda x, y: x + " " + y,
                [
                    f"{float(line.split()[0]) * BOHR2ANG:.10f}",
                    f"{float(line.split()[1]) * BOHR2ANG:.10f}",
                    f"{float(line.split()[2]) * BOHR2ANG:.10f}",
                    f"{str(line.split()[3].lower()).capitalize()}",
                ],
            )
        )

    # get path from args without the filename of the ensemble (last element of path)
    if os.path.isfile(path):
        outpath = functools.reduce(
            lambda x, y: os.path.join(x, y), list(os.path.split(path))[::-1][1:][::-1]
        )
    # or just use the given path if it is not a file path
    else:
        outpath = path

    # write converted coordinates to xyz outfile if wanted
    if writexyz:
        with open(os.path.join(outpath, outfile), "w") as out:
            out.write(str(len(xyzatom)) + "\n")
            for line in xyzatom:
                out.write(line)
    return xyzatom, len(xyzatom), os.path.join(outpath, outfile)


def check_for_float(line: str) -> float | None:
    """
    Go through line and check for float, return first float.

    :param line: The line to check.
    :return: The first float found, or None.
    """
    elements = line.strip().split()
    value = None
    for element in elements:
        try:
            value = float(element)
        except ValueError:
            value = None
        if value is not None:
            break
    return value


def timeit(f: Callable[..., None]) -> Callable[..., float]:
    """
    Time function execution.

    .. deprecated::
        This decorator is deprecated. Timing should be handled at the CLI level.
        TODO: Remove in a follow-up change.

    :param f: The function to time (should have no return value).
    :return: A wrapper that returns the execution time.
    """

    @functools.wraps(f)
    def wrapper(*args, **kwargs) -> float:
        start = time.perf_counter()
        f(*args, **kwargs)
        end = time.perf_counter()
        return end - start

    return wrapper


def get_time(time: float) -> tuple[int, int, int]:
    """
    Calculate seconds, minutes, hours from time in seconds.

    :param time: The time in seconds.
    :return: Tuple of seconds, minutes, hours.
    """
    time_taken = timedelta(seconds=int(time))
    hours, r = divmod(time_taken.seconds, 3600)
    minutes, seconds = divmod(r, 60)
    if time_taken.days:
        hours += time_taken.days * 24
    return seconds, minutes, hours


def h1(text: str) -> str:
    """
    Create a formatted header of type 1.

    :param text: The text to be formatted.
    :return: The formatted header.
    """
    return "\n" + f" {text} ".center(PLENGTH, "-") + "\n"


def h2(text: str) -> str:
    """
    Create a formatted header of type 2.

    :param text: The text to be formatted.
    :return: The formatted header.
    """
    return f"""
{"-" * PLENGTH}
{text.center(PLENGTH, " ")}
{"-" * PLENGTH}
    """
