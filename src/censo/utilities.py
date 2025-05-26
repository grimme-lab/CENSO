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

from .params import BOHR2ANG, PLENGTH
from .logging import setup_logger

logger = setup_logger(__name__)


class Factory[T]:
    """
    Abstract object factory class.
    """

    _registry: dict[str, Callable[..., T]] = {}

    @classmethod
    def register_builder(cls, name: str, builder: Callable[..., T]) -> None:
        """
        Registers a builder.

        Args:
            name (str): name of the builder.
            builder (Callable[..., T]): type of the builder.
        """
        cls._registry[name] = builder

    @classmethod
    def create(cls, name: str, *args, **kwargs) -> T:
        """
        Generic factory method
        """
        if name not in cls._registry:
            raise TypeError(f"No type was found for '{name}' in {list(cls._registry)}.")
        builder = cls._registry[name]
        return builder(*args, **kwargs)


class DataDump(BaseModel):
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
    Creates a range of floats, adding 'step' to 'start' while it's less or equal than 'end'.

    Args:
        start (float): The start of the range.
        end (float): The end of the range.
        step (float, optional): The step size. Defaults to 1.

    Returns:
        list[float]: The list of floats.
    """
    result = []
    current = start
    while current <= end:
        result.append(current)
        current += step
    return result


def t2x(
    path: str, writexyz: bool = False, outfile: str = "original.xyz"
) -> tuple[list, int, str]:
    """
    convert TURBOMOLE coord file to xyz data and/or write *.xyz output

     - path [abs. path] either to dir or file directly
     - writexyz [bool] default=False, directly write to outfile
     - outfile [filename] default = 'original.xyz' filename of xyz file which
                        is written into the same directory as
     returns:
     - coordxyz --> list of strings including atom x y z information
     - number of atoms
    """
    # read lines from coord file
    with open(path, "r", encoding=Config.CODING, newline=None) as f:
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
        with open(os.path.join(outpath, outfile), "w", encoding=Config.CODING) as out:
            out.write(str(len(xyzatom)) + "\n")
            for line in xyzatom:
                out.write(line)
    return xyzatom, len(xyzatom), os.path.join(outpath, outfile)


def check_for_float(line: str) -> float | None:
    """Go through line and check for float, return first float"""
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
    time function execution
    timed function should have no return value, since it is lost in the process
    calling a decorated function returns the time spent for it's execution
    """

    @functools.wraps(f)
    def wrapper(*args, **kwargs) -> float:
        start = time.perf_counter()
        f(*args, **kwargs)
        end = time.perf_counter()
        return end - start

    return wrapper


def h1(text: str) -> str:
    """
    Creates a formatted header of type 1:
        ---- text ----

    Args:
        text: The text to be formatted.

    Returns:
        The formatted header.
    """
    return "\n" + f" {text} ".center(PLENGTH, "-") + "\n"


def h2(text: str) -> str:
    """
    Creates a formatted header of type 2:
        ----------
           text
        ----------

    Args:
        text: The text to be formatted.

    Returns:
        The formatted header.
    """
    return f"""
{'-' * PLENGTH}
{text.center(PLENGTH, " ")}
{'-' * PLENGTH}
    """
