"""
Utility functions which are used in the CENSO modules. From creating folders to
printout routines.
"""

import functools
import hashlib
import os
import time
import re
from collections import OrderedDict
from collections.abc import Callable, Sequence
from pydantic import BaseModel
from typing import Any
import math

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
    data: dict[str, dict[str, float]] = {}
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


def average(x: list[int | float]):
    assert len(x) > 0
    return float(sum(x)) / len(x)


def pearson_def(x: list[int | float], y: list[int | float]):
    # Pad with last value
    if len(x) > len(y):
        while len(x) > len(y):
            y.append(y[-1])
    elif len(x) < len(y):
        while len(x) < len(y):
            x.append(x[-1])

    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    try:
        return diffprod / math.sqrt(xdiff2 * ydiff2)
    except ZeroDivisionError:
        return 1.0


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


def do_md5(path):
    """
    Calculate md5 of file to identifly if restart happend on the same file!
    Input is buffered into smaller sizes to ease on memory consumption.
    Hashes entire content of ensemble input file to compare later
    """
    BUF_SIZE = 65536
    md5 = hashlib.md5()
    if os.path.isfile(path):
        with open(path, "rb") as f:
            while True:
                data = f.read(BUF_SIZE)
                if not data:
                    break
                md5.update(data)
        return md5.hexdigest()
    else:
        raise FileNotFoundError


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


def od_insert(
    od: OrderedDict[str, any], key: str, value: any, index: int
) -> OrderedDict[str, any]:
    """
    Insert a new key/value pair into an OrderedDict at a specific position.
    If it was a normal dict:
        od[key] = value, with insertion before the 'index'th key.

    Args:
        od: The OrderedDict to insert into.
        key: The key to insert.
        value: The value associated with the key.
        index: The index before which to insert the key/value pair.

    Returns:
        The updated OrderedDict.
    """
    # FIXME - somehow this doesn't work reliably, no idea why but sometimes the value is not inserted
    items: list[tuple[str, any]] = list(od.items())
    items.insert(index, (key, value))
    return OrderedDict(items)


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
