"""
Utility functions which are used in the CENSO modules. From creating folders to
printout routines.
"""

import functools
import hashlib
import json
import os
import time
import re
from builtins import print as print_orig
from collections import OrderedDict
from collections.abc import Callable

from .params import BOHR2ANG, PLENGTH, Config
from .logging import setup_logger

logger = setup_logger(__name__)


class Factory:
    """
    Generic object factory class.
    """

    __builders: dict[str, type] = {}

    @classmethod
    def register_builder(cls, name: str, builder: type) -> None:
        """
        Registers a builder.

        Args:
            name (str): name of the builder.
            builder (type): type of the builder.
        """
        cls.__builders[name] = builder

    @classmethod
    def create(cls, name: str, *args, **kwargs) -> object:
        """
        Generic factory method
        """
        builder = cls.__builders.get(name, None)

        if builder is not None:
            return builder(*args, **kwargs)
        raise TypeError(f"No type was found for '{name}' in {list(cls.__builders)}.")


class DfaHelper:
    _dfa_dict: dict

    @classmethod
    def set_dfa_dict(cls, dfadict_path: str):
        with open(dfadict_path, "r") as f:
            cls._dfa_dict = json.load(f)

    @classmethod
    def get_funcs(cls, prog: str):
        """
        Returns all functionals available for a given qm program.

        Args:
            prog (str): The qm program name.

        Returns:
            list[str]: The list of functionals.
        """
        return [
            func
            for func, v in cls._dfa_dict["functionals"].items()
            if v[prog.lower()] is not None
        ]

    @classmethod
    def get_name(cls, func: str, prog: str):
        """
        Returns the name of a certain functional in the given qm program. If name could not
        be found, the string passed as func will be returned instead.

        Args:
            func (str): The functional.
            prog (str): The qm program.

        Returns:
            str: The name of the functional.
        """
        func = func.lower()
        prog = prog.lower()
        if func in cls._dfa_dict["functionals"].keys():
            name = cls._dfa_dict["functionals"][func][prog]
        else:
            logger.warning(
                f"Functional {func} not found for program {prog}. Applying name literally."
            )
            name = func
        return name

    @classmethod
    def get_disp(cls, func: str):
        """
        Returns the dispersion correction of a given functional. If dispersion correction
        cannot be determined, apply none.

        Args:
            func (str): The functional.

        Returns:
            str: The dispersion correction name.
        """
        func = func.lower()
        if func in cls._dfa_dict["functionals"].keys():
            disp = cls._dfa_dict["functionals"][func]["disp"]
        else:
            logger.warning(
                f"Could not determine dispersion correction for {func}. Applying none."
            )
            disp = "novdw"
        return disp

    @classmethod
    def get_type(cls, func: str):
        """
        Returns the type of a certain functional. If the type cannot be determined, it
        is assumed to be a GGA.

        Args:
            func (str): The functional.

        Returns:
            str: The type of the functional.
        """
        func = func.lower()
        if func in cls._dfa_dict["functionals"].keys():
            rettype = cls._dfa_dict["functionals"][func]["type"]
        else:
            logger.warning(
                f"Could not determine functional type for {func}. Assuming GGA."
            )
            rettype = "GGA"
        return rettype

    @classmethod
    def functionals(cls) -> dict[str, dict]:
        return cls._dfa_dict["functionals"]


class SolventHelper:
    """
    Helper class to manage solvent lookup.
    """

    @classmethod
    def set_solvent_dict(cls, solvent_dict_path: str) -> None:
        """
        Load the solvents lookup dict.

        Args:
            solvent_dict_path (str): The path to the solvents lookup dict.
        """
        with open(solvent_dict_path, "r") as f:
            cls._solv_dict = json.load(f)

    @classmethod
    def get_solvent(cls, sm: str, name: str) -> str | None:
        """
        Try to lookup the solvent model keyword for the given solvent name. If it is not found, return None.

        Args:
            sm (str): The solvent model.
            name (str): The solvent name.

        Returns:
            str | None: The solvent model keyword or None if not found.
        """
        available_solvent_names_dict = cls.get_solvents_dict(sm)
        if name not in available_solvent_names_dict:
            return None
        return available_solvent_names_dict[name]

    @classmethod
    def get_solvents_dict(cls, sm: str) -> dict:
        """
        Get all available solvent names for a specified solvent model with the respective internal keyword.

        Args:
            sm (str): The solvent model.

        Returns:
            dict: The solvent names mapping onto the solvent keyword in the model.
        """
        return {
            name: sm_keys[sm]
            for name, sm_keys in cls._solv_dict.items()
            if sm in sm_keys
        }


def print(*args, **kwargs):
    """
    patch print to always flush
    """
    sep = " "
    end = "\n"
    file = None
    flush = True
    for key, value in kwargs.items():
        if key == "sep":
            sep = value
        elif key == "end":
            end = value
        elif key == "file":
            file = value
        elif key == "flush":
            flush = value
    print_orig(*args, sep=sep, end=end, file=file, flush=flush)


def format_data(
    headers: list[str],
    rows: list[list[str]],
    units: list[str] = None,
    sortby: int = 0,
    padding: int = 6,
) -> list[str]:
    """
    Generates a formatted table based on the given headers, rows, units, and sortby index.

    Args:
        headers (list[str]): The list of column headers.
        rows (list[list[str]]): The list of rows, where each row is a list of values.
        units (list[str], optional): The list of units for each column. Defaults to None.
        sortby (int, optional): The index of the column to sort by. Defaults to 0. In case of a string column,
                                use natural sorting.

    Returns:
        list[str]: The list of formatted lines representing the table.

    """

    def natural_sort_key(s):
        """
        Natural sorting key for strings.
        """
        return [int(text) if text.isdigit() else text for text in re.split(r"(\d+)", s)]

    lines = []

    # First, determine the maximium width for each column
    ncols = len(headers)
    if units is not None:
        maxcolw = [
            max(
                [
                    len(headers[i]),
                    max(len(rows[j][i]) for j in range(len(rows))),
                    len(units[i]),
                ]
            )
            for i in range(ncols)
        ]
    else:
        maxcolw = [
            max(len(headers[i]), max(len(rows[j][i]) for j in range(len(rows))))
            for i in range(ncols)
        ]

    # add table header
    lines.append(
        " ".join(f"{headers[i]:^{width + padding}}" for i, width in enumerate(maxcolw))
        + "\n"
    )

    # Add units
    if units is not None:
        lines.append(
            " ".join(
                f"{units[i]:^{width + padding}}" for i, width in enumerate(maxcolw)
            )
            + "\n"
        )

    # TODO - draw an arrow if conformer is the best in current ranking
    # ("    <------\n" if self.key(conf) == self.key(self.core.conformers[0]) else "\n")

    # Sort rows lexicographically if column sorted by is a number
    if rows[0][sortby].replace(".", "", 1).isdigit():
        rows = sorted(rows, key=lambda x: x[sortby])
    # Otherwise use natural sorting
    else:
        rows = sorted(rows, key=lambda x: natural_sort_key(x[sortby]))

    # add a line for every row
    for row in rows:
        lines.append(
            " ".join(f"{row[i]:^{width + padding}}" for i, width in enumerate(maxcolw))
            + "\n"
        )

    # Remove leading whitespace
    start = min(len(line) - len(line.lstrip()) for line in lines)
    for i in range(len(lines)):
        lines[i] = lines[i][start:]

    return lines


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


def timeit(f) -> Callable:
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
