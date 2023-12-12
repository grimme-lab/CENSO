"""
Utility functions which are used in the CENSO modules. From creating folders to
printout routines.
"""
import functools
import hashlib
import json
import logging
import os
import sys
import time
from builtins import print as print_orig
from collections import OrderedDict
from typing import Any, Callable, Tuple, Union, List, Dict

from .params import CODING, BOHR2ANG

__logpath: str = os.path.join(os.getcwd(), "censo.log")
__loglevel = logging.INFO


class DfaHelper:
    __dfa_dict: Dict

    @classmethod
    def set_dfa_dict(cls, dfadict_path: str):
        with open(dfadict_path, "r") as f:
            cls.__dfa_dict = json.load(f)

    @classmethod
    def find_func(cls, part: str, prog=None):
        """
        return all functionals available for a certain part and (optionally) program
        """
        # TODO - turn into filter using filterfunction defined within find_func
        tmp = []
        for k, v in cls.__dfa_dict["functionals"].items():
            if part in v["part"]:
                if prog is None:
                    tmp.append(k)
                else:
                    if v[prog] != "":
                        tmp.append(k)

        return tmp

    @classmethod
    def get_name(cls, func: str, prog: str):
        """
        return the name of a certain functional in the given qm program
        """
        return cls.__dfa_dict["functionals"][func][prog]

    @classmethod
    def get_disp(cls, func: str):
        """
        return the dispersion correction of a certain functional
        """
        return cls.__dfa_dict["functionals"][func]['disp']

    @classmethod
    def get_type(cls, func: str):
        """
        return the type of a certain functional
        """
        return cls.__dfa_dict["functionals"][func]["type"]

    @classmethod
    def functionals(cls) -> Dict[str, Dict]:
        return cls.__dfa_dict["functionals"]


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


def format_data(headers: List[str], rows: List[List[Any]], units: List[str] = None, sortby: int = 0) -> List[str]:
    """
    Generates a formatted table based on the given headers, rows, units, and sortby index.

    Args:
        headers (List[str]): The list of column headers.
        rows (List[List[Any]]): The list of rows, where each row is a list of values.
        units (List[str], optional): The list of units for each column. Defaults to None.
        sortby (int, optional): The index of the column to sort by. Defaults to 0.

    Returns:
        List[str]: The list of formatted lines representing the table.

    """
    lines = []

    # determine column width 'collen' of column with header 'header' 
    # by finding the length of the maximum width entry
    # for each column (header)
    collens = {
        header: collen for header, collen in zip(
            headers,
            (max(len(header), max(len(row) for row in rows)) for header in headers)
        )
    }

    # add table header
    lines.append(" ".join(f"{header:^{collen + 6}}" for header, collen in collens.items()))
    lines[0] += "\n"
    if units is not None:
        lines.append(" ".join(f"{unit:^{collen + 6}}" for unit, collen in zip(units, collens.values())))
        lines[1] += "\n"

    # TODO - draw an arrow if conformer is the best in current ranking
    # ("    <------\n" if self.key(conf) == self.key(self.core.conformers[0]) else "\n")

    # add a line for every row, sorted by the 'sortby'th column
    for row in sorted(rows, key=lambda x: x[sortby]):
        lines.append(" ".join(f"{row:^{collen + 6}}" for row, collen in zip(row, collens.values())))
        lines[-1] += "\n"

    return lines


def frange(start, end, step=1):
    """
    range with floats
    """
    result = []
    current = start
    while current < end:
        result.append(current)
        current += step
    return result


def t2x(path: str, writexyz: bool = False, outfile: str = "original.xyz") -> Tuple[list, int, str]:
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
    with open(path, "r", encoding=CODING, newline=None) as f:
        coord = f.readlines()

    # read coordinates with atom labels directly into a string
    # and append the string to a list to be written/returned later
    xyzatom = []
    for line in coord:
        if "$end" in line:  # stop at $end ...
            break
        xyzatom.append(functools.reduce(
            lambda x, y: x + " " + y,
            [
                f"{float(line.split()[0]) * BOHR2ANG:.10f}",
                f"{float(line.split()[1]) * BOHR2ANG:.10f}",
                f"{float(line.split()[2]) * BOHR2ANG:.10f}",
                f"{str(line.split()[3].lower()).capitalize()}",
            ]
        )
        )

    # get path from args without the filename of the ensemble (last element of path)
    if os.path.isfile(path):
        outpath = functools.reduce(
            lambda x, y: os.path.join(x, y),
            list(os.path.split(path))[::-1][1:][::-1]
        )
    # or just use the given path if it is not a file path
    else:
        outpath = path

    # write converted coordinates to xyz outfile if wanted
    if writexyz:
        with open(os.path.join(outpath, outfile), "w", encoding=CODING) as out:
            out.write(str(len(xyzatom)) + "\n")
            for line in xyzatom:
                out.write(line)
    return xyzatom, len(xyzatom), os.path.join(outpath, outfile)


def check_for_float(line: str) -> Union[float, None]:
    """ Go through line and check for float, return first float"""
    elements = line.strip().split()
    value = None
    for element in elements:
        try:
            value = float(element)
        except ValueError:
            value = None
        if value:
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


def od_insert(od: OrderedDict[str, Any], key: str, value: Any, index: int) -> OrderedDict[str, Any]:
    """
    Insert a new key/value pair into an OrderedDict at a specific position.

    Args:
        od: The OrderedDict to insert into.
        key: The key to insert.
        value: The value associated with the key.
        index: The index before which to insert the key/value pair.

    Returns:
        The updated OrderedDict.
    """
    items: List[Tuple[str, Any]] = list(od.items())
    items.insert(index, (key, value))
    return OrderedDict(items)


def setup_logger(name: str, silent: bool = True) -> logging.Logger:
    """
    Initializes and configures a logger with the specified name.

    Args:
        name (str): The name of the logger.
        silent (bool, optional): Whether to print logpath or not. Defaults to True.

    Returns:
        logging.Logger: The configured logger instance.
    """
    global __logpath, __loglevel

    if not silent:
        print(f"LOGFILE CAN BE FOUND AT: {__logpath}")

    # Create a logger instance with the specified name
    logger = logging.getLogger(name)
    logger.setLevel(__loglevel)

    # Create a FileHandler to log messages to the logpath file
    handler = logging.FileHandler(__logpath)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)

    # Define the log message format
    formatter = logging.Formatter('{asctime:24s}-{name:^20s}-{levelname:^10s}- {message}', style="{")
    stream_formatter = logging.Formatter('{levelname:^10s}- {message}', style="{")
    handler.setFormatter(formatter)
    stream_handler.setFormatter(stream_formatter)

    # Add the FileHandler and StreamHandler to the logger
    logger.addHandler(handler)
    logger.addHandler(stream_handler)

    return logger