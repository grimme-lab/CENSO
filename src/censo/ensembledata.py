"""
stores ensembledata and conformers
functionality for program setup
"""

import os
import re
import json
from collections.abc import Callable
from math import exp

from .datastructure import MoleculeData
from .logging import setup_logger
from .params import Params
from .utilities import check_for_float, print, t2x

logger = setup_logger(__name__)


class EnsembleData:
    """ """

    def __init__(self, input_file: str = None):
        """
        Setup an EnsembleData object, which contains a list of conformers, read from
        input_file. If input_file is not passed here, conformers can be read using
        read_input.

        Args:
            input_file (str, optional): Path to the ensemble input file. Defaults to None.
            If this is provided, the charge and unpaired electron count will be assumed to be 0 and all conformers will be read from the input file.
        """
        # contains run-specific info that may change during runtime
        # initialized in EnsembleData.read_input
        self.runinfo = {
            "charge": None,
            "unpaired": None,
        }

        # stores the conformers with all info
        # NOTE: this is deliberately chosen to be a list since lists are ordered
        self.__conformers: list[MoleculeData] = []

        # stores the conformers which were sorted out
        self.rem: list[MoleculeData] = []

        if input_file is not None:
            self.read_input(input_file, charge=0, unpaired=0)

        # A list containing all parts run using this ensemble
        self.results: list = []

    @property
    def conformers(self):
        """
        Returns the conformers list. Includes a check wether there are any conformers left.
        """
        # TODO - no checks for now
        return self.__conformers

    @conformers.setter
    def conformers(self, confs):
        assert all(isinstance(conf, MoleculeData) for conf in confs)
        self.__conformers = confs

    def read_output(self, outpath: str) -> None:
        """
        Read json output file of a previous execution. Will try to load results into current conformer ensemble, matching
        based on names. If a conformer name does not exist in the current ensemble it will be ignored. If a conformer
        does not exist in the output data RuntimeError will be raised.

        Args:
            outpath (str): Path to the output file.

        Returns:
            None
        """

        with open(outpath, "r") as file:
            data = json.load(file)

        # Check if all conformers from the current ensemble are also found in the output data
        for partname, results in data.items():
            # First level of json output is the part name (loop iterates through all parts in the json file)
            # .values() are the results of each part
            if not all(conf.name in results for conf in self.conformers):
                raise RuntimeError(
                    "Not all conformers from the current ensemble are found in the output data."
                )

            # Update results dict for the conformers
            for conf in self.conformers:
                conf.results.setdefault(partname, {}).update(results[conf.name])

        logger.info(f"Reloaded results from {outpath}.")

    def read_input(
        self,
        input_path: str,
        charge: int = None,
        unpaired: int = None,
        nconf: int = None,
    ) -> None:
        """
        Read ensemble input file. Should be a file in xyz-file format with all the conformers in consecutive order.
        If command line arguments are given, the priorities of parameters will be:
            1. method arguments,
            2. cml arguments,
            3. defaults.

        Args:
            input_path (str): Path to the ensemble input file.
            charge (int, optional): Charge of the system. Defaults to None. Overwrites preexisting values.
            unpaired (int, optional): Number of unpaired electrons. Defaults to None. Overwrites preexisting values.
            nconf (int, optional): Number of conformers to consider. Defaults to None, so all conformers are read.

        Returns:
            None

        Raises:
            RuntimeError: If the charge or the number of unpaired electrons is not defined.
        """
        # If $coord in file => tm format, needs to be converted to xyz
        with open(input_path, "r") as inp:
            lines = inp.readlines()
            if any("$coord" in line for line in lines):
                _, nat, input_path = t2x(
                    input_path, writexyz=True, outfile="converted.xyz"
                )
            else:
                nat = int(lines[0].split()[0])

        # Set charge and unpaired via funtion args
        self.runinfo["charge"] = charge
        self.runinfo["unpaired"] = unpaired

        if self.runinfo["charge"] is None or self.runinfo["unpaired"] is None:
            raise RuntimeError("Charge or number of unpaired electrons not defined.")

        self.__setup_conformers(input_path, maxconf=nconf)

        # Print information about read ensemble
        print(
            f"Read {len(self.conformers)} conformers.\n",
            "Number of atoms:".ljust(Params.DIGILEN // 2, " ") + f"{nat}" + "\n",
            "Charge:".ljust(Params.DIGILEN // 2, " ")
            + f"{self.runinfo['charge']}"
            + "\n",
            "Unpaired electrons:".ljust(Params.DIGILEN // 2, " ")
            + f"{self.runinfo['unpaired']}"
            + "\n",
            sep="",
        )

    def __setup_conformers(self, input_path: str, maxconf: int = None) -> None:
        """
        open ensemble input
        split into conformers
        create MoleculeData objects out of coord input
        read out energy from xyz file if possible

        Args:
            input_path (str): Path to the ensemble input file.
            maxconf (int, optional): Maximum number of conformers to consider. Defaults to None, so all conformers are read.

        Returns:
            None
        """
        # open ensemble input
        with open(input_path, "r") as file:
            lines = file.readlines()

        # Get rid of unnecessary empty lines
        # Basically this filters out all only-whitespace lines except the comment lines after the number of atoms is declared
        lines = list(
            filter(
                lambda line: not (
                    bool(re.match(r"^\s*$", line))  # matches only whitespace chars
                    and len(lines[lines.index(line) - 1].split()) != 1
                ),
                lines,
            )
        )

        # assuming consecutive xyz-file format
        # (every conf geometry is separated by a line with split length of 4 followed by a line of split length 1)
        split_indices = [
            i
            for i in len(lines)
            if i == 0 or (len(lines[i].split()) == 1 and len(lines[i - 1].split()) == 4)
        ]

        # Check the number of conformers found in the input file with the max number of conformers to be read
        nconf = len(split_indices)
        if maxconf is not None and maxconf < len(split_indices):
            nconf = maxconf

        for i in range(nconf):
            # Check whether the names are stored in the ensemble file,
            # use those if possible because of crest rotamer files
            conf_index = split_indices[i]
            if "CONF" in lines[conf_index + 1]:
                confname = next(s for s in lines[conf_index + 1].split() if "CONF" in s)
            else:
                # Start counting from 1
                confname = f"CONF{i + 1}"

            # Determine end of geometry definition for this conf
            # which is either the next conf definition or EOF
            conf_end_index = (
                split_indices[i + 1] if i + 1 < len(split_indices) else len(lines)
            )

            # Don't use the property here, instead access the attribute directly
            # Create a new conformer object and append it to the ensemble
            self.__conformers.append(
                MoleculeData(
                    confname,
                    lines[conf_index:conf_end_index],
                )
            )

            # get precalculated energies if possible
            # precalculated energy set to 0.0 if it cannot be found
            self.conformers[i].xtb_energy = (
                check_for_float(lines[conf_index + 1]) or 0.0
            )

            # also works if xtb_energy is None for some reason (None is put first)
            self.conformers.sort(key=lambda x: x.xtb_energy)

    def remove_conformers(self, confnames: list[str]) -> None:
        """
        Remove the conformers with the names listed in 'confnames' from further consideration.
        The removed conformers will be stored in self.rem.

        Args:
            confnames (list[str]): A list of conformer names.

        Returns:
            None
        """
        for confname in confnames:
            remove = next(c for c in self.conformers if c.name == confname)

            # pop item from conformers and insert this item at index 0 in rem
            self.rem.insert(0, self.conformers.pop(self.conformers.index(remove)))

            # Log removed conformers
            logger.debug(f"Removed {remove.name}.")

    def dump(self, filename: str) -> None:
        """
        dump the conformers to a file
        """
        with open(os.path.join(f"{os.getcwd()}", f"{filename}.xyz"), "w") as file:
            for conf in self.conformers:
                file.writelines(conf.geom.toxyz())
