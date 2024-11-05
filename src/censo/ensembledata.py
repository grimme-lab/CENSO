"""
stores ensembledata and conformers
functionality for program setup
"""

import os
import re
import json

from .datastructure import MoleculeData
from .logging import setup_logger
from .params import DIGILEN
from .utilities import check_for_float, print, t2x, Factory

logger = setup_logger(__name__)


class EnsembleData:
    """
    Class to store conformer rotamer ensembles for use in CENSO.
    """

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

        # A list containing all part references in order of execution or loading
        self.results = []

        if input_file is not None:
            self.read_input(input_file, charge=0, unpaired=0)

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
        Read json output file of a previous execution. Will try to load data into current conformer ensemble, matching
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
        if not all(conf.name in data["results"] for conf in self.conformers):
            raise RuntimeError(
                "Not all conformers from the current ensemble are found in the output data."
            )

        # Create a part instance and load in the results
        part = Factory.create(data["partname"], self)
        part.data.update(data)

        logger.info(f"Reloaded results from {outpath}.")

        self.results.append(part)

    def read_input(
        self,
        input_path: str,
        charge: int = None,
        unpaired: int = None,
        nconf: int = None,
        append: bool = False,
    ) -> None:
        """
        Read ensemble input file. Should be a file in xyz-file format with all the conformers in consecutive order.

        Args:
            input_path (str): Path to the ensemble input file.
            charge (int, optional): Charge of the system. Defaults to None. Overwrites preexisting values.
            unpaired (int, optional): Number of unpaired electrons. Defaults to None. Overwrites preexisting values.
            nconf (int, optional): Number of conformers to consider. Defaults to None, so all conformers are read.
            append (bool, optional): If True, the conformers will be appended to the existing ensemble. Defaults to False.

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

        confs = self.__setup_conformers(input_path)
        if len(confs) == 0:
            logger.warning("Input file is empty!")

        if nconf is None:
            nconf = len(confs)

        if append:
            self.conformers.append(confs[:nconf])
        else:
            self.conformers = confs[:nconf]

        try:
            self.conformers.sort(key=lambda x: x.xtb_energy)
        except TypeError:
            # Only sort if all conformers have a defined precalculated energy
            pass

        # Print information about read ensemble
        print(
            f"Read {len(self.conformers)} conformers.\n",
            "Number of atoms:".ljust(DIGILEN // 2, " ") + f"{nat}" + "\n",
            "Charge:".ljust(DIGILEN // 2, " ") + f"{self.runinfo['charge']}" + "\n",
            "Unpaired electrons:".ljust(DIGILEN // 2, " ")
            + f"{self.runinfo['unpaired']}"
            + "\n",
            sep="",
        )

    def __setup_conformers(self, input_path: str) -> list[MoleculeData]:
        """
        open ensemble input
        split into conformers
        create MoleculeData objects out of coord input
        read out energy from xyz file if possible
        In principle this can also read xyz-files with molecules of different sizes.

        Args:
            input_path (str): Path to the ensemble input file.

        Returns:
            list[MoleculeData]: A list of MoleculeData objects.
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
        #
        # 14                               <-- split_index refers to this line (this is line 0 for the first conf)
        # CONF12 -22.521386
        # H x.xxxxxxxx x.xxxxxxx x.xxxxxx
        # ...
        split_indices = [
            i
            for i in range(len(lines))
            if i == 0 or (len(lines[i].split()) == 1 and len(lines[i - 1].split()) == 4)
        ]

        conformers = []
        for i, split_index in enumerate(split_indices):
            # Check whether the names are stored in the ensemble file,
            # use those if possible because of crest rotamer files
            if "CONF" in lines[split_index + 1]:
                confname = next(
                    s for s in lines[split_index + 1].split() if "CONF" in s
                )
            else:
                # Start counting from 1
                confname = f"CONF{i + 1}"

            # Determine end of geometry definition for this conf
            # which is either the next conf definition or EOF
            conf_end_index = (
                split_indices[i + 1] if i + 1 < len(split_indices) else len(lines)
            )

            # Create a new conformer object and append it to the ensemble
            conformers.append(
                MoleculeData(
                    confname,
                    lines[split_index + 2 : conf_end_index],
                )
            )

            # get precalculated energies if possible
            # precalculated energy set to 0.0 if it cannot be found
            conformers[i].xtb_energy = check_for_float(lines[split_index + 1]) or 0.0

        return conformers

    def remove_conformers(self, confnames: list[str]) -> None:
        """
        Remove the conformers with the names listed in 'confnames' from further consideration.
        The removed conformers will be stored in self.rem.

        Args:
            confnames (list[str]): A list of conformer names.

        Returns:
            None
        """
        if len(confnames) > 0:
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
