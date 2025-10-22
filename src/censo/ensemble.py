"""
stores ensembledata and conformers
functionality for program setup
"""

from pathlib import Path
import re
import json
from typing import Callable, Iterator
from math import exp

from .molecules import MoleculeData, Contributions
from .logging import setup_logger
from .params import AU2J, KB
from .utilities import check_for_float, printf, t2x

logger = setup_logger(__name__)


class EnsembleData:
    """
    Class to store conformer rotamer ensembles for use in CENSO.
    """

    def __init__(self):
        """
        Setup an EnsembleData object. Conformers can be read using read_input.
        """
        # stores the conformers with all info
        self.__conformers: list[MoleculeData] = []

        # stores the conformers which were sorted out
        self.__rem: list[MoleculeData] = []

        # String containing the constraints in xtb format to be applied to geometry optimizations
        # TODO: and xtb_rrho?
        self.constraints: str | None = None

    def __iter__(self) -> Iterator[MoleculeData]:
        return iter(self.__conformers)

    @property
    def conformers(self):
        """
        Returns the conformers list. Includes a check whether there are any conformers left.

        :return: List of conformers.
        """
        # TODO - no checks for now
        return self.__conformers

    @conformers.setter
    def conformers(self, confs: list[MoleculeData]):
        """
        Set the conformers list.

        :param confs: List of MoleculeData instances.
        :type confs: list[MoleculeData]
        :raises ValueError: If not all objects are MoleculeData instances.
        """
        if not all(isinstance(conf, MoleculeData) for conf in confs):
            raise ValueError("All objects need to be MoleculeData instances.")
        self.__conformers = confs

    @property
    def rem(self) -> list[MoleculeData]:
        """
        Returns the list of removed conformers.

        :return: List of removed conformers.
        :rtype: list[MoleculeData]
        """
        return self.__rem

    def read_output(self, outpath: str | Path) -> None:
        """
        Read json output file of a previous execution. Will try to load data into current conformer ensemble, matching
        based on names. If a conformer name does not exist in the current ensemble it will be ignored. If a conformer
        does not exist in the output data RuntimeError will be raised.

        :param outpath: Path to the output file.
        :type outpath: str | Path
        :return: None
        :rtype: None
        """

        data = json.loads(Path(outpath).read_text())

        # Check if all conformers from the current ensemble are also found in the output data
        if not all(conf.name in data for conf in self.conformers):
            raise RuntimeError(
                "Not all conformers from the current ensemble are found in the output data."
            )

        for conf in self.conformers:
            contributions = Contributions(
                data[conf.name]["energy"],
                data[conf.name]["gsolv"],
                data[conf.name]["grrho"],
            )
            conf.update(contributions)

        logger.info(f"Reloaded results from {outpath}.")

    def read_input(
        self,
        input_path: str,
        charge: int = 0,
        unpaired: int = 0,
        nconf: int | None = None,
        append: bool = False,
    ) -> None:
        """
        Read ensemble input file. Should be a file in xyz-file format with all the conformers in consecutive order.

        :param input_path: Path to the ensemble input file.
        :type input_path: str
        :param charge: Sets the charge of all molecules to this value. Defaults to 0.
        :type charge: int
        :param unpaired: Sets the unpaired electrons of all molecules to this value. Defaults to 0.
        :type unpaired: int
        :param nconf: Number of conformers to consider. Defaults to None, so all conformers are read.
        :type nconf: int | None
        :param append: If True, the conformers will be appended to the existing ensemble. Defaults to False.
        :type append: bool
        :return: None
        :rtype: None
        """
        # If $coord in file => tm format, needs to be converted to xyz
        with open(input_path) as inp:
            lines = inp.readlines()
            if any("$coord" in line for line in lines):
                _, _, input_path = t2x(
                    input_path, writexyz=True, outfile="converted.xyz"
                )

        confs = self.__setup_conformers(input_path)
        if len(confs) == 0:
            logger.warning("Input file is empty!")

        if nconf is None:
            nconf = len(confs)

        if append:
            self.conformers.extend(confs[:nconf])
        else:
            self.conformers = confs[:nconf]

        try:
            self.conformers.sort(key=lambda x: x.xtb_energy)
        except TypeError:
            # Only sort if all conformers have a defined precalculated energy
            pass

        for conf in self.conformers:
            conf.charge = charge
            conf.unpaired = unpaired

        # Print information about read ensemble
        printf(
            f"Read {len(self.conformers)} conformers.\n\n",
        )

    def get_populations(self, temperature: float) -> dict[str, float]:
        """
        Calculate populations for boltzmann distribution of ensemble at given
        temperature given values for free enthalpy.

        :param temperature: Temperature for Boltzmann distribution.
        :type temperature: float
        :return: Dictionary of conformer names to populations.
        :rtype: dict[str, float]
        """
        # find lowest gtot value
        minfree: float = min(conf.gtot for conf in self.conformers)

        # calculate boltzmann factors
        bmfactors = {
            conf.name: conf.degen
            * exp(-(conf.gtot - minfree) * AU2J / (KB * temperature))
            for conf in self.conformers
        }

        # calculate partition function from boltzmann factors
        bsum: float = sum(bmfactors.values())

        # Set Boltzmann populations
        populations: dict[str, float] = {}
        for conf in self.conformers:
            populations[conf.name] = bmfactors[conf.name] / bsum

        return populations

    def update_contributions(self, contributions_dict: dict[str, Contributions]):
        """
        Update contributions for all conformers in the ensemble. Convenience wrapper to call update on all conformers with correct mapping.

        :param contributions_dict: Dictionary mapping conformer names to Contributions.
        :type contributions_dict: dict[str, Contributions]
        :raises AssertionError: If contributions dict is incomplete.
        """
        assert all(
            conf.name in contributions_dict for conf in self.conformers
        ), "Attempted to pass incomplete contributions dict to update_contributions."
        for conf in self.conformers:
            conf.update(contributions_dict[conf.name])

    def __setup_conformers(self, input_path: str) -> list[MoleculeData]:
        """
        Open ensemble input, split into conformers, create MoleculeData objects out of coord input, read out energy from xyz file if possible.
        In principle this can also read xyz-files with molecules of different sizes.

        :param input_path: Path to the ensemble input file.
        :type input_path: str
        :return: A list of MoleculeData objects.
        :rtype: list[MoleculeData]
        """
        # open ensemble input
        with open(input_path) as file:
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

    def remove_conformers(self, cond: Callable[[MoleculeData], bool]) -> None:
        """
        Remove conformers from further consideration if 'cond' evaluates to True.
        The removed conformers will be stored in self.rem.

        :param cond: Condition to check for the conf objects.
        :type cond: Callable[[MoleculeData], bool]
        :return: None
        :rtype: None
        """
        filtered = list(filter(cond, self.conformers))
        for conf in filtered:
            # pop item from conformers and insert this item at index 0 in rem
            self.rem.insert(0, self.conformers.pop(self.conformers.index(conf)))

            # Log removed conformers
            logger.debug(f"Removed {conf.name}.")

    def dump_xyz(self, file: Path):
        """
        Dump the current ensemble in xyz-format.

        :param file: Path to write the xyz file.
        :type file: Path
        :return: None
        :rtype: None
        """
        text = "".join(sum([conf.geom.toxyz() for conf in self], []))
        file.write_text(text)

    def dump_rem_xyz(self, file: Path):
        """
        Dump the conformers removed via 'remove_conformers' in xyz-format.

        :param file: Path to write the xyz file.
        :type file: Path
        :return: None
        :rtype: None
        """
        text = "".join(sum([conf.geom.toxyz() for conf in self.rem], []))
        file.write_text(text)

    def dump_json(self, file: Path):
        """
        Dump the ensemble with most recent rankings and values in json-format.

        :param file: Path to write the json file.
        :type file: Path
        :return: None
        :rtype: None
        """
        dump = {
            conf.name: {
                "energy": conf.energy,
                "gsolv": conf.gsolv,
                "grrho": conf.grrho,
            }
            for conf in self.conformers
        }
        file.write_text(json.dumps(dump, indent=4))
