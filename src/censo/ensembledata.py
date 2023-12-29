"""
stores ensembledata and conformers
functionality for program setup
"""

from argparse import Namespace
import os
from math import exp
from collections.abc import Callable

from .params import AU2J, KB, AU2KCAL
from .datastructure import MoleculeData
from .utilities import (
    check_for_float,
    do_md5,
    t2x,
    setup_logger,
)

logger = setup_logger(__name__)


class EnsembleData:
    """ """

    def __init__(self, workdir: str, args: Namespace = None):
        """
        Setup a EnsembleData object using the args from the command line
        workdir is the directory where the CENSO run should be executed in
        """

        # current working directory
        self.workdir: str = workdir

        # if args are given set accordingly, otherwise assume CENSO is used without commandline
        self.args: Namespace = args

        # contains run-specific info that may change during runtime
        # initialized in EnsembleData.read_input
        self.runinfo = {
            "nconf": None,
            "nat": None,  # TODO - maybe remove?
            "maxconf": None,  # TODO - probably remove
            "md5": None,
            "consider_unconverged": None,  # TODO - probably remove
            "charge": None,
            "unpaired": None,
        }

        # stores the conformers with all info
        # NOTE: this is deliberately chosen to be a list since lists are ordered
        self.__conformers: list[MoleculeData] = []

        # stores the conformers which were sorted out
        self.rem: list[MoleculeData] = []

        # absolute path to ensemble input file
        self.ensemble_path: str

    @property
    def conformers(self):
        """
        Returns the conformers list. Includes a check wether there are any conformers left.
        """
        assert len(self.conformers) > 0
        return self.__conformers

    @conformers.setter
    def conformers(self, confs):
        assert all(isinstance(conf, MoleculeData) for conf in confs)
        self.__conformers = confs

    def read_input(
        self,
        ensemble_path: str,
        charge: int = None,
        unpaired: int = None,
        nconf: int = None,
    ) -> None:
        """
        read ensemble input file (e.g. crest_conformers.xyz)
        """

        self.ensemble_path = ensemble_path

        # store md5 hash for quick comparison of inputs later
        self.runinfo["md5"] = do_md5(self.ensemble_path)

        # if $coord in file => tm format, needs to be converted to xyz
        with open(self.ensemble_path, "r") as inp:
            lines = inp.readlines()
            if any(["$coord" in line for line in lines]):
                _, self.runinfo["nat"], self.ensemble_path = t2x(
                    self.ensemble_path, writexyz=True, outfile="converted.xyz"
                )
            else:
                self.runinfo["nat"] = int(lines[0].split()[0])

        # set charge and unpaired via funtion args or cml args
        if self.args is not None:
            self.runinfo["charge"] = charge or self.args.charge
            self.runinfo["unpaired"] = unpaired or self.args.unpaired
        else:
            self.runinfo["charge"] = charge
            self.runinfo["unpaired"] = unpaired

        if self.runinfo["charge"] is None or self.runinfo["unpaired"] is None:
            raise RuntimeError(
                "Charge or number of unpaired electrons not defined.")

        self.setup_conformers(nconf)

    def setup_conformers(self, maxconf: int) -> None:
        """
        open ensemble input
        split into conformers
        create MoleculeData objects out of coord input
        read out energy from xyz file if possible
        """
        # open ensemble input
        with open(self.ensemble_path, "r") as file:
            lines = file.readlines()
            nat = self.runinfo["nat"]

            # check for correct line count in input
            # assuming consecutive xyz-file format (every conf therefore needs nat + 2 lines)
            if len(lines) % (nat + 2) != 0:
                # check if the last lines contain any alphanumeric characters
                divisor = len(lines) // (nat + 2)
                if any([line.isalnum() for line in lines[divisor:]]):
                    raise RuntimeError("Could not read ensemble input file.")
                else:
                    # if not you can trim the last lines
                    lines = lines[divisor:]

            if self.args is not None and self.args.nconf is not None:
                nconf = int(len(lines) / (nat + 2))
                if self.args.nconf > nconf:
                    global logger
                    logger.warning(
                        f"Provided nconf is larger than number of conformers in input file. Setting to "
                        f"the max. amount automatically."
                    )
                else:
                    nconf = self.args.nconf
            else:
                nconf = maxconf or int(len(lines) / (nat + 2))

            self.runinfo["nconf"] = nconf

            # get precalculated energies if possible
            for i in range(nconf):
                self.conformers.append(
                    MoleculeData(
                        f"CONF{i + 1}", lines[2 + i *
                                              (nat + 2):(i + 1) * (nat + 2)]
                    )
                )

                # precalculated energy set to 0.0 if it cannot be found
                self.conformers[i].xtb_energy = (
                    check_for_float(lines[i * (nat + 2) + 1]) or 0.0
                )

            # also works if xtb_energy is None for some reason (None is put first)
            self.conformers.sort(key=lambda x: x.xtb_energy)

    def update_conformers(
        self,
        target: Callable[[MoleculeData], float],
        threshold: float,
        additional_filter: Callable[[MoleculeData], bool] = None,
        boltzmann: bool = False,
    ) -> list[str]:
        """
        Update the conformers based on a target function, a threshold, and an additional filter.
        Logic for the additional filter is "and".
        Returns a list of the names of the removed conformers.

        Args:
            target (Callable[[MoleculeData], float]): A function that takes a MoleculeData object as input and returns a
                                                      float.
            threshold (float): The threshold value.
            additional_filter (Callable[[MoleculeData], bool], optional): An optional function that takes a MoleculeData
                                                                          object as input and returns a boolean.
                                                                          Defaults to None.

        Returns:
            None
        """
        if not boltzmann:
            # pick the free enthalpy of the lowest conformer
            limit = min([target(conf) for conf in self.conformers])

            # filter out all conformers above threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf
                for conf in filter(
                    (additional_filter or True)
                    and (lambda x: target(x) - limit > threshold),
                    self.conformers,
                )
            ]
        elif boltzmann and 0.0 <= threshold <= 1.0:
            # Sort and iterate through the conformers by target (should be the Boltzmann population) in reverse order
            # Therefore, the conformer with the highest population is first
            s = 0.0
            filtered = []
            for conf in sorted(self.conformers, key=target, reverse=True):
                if s > threshold:
                    # The conformer is above the population threshold and should be removed
                    filtered.append(conf)
                else:
                    # The population of the conformer is appended
                    s += target(conf)
        else:
            raise RuntimeError(
                "Invalid filter settings for updating conformer list.")

        # move the sorted out conformers to rem list
        for conf in filtered:
            # pop item from conformers and insert this item at index 0 in rem
            self.rem.insert(0, self.conformers.pop(
                self.conformers.index(conf)))

            # Log removed conformers
            logger.debug(f"Removed {conf.name}.")

        # TODO - signal CENSO to stop if there are no conformers left

        return [conf.name for conf in filtered]

    def dump_ensemble(self, part: str) -> None:
        """
        dump the conformers to a file
        """
        with open(
            os.path.join(f"{self.workdir}", f"censo_ensemble_{part}.xyz"), "w"
        ) as file:
            for conf in self.conformers:
                file.writelines(conf.geom.toxyz())

    def calc_boltzmannweights(self, temp: float, part: str) -> None:
        """
        Calculate weights for boltzmann distribution of ensemble at given temperature and part name to search results,
        given values for free enthalpy
        """
        # find lowest gtot value
        if all(["gtot" in conf.results[part].keys() for conf in self.conformers]):
            minfree: float = min(
                [conf.results[part]["gtot"] for conf in self.conformers]
            )

            # calculate boltzmann factors
            bmfactors = {
                id(conf): conf.degen
                * exp(-(conf.results[part]["gtot"] - minfree) * AU2J / (KB * temp))
                for conf in self.conformers
            }
        else:
            # NOTE: if anything went wrong in the single-point calculation ("success": False),
            # this should be handled before coming to this step
            # since then the energy might be 'None'
            gtot_replacement = False
            for jt in ["xtb_opt", "sp"]:
                if all(jt in conf.results[part].keys() for conf in self.conformers):
                    minfree: float = min(
                        [conf.results[part][jt]["energy"]
                            for conf in self.conformers]
                    )

                    # calculate boltzmann factors
                    bmfactors = {
                        id(conf): conf.degen
                        * exp(
                            -(conf.results[part][jt]["energy"] - minfree)
                            * AU2J
                            / (KB * temp)
                        )
                        for conf in self.conformers
                    }
                    gtot_replacement = True
                    break

            if not gtot_replacement:
                raise RuntimeError(
                    f"Could not determine Boltzmann factors for {part}.")

        # calculate partition function from boltzmann factors
        bsum: float = sum(bmfactors.values())

        # Store Boltzmann populations in results and also in a special list for convenience
        for conf in self.conformers:
            conf.results[part]["bmw"] = bmfactors[id(conf)] / bsum
            conf.bmws.append(bmfactors[id(conf)] / bsum)
