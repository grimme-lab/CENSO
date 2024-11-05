"""
Contains boilerplate class for calculating ensemble properties.
"""

from ..logging import setup_logger
from ..part import CensoPart
from ..utilities import timeit, SolventHelper
from ..datastructure import MoleculeData
from ..ensembleopt import EnsembleOptimizer

logger = setup_logger(__name__)


class PropertyCalculator(CensoPart):
    """
    Boilerplate class for all property calculations.
    """

    _grid = ""

    @timeit
    @CensoPart._create_dir
    def __call__(self, using_part: CensoPart = None) -> None:
        """
        Boilerplate run logic for any ensemble property calculation. The 'property' method should be implemented for every
        class respectively.

        Running a property calculation requires some kind of ensemble energetic ranking beforehand.

        It is possible to pass a specific part output to determine the Boltzmann populations
        """
        # print instructions
        self._print_info()

        # Set energy values to use later
        self._set_energy(using_part=using_part)
        for conf in self._ensemble.conformers:
            self.data["results"][conf.name]["gtot"] = self._gtot(conf)

        # Calculate Boltzmann populations
        self._update_results(self._calc_boltzmannweights())

        # Perform the property calculations
        self._property()

        # Write out the results
        self._write_results()

        # DONE

    def _output(self) -> None:
        """
        Implements printouts and writes for any output data.
        Necessary to implement for each part.
        """
        # Write out results
        self._write_results()

    def _property(self):
        raise NotImplementedError

    def _write_results(self):
        raise NotImplementedError

    def _gtot(self, conf: MoleculeData) -> float:
        return (
            self.data["results"][conf.name]["energy"]
            + self.data["results"][conf.name]["gsolv"]
            + self.data["results"][conf.name]["grrho"]
        )

    def _setup_prepinfo_rrho(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self.name
        prepinfo["charge"] = self._ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self._ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        prepinfo["xtb_rrho"] = {
            "gfnv": self.get_settings()["gfnv"],
        }
        # Only lookup solvent if solvation should be used
        if not self.get_general_settings()["gas-phase"]:
            prepinfo["xtb_rrho"]["solvent_key_xtb"] = SolventHelper.get_solvent(
                self.get_general_settings()["sm_rrho"],
                self.get_general_settings()["solvent"],
            )

        return prepinfo

    def _set_energy(self, using_part: CensoPart = None):
        """
        Looks through results to set energy values.
        Order of preference:
            refinement -> optimization -> screening -> prescreening

        If None of these are found, raise RuntimeError.
        """
        if using_part is None:
            # Determine the smallest usable optimization results
            # First filter ensemble optimizations
            opts = filter(
                lambda part: issubclass(type(part), EnsembleOptimizer),
                self._ensemble.results,
            )
            opts = sorted(opts, key=lambda part: part.data["nconf_out"])

            # Get the results with the smallest outputs
            opts_iter = iter(opts)
            res = next(opts_iter, None)
            if res is None:
                raise RuntimeError(
                    "Calculating an ensemble property requires some kind of energetic ensemble ranking performed beforehand."
                )

            smallest_results = []
            while res.data["nconf_out"] == opts[0].data["nconf_out"]:
                smallest_results.append(res)
                try:
                    res = next(opts_iter)
                except StopIteration:
                    break

            # Get the highest (assumed) quality part from those
            if len(smallest_results) == 1:
                using_part = smallest_results[0]
            else:
                # This will put the highest quality part at the top (highest part number)
                smallest_results.sort(
                    lambda part: self._part_nos[part.name], reverse=True
                )
                using_part = smallest_results[0]

        # Get the index of this results from the ensemble results
        assert using_part is not None
        using_part = self._ensemble.results.index(using_part)

        energy_values = {
            "prescreening": lambda conf: {
                "energy": self._ensemble.results[using_part].data["results"][conf.name][
                    "sp"
                ]["energy"],
                "gsolv": (
                    self._ensemble.results[using_part].data["results"][conf.name][
                        "xtb_gsolv"
                    ]["gsolv"]
                    if "xtb_gsolv"
                    in self._ensemble.results[using_part].data["results"][conf.name]
                    else 0.0
                ),
                "grrho": 0.0,
            },
            "screening": lambda conf: {
                "energy": (
                    self._ensemble.results[using_part].data["results"][conf.name][
                        "gsolv"
                    ]["energy_gas"]
                    if "gsolv"
                    in self._ensemble.results[using_part].data["results"][conf.name]
                    else self._ensemble.results[using_part].data["results"][conf.name][
                        "sp"
                    ]["energy"]
                ),
                "gsolv": (
                    self._ensemble.results[using_part].data["results"][conf.name][
                        "gsolv"
                    ]["gsolv"]
                    if "gsolv"
                    in self._ensemble.results[using_part].data["results"][conf.name]
                    else 0.0
                ),
                "grrho": self._ensemble.results[using_part]
                .data["results"][conf.name]
                .get("xtb_rrho", {"energy": 0.0})["energy"],
            },
            "optimization": lambda conf: {
                "energy": self._ensemble.results[using_part].data["results"][conf.name][
                    "xtb_opt"
                ]["energy"],
                "gsolv": 0.0,
                "grrho": self._ensemble.results[using_part]
                .data["results"][conf.name]
                .get("xtb_rrho", {"energy": 0.0})["energy"],
            },
            "refinement": lambda conf: {
                "energy": (
                    self._ensemble.results[using_part].data["results"][conf.name][
                        "gsolv"
                    ]["energy_gas"]
                    if "gsolv"
                    in self._ensemble.results[using_part].data["results"][conf.name]
                    else self._ensemble.results[using_part].data["results"][conf.name][
                        "sp"
                    ]["energy"]
                ),
                "gsolv": (
                    self._ensemble.results[using_part].data["results"][conf.name][
                        "gsolv"
                    ]["gsolv"]
                    if "gsolv"
                    in self._ensemble.results[using_part].data["results"][conf.name]
                    else 0.0
                ),
                "grrho": self._ensemble.results[using_part]
                .data["results"][conf.name]
                .get("xtb_rrho", {"energy": 0.0})["energy"],
            },
        }

        for conf in self._ensemble.conformers:
            self.data["results"].setdefault(
                conf.name, energy_values[self._ensemble.results[using_part].name](conf)
            )
