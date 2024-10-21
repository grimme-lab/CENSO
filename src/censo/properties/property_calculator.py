"""
Contains boilerplate class for calculating ensemble properties.
"""

from ..ensembledata import EnsembleData
from ..logging import setup_logger
from ..part import CensoPart
from ..utilities import timeit, SolventHelper
from ..datastructure import MoleculeData
from ..parallel import execute

logger = setup_logger(__name__)


class PropertyCalculator(CensoPart):
    """
    Boilerplate class for all property calculations.
    """

    _grid = ""

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self, ncores: int) -> None:
        """
        Boilerplate run logic for any ensemble property calculation. The 'property' method should be implemented for every
        class respectively.

        Running a property calculation requires some kind of ensemble energetic ranking beforehand.
        TODO - maybe add option to turn that behaviour off?
        """
        # print instructions
        self.print_info()

        # Set energy values to use later
        self._set_energy()
        for conf in self.ensemble.conformers:
            conf.results[self._name]["gtot"] = self.gtot(conf)

        # Calculate Boltzmann populations
        self.ensemble.calc_boltzmannweights(
            self.get_general_settings()["temperature"], self._name
        )

        # Perform the property calculations
        self.property(ncores)

        # DONE

    def property(self, ncores: int):
        raise NotImplementedError

    def gtot(self, conformer: MoleculeData) -> float:
        resdict = conformer.results[self._name]
        return resdict["energy"] + resdict["gsolv"] + resdict["grrho"]

    def setup_prepinfo_rrho(self) -> dict[str, dict]:
        prepinfo = {}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
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

    def _set_energy(self):
        """
        Looks through results to set energy values.
        Order of preference:
            refinement -> optimization -> screening -> prescreening

        If None of these are found, raise RuntimeError.
        """
        using_part = next(
            (
                partname
                for partname in [
                    "refinement",
                    "optimization",
                    "screening",
                    "prescreening",
                ]
                if all(partname in conf.results for conf in self.ensemble.conformers)
            ),
            None,
        )

        if using_part is None:
            raise RuntimeError(
                "Calculating an ensemble property requires some kind of energetic ensemble ranking performed beforehand."
            )

        energy_values = {
            "optimization": lambda conf: {
                "energy": conf.results[using_part]["xtb_opt"]["energy"],
                "gsolv": 0.0,
                "grrho": conf.results[using_part].get("xtb_rrho", {"energy": 0.0})[
                    "energy"
                ],
            },
            "screening": lambda conf: {
                "energy": (
                    conf.results[using_part]["gsolv"]["energy_gas"]
                    if "gsolv" in conf.results[using_part]
                    else conf.results[using_part]["sp"]["energy"]
                ),
                "gsolv": (
                    conf.results[using_part]["gsolv"]["gsolv"]
                    if "gsolv" in conf.results[using_part]
                    else 0.0
                ),
                "grrho": conf.results[using_part].get("xtb_rrho", {"energy": 0.0})[
                    "energy"
                ],
            },
            "refinement": lambda conf: {
                "energy": (
                    conf.results[using_part]["gsolv"]["energy_gas"]
                    if "gsolv" in conf.results[using_part]
                    else conf.results[using_part]["sp"]["energy"]
                ),
                "gsolv": (
                    conf.results[using_part]["gsolv"]["gsolv"]
                    if "gsolv" in conf.results[using_part]
                    else 0.0
                ),
                "grrho": conf.results[using_part].get("xtb_rrho", {"energy": 0.0})[
                    "energy"
                ],
            },
            "prescreening": lambda conf: {
                "energy": conf.results[using_part]["sp"]["energy"],
                "gsolv": (
                    conf.results[using_part]["xtb_gsolv"]["gsolv"]
                    if "xtb_gsolv" in conf.results[using_part]
                    else 0.0
                ),
                "grrho": 0.0,
            },
        }

        for conf in self.ensemble.conformers:
            conf.results.setdefault(self._name, energy_values[using_part](conf))
