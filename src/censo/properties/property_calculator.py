"""
Contains boilerplate class for calculating ensemble properties.
"""

from ..ensembledata import EnsembleData
from ..logging import setup_logger
from ..part import CensoPart
from ..utilities import timeit, SolventHelper
from ..datastructure import MoleculeData
from ..ensembleopt.optimizer import EnsembleOptimizer
from ..ensembleopt.prescreening import Prescreening
from ..ensembleopt.screening import Screening
from ..ensembleopt.optimization import Optimization
from ..ensembleopt.refinement import Refinement

logger = setup_logger(__name__)


class PropertyCalculator(CensoPart):
    """
    Boilerplate class for all property calculations.
    """

    _grid = ""

    @timeit
    @CensoPart._create_dir
    def __call__(self, ensemble: EnsembleData) -> None:
        """
        Boilerplate run logic for any ensemble property calculation. The 'property' method should be implemented for every
        class respectively.

        Running a property calculation requires some kind of ensemble energetic ranking beforehand.
        """
        # print instructions
        self._print_info()

        # Set energy values to use later
        self._set_energy()
        for conf in self.ensemble.conformers:
            self.results[conf.name]["gtot"] = self._gtot(conf)

        # Calculate Boltzmann populations
        self.results.update(self._calc_boltzmannweights())

        # Perform the property calculations
        self._property()

        # DONE

    def _property(self):
        raise NotImplementedError

    def _gtot(self, conformer: MoleculeData) -> float:
        return (
            self.results[conformer.name]["energy"]
            + self.results[conformer.name]["gsolv"]
            + self.results[conformer.name]["grrho"]
        )

    def _setup_prepinfo_rrho(self) -> dict[str, dict]:
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
        Uses the most recent ensemble optimization part as reference.

        If None are found, raise RuntimeError.
        """
        using_part = next(
            (
                p
                for p in self.ensemble.results[::-1]
                if issubclass(type(p), EnsembleOptimizer)
            ),
            (None, None),
        )

        if using_part is None:
            raise RuntimeError(
                "Calculating an ensemble property requires some kind of energetic ensemble ranking performed beforehand."
            )

        energy_values = {
            Optimization: lambda conf: {
                "energy": using_part.results[conf.name]["xtb_opt"]["energy"],
                "gsolv": 0.0,
                "grrho": using_part.results[conf.name].get("xtb_rrho", {"energy": 0.0})[
                    "energy"
                ],
            },
            Screening: lambda conf: {
                "energy": (
                    using_part.results[conf.name]["gsolv"]["energy_gas"]
                    if "gsolv" in using_part.results[conf.name]
                    else using_part.results[conf.name]["sp"]["energy"]
                ),
                "gsolv": (
                    using_part.results[conf.name]["gsolv"]["gsolv"]
                    if "gsolv" in using_part.results[conf.name]
                    else 0.0
                ),
                "grrho": using_part.results[conf.name].get("xtb_rrho", {"energy": 0.0})[
                    "energy"
                ],
            },
            Refinement: lambda conf: {
                "energy": (
                    using_part.results[conf.name]["gsolv"]["energy_gas"]
                    if "gsolv" in using_part.results[conf.name]
                    else using_part.results[conf.name]["sp"]["energy"]
                ),
                "gsolv": (
                    using_part.results[conf.name]["gsolv"]["gsolv"]
                    if "gsolv" in using_part.results[conf.name]
                    else 0.0
                ),
                "grrho": using_part.results[conf.name].get("xtb_rrho", {"energy": 0.0})[
                    "energy"
                ],
            },
            Prescreening: lambda conf: {
                "energy": using_part.results[conf.name]["sp"]["energy"],
                "gsolv": (
                    using_part.results[conf.name]["xtb_gsolv"]["gsolv"]
                    if "xtb_gsolv" in using_part.results[conf.name]
                    else 0.0
                ),
                "grrho": 0.0,
            },
        }

        for conf in self.ensemble.conformers:
            self.results.setdefault(conf.name, energy_values[type(using_part)](conf))
