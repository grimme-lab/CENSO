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
            assert prepinfo["xtb_rrho"]["solvent_key_xtb"] is not None

        return prepinfo

    def _set_energy(self):
        """
        Looks through results to set energy values.
        Order of preference:
            refinement -> optimization -> screening -> prescreening

        If None of these are found, raise RuntimeError.
        """
        using_part = None
        for partname in ["prescreening", "screening", "optimization", "refinement"]:
            # This way, the most high-level partname should get stuck in using_part
            if all(
                partname in conf.results.keys() for conf in self.ensemble.conformers
            ):
                using_part = partname

        if using_part is None:
            raise RuntimeError(
                "Calculating an ensemble property requires some kind of energetic ensemble ranking performed beforehand."
            )

        if using_part == "optimization":
            for conf in self.ensemble.conformers:
                conf.results[self._name]["energy"] = conf.results[using_part][
                    "xtb_opt"
                ]["energy"]
                conf.results[self._name]["gsolv"] = 0.0
                conf.results[self._name]["grrho"] = conf.results[using_part].get(
                    "xtb_rrho", {"energy": 0.0}
                )["energy"]
        elif using_part in ["screening", "refinement"]:
            if all(
                "gsolv" in conf.results[using_part].keys()
                for conf in self.ensemble.conformers
            ):
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["energy"] = conf.results[using_part][
                        "gsolv"
                    ]["energy_gas"]
                    conf.results[self._name]["gsolv"] = conf.results[using_part][
                        "gsolv"
                    ]["gsolv"]
                    conf.results[self._name]["grrho"] = conf.results[using_part].get(
                        "xtb_rrho", {"energy": 0.0}
                    )["energy"]
            else:
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["energy"] = conf.results[using_part]["sp"][
                        "energy"
                    ]
                    conf.results[self._name]["gsolv"] = 0.0
                    conf.results[self._name]["grrho"] = conf.results[using_part].get(
                        "xtb_rrho", {"energy": 0.0}
                    )["energy"]
        elif using_part == "prescreening":
            if all(
                "xtb_gsolv" in conf.results[using_part].keys()
                for conf in self.ensemble.conformers
            ):
                for conf in self.ensemble.conformers:
                    conf.results[self._name]["energy"] = conf.results[using_part]["sp"][
                        "energy"
                    ]
                    conf.results[self._name]["gsolv"] = conf.results[using_part][
                        "xtb_gsolv"
                    ]["gsolv"]
                    conf.results[self._name]["grrho"] = 0.0
