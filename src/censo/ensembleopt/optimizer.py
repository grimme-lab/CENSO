from ..ensembledata import EnsembleData
from ..logging import setup_logger
from ..params import AU2KCAL, DIGILEN, PLENGTH
from ..part import CensoPart
from ..utilities import DfaHelper, SolventHelper, format_data, h1, print, timeit

logger = setup_logger(__name__)


class EnsembleOptimizer(CensoPart):
    """
    Boilerplate class for all ensemble optimization steps.
    """

    _grid = ""

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self, ncores: int, cut: bool = True) -> None:
        """
        Boilerplate run logic for any ensemble optimization step. The 'optimize' method should be implemented for every
        class respectively.
        """
        # print instructions
        self.print_info()

        # Print information about ensemble before optimization
        self.print_update()

        # Perform the actual optimization logic
        self.optimize(ncores, cut=cut)

        # Print comparison with previous parts
        self.print_comparison()

        # Print information about ensemble after optimization
        self.print_update()

        # dump ensemble
        self.ensemble.dump_ensemble(f"{self._part_no}_{self._name.upper()}")

        # DONE

    def optimize(self, ncores: int, cut: bool = True):
        raise NotImplementedError

    def setup_prepinfo(self, jobtype: list[str]) -> dict[str, dict]:
        """
        Sets up lookup information to be used by the processor in parallel execution. Returns a dictionary
        containing all information for all jobtypes provided.

        Args:
            jobtype (list[str]): list of jobtypes to be run.

        Returns:
            dict[str, dict]: dictionary containing all information for all jobtypes provided.
        """
        prepinfo = {jt: {} for jt in jobtype}

        prepinfo["partname"] = self._name
        prepinfo["charge"] = self.ensemble.runinfo.get("charge")
        prepinfo["unpaired"] = self.ensemble.runinfo.get("unpaired")
        prepinfo["general"] = self.get_general_settings()

        if "sp" in jobtype or "gsolv" in jobtype:
            prepinfo["sp"] = {
                "func_name": DfaHelper.get_name(
                    self.get_settings()["func"], self.get_settings()["prog"]
                ),
                "func_type": DfaHelper.get_type(self.get_settings()["func"]),
                "disp": DfaHelper.get_disp(self.get_settings()["func"]),
                "basis": self.get_settings()["basis"],
                "grid": self._grid,
                "template": self.get_settings()["template"],
                "gcp": True,
            }

            # Add the solvent key if a solvent model exists in the part settings
            # NOTE: 'sm' in key catches also cases like NMR (sm_s and sm_j)
            # Only look up solvent if solvation is used
            if (
                "sm" in self.get_settings()
                and not self.get_general_settings()["gas-phase"]
            ):
                prepinfo["sp"]["sm"] = self.get_settings()["sm"]
                prepinfo["sp"]["solvent_key_prog"] = SolventHelper.get_solvent(
                    self.get_settings()["sm"], self.get_general_settings()["solvent"]
                )
                assert prepinfo["sp"]["solvent_key_prog"] is not None

            if (
                self.get_settings()["prog"] == "tm"
                and prepinfo["sp"]["disp"] == "d4"
                and prepinfo["sp"]["gcp"]
            ):
                logger.warning(
                    "Due to a bug in TURBOMOLE it is currently not possible to use GCP "
                    "together with the D4 correction. GCP will be disabled."
                )
                prepinfo["sp"]["gcp"] = False

        # TODO - this doesn't look very nice
        if "xtb_gsolv" in jobtype:
            prepinfo["xtb_sp"] = {
                "gfnv": self.get_settings()["gfnv"],
                "solvent_key_xtb": SolventHelper.get_solvent(
                    self.get_general_settings()["sm_rrho"],
                    self.get_general_settings()["solvent"],
                ),
            }
            # gsolv implies that solvation should be used, so no check here
            assert prepinfo["xtb_sp"]["solvent_key_xtb"] is not None

        if "xtb_rrho" in jobtype:
            prepinfo["xtb_rrho"] = {
                "gfnv": self.get_settings()["gfnv"],
            }
            # Only look up solvent if solvation is used
            if not self.get_general_settings()["gas-phase"]:
                prepinfo["xtb_rrho"]["solvent_key_xtb"] = SolventHelper.get_solvent(
                    self.get_general_settings()["sm_rrho"],
                    self.get_general_settings()["solvent"],
                )
                assert prepinfo["xtb_rrho"]["solvent_key_xtb"] is not None

        for jt in ["xtb_opt", "opt"]:
            if jt in jobtype:
                prepinfo[jt] = {
                    "func_name": DfaHelper.get_name(
                        self.get_settings()["func"], self.get_settings()["prog"]
                    ),
                    "func_type": DfaHelper.get_type(self.get_settings()["func"]),
                    "disp": DfaHelper.get_disp(self.get_settings()["func"]),
                    "basis": self.get_settings()["basis"],
                    "grid": self._grid,
                    "template": self.get_settings()["template"],
                    "gcp": True,
                    "optcycles": self.get_settings()["optcycles"],
                    "hlow": self.get_settings()["hlow"],
                    "optlevel": self.get_settings()["optlevel"],
                    "macrocycles": self.get_settings()["macrocycles"],
                    "constraints": self.constraints,
                    # this is set to a path if constraints should be used, otherwise None
                }

                # Only look up solvent if solvation is used
                if not self.get_general_settings()["gas-phase"]:
                    prepinfo[jt]["sm"] = self.get_settings()["sm"]
                    prepinfo[jt]["solvent_key_prog"] = SolventHelper.get_solvent(
                        self.get_settings()["sm"],
                        self.get_general_settings()["solvent"],
                    )
                    assert prepinfo[jt]["solvent_key_prog"] is not None

                if (
                    self.get_settings()["prog"] == "tm"
                    and prepinfo["sp"]["disp"] == "d4"
                    and prepinfo["sp"]["gcp"]
                ):
                    logger.warning(
                        "Due to a bug in TURBOMOLE it is currently not possible to use GCP "
                        "together with the D4 correction. GCP will be disabled."
                    )
                    prepinfo["sp"]["gcp"] = False

                break

        return prepinfo

    def print_update(self) -> None:
        print("\n")
        print(
            "Number of conformers:".ljust(DIGILEN // 2, " ")
            + f"{len(self.ensemble.conformers)}"
        )
        print(
            "Highest ranked conformer:".ljust(DIGILEN // 2, " ")
            + f"{self.ensemble.conformers[0].name}"
        )
        print("\n")

    def print_comparison(self) -> None:
        print(h1(f"{self._name.upper()} RANKING COMPARISON"))

        headers = ["CONF#"]

        parts = list(self.ensemble.conformers[0].results.keys())

        headers.extend([f"ΔGtot {part}" for part in parts])

        # column units
        units = [
            "",
        ]

        units.extend(["[kcal/mol]" for _ in range(len(parts))])

        # variables for printmap
        gtotmin = {part: 0.0 for part in parts}
        for part in parts:
            gtotmin[part] = min(
                conf.results[part]["gtot"] for conf in self.ensemble.conformers
            )

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
        }
        for header, part in zip(headers[1:], parts):
            # Same lambda bullshittery as in parallel.py/dqp, python needs the lambda kwargs or it will
            # use the same values for every lambda call
            printmap[header] = (
                lambda conf, partl=part, headerl=header: f"{(conf.results[partl]['gtot'] - gtotmin[partl]) * AU2KCAL:.2f}"
            )

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # Print everything
        for line in lines:
            print(line, flush=True, end="")

        print("".ljust(int(PLENGTH), "-") + "\n")
