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

    @classmethod
    def _validate(cls, tovalidate: dict[str, any]) -> None:
        """
        Validates the type of each setting in the given dict. Also potentially validate if the setting is allowed by
        checking with cls._options.
        This is the part-specific version of the method. It will run the general validation first and then
        check part-specific logic.

        Args:
            tovalidate (dict[str, any]): The dict containing the settings to be validated.

        Returns:
            None

        Raises:
            ValueError: If the setting is not allowed or the value is not within the allowed options.
        """
        # General validation
        super()._validate(tovalidate)

        # Part-specific validation
        # NOTE: tovalidate is always complete
        # Check availability of func for prog
        func = tovalidate["func"]
        if func not in cls._options["func"]["options"][tovalidate["prog"]]:
            raise ValueError(
                f"Functional {func} is not available for {tovalidate['prog']}. "
                "Check spelling w.r.t. CENSO functional naming convention (case insensitive)."
            )

        # Check sm availability for prog
        sm = tovalidate.get("sm", None)
        if (
            sm is not None
            and sm not in cls._options["sm"]["options"][tovalidate["prog"]]
        ):
            raise ValueError(
                f"Solvent model {sm} not available for {tovalidate['prog']}."
            )

        # Check solvent availability for sm
        if (
            sm is not None
            and cls.get_general_settings()["solvent"]
            not in CensoPart._options["solvent"]["options"][sm]
        ):
            raise ValueError(
                f"Solvent {cls.get_general_settings()['solvent']} is not available for {sm}. "
                "Please create an issue on GitHub if you think this is incorrect."
            )

        # dummy/template functionality not implemented yet for TM
        if tovalidate["prog"] == "tm" and (
            func == "dummy" or tovalidate.get("template", False)
        ):
            raise NotImplementedError(
                "Dummy and template functionality is not implemented yet for use with TURBOMOLE."
            )

    @timeit
    @CensoPart._create_dir
    def _run(self, cut: bool = True) -> None:
        """
        Boilerplate run logic for any ensemble optimization step. The 'optimize' method should be implemented for every
        class respectively.
        """
        # print instructions
        self._print_info()

        # Print information about ensemble before optimization
        self._print_update()

        # Perform the actual optimization logic
        self._optimize(cut=cut)

        # Print comparison with previous parts
        self._print_comparison()

        # Print information about ensemble after optimization
        self._print_update()

        # dump ensemble
        self.ensemble.dump(f"{self._part_nos[self._name]}_{self._name.upper()}")

        # DONE

    def _optimize(self, cut: bool = True):
        raise NotImplementedError

    def _cut_conformers(self) -> None:
        """
        Cut down the conformer ensemble based on a given threshold (and threshold type).
        The threshold can be either a kcal/mol value if cutting by Gtot or a population
        threshold when cutting based on populations.
        """
        filtered = []
        threshold = self.get_settings()["threshold"]

        # Refinement cuts based on Boltzmann population
        if self._name == "refinement" and 0.0 <= threshold <= 1.0:
            # Sort and iterate through the conformers by target (should be the Boltzmann population) in reverse order
            # Therefore, the conformer with the highest population is first
            s = 0.0
            for conf in sorted(
                self.ensemble.conformers,
                key=lambda conf: self.results[conf.name]["bmw"],
                reverse=True,
            ):
                if s > threshold:
                    # The conformer is above the population threshold and should be removed
                    filtered.append(conf)
                else:
                    # The population of the conformer is appended
                    s += target(conf)
        # The rest cuts based on gtot
        else:
            # pick the free enthalpy of the lowest conformer
            limit = min(
                self.results[conf.name]["gtot"] for conf in self.ensemble.conformers
            )

            # filter out all conformers above threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = list(
                filter(
                    self.results[conf.name]["gtot"] - limit > threshold,
                    self.ensemble.conformers,
                )
            )

        # move the sorted out conformers to rem list
        for conf in filtered:
            # pop item from conformers and insert this item at index 0 in rem
            self.rem.insert(0, self.conformers.pop(self.conformers.index(conf)))

            # Log removed conformers
            logger.debug(f"Removed {conf.name}.")

        self.ensemble.remove_conformers([conf.name for conf in filtered])

    def _setup_prepinfo(self, jobtype: list[str]) -> dict[str, dict]:
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

            if (
                self.get_settings()["prog"] == "tm"
                and prepinfo["sp"]["disp"] == "d4"
                and prepinfo["sp"]["gcp"]
            ):
                # Basis sets including the following naming patterns should definitely use GCP
                gcp_basis_patterns = ["sv", "dz", "tz", "mini", "6-31g(d)"]
                if any(
                    pattern in prepinfo["sp"]["basis"] for pattern in gcp_basis_patterns
                ):
                    logger.warning(
                        "Due to a bug in TURBOMOLE it is currently not possible to use GCP "
                        "together with the D4 correction. Switching to D3."
                    )
                    prepinfo["sp"]["disp"] = DfaHelper.get_disp(
                        self.get_settings()["func"].replace("d4", "d3")
                    )
                else:
                    logger.warning(
                        "Due to a bug in TURBOMOLE it is currently not possible to use GCP "
                        "together with the D4 correction. Turning off GCP."
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

                if (
                    self.get_settings()["prog"] == "tm"
                    and prepinfo[jt]["disp"] == "d4"
                    and prepinfo[jt]["gcp"]
                ):
                    # Basis sets including the following naming patterns should definitely use GCP
                    gcp_basis_patterns = ["sv", "dz", "tz", "mini", "6-31g(d)"]
                    if any(
                        pattern in prepinfo[jt]["basis"]
                        for pattern in gcp_basis_patterns
                    ):
                        logger.warning(
                            "Due to a bug in TURBOMOLE it is currently not possible to use GCP "
                            "together with the D4 correction. Switching to D3."
                        )
                        prepinfo[jt]["disp"] = DfaHelper.get_disp(
                            self.get_settings()["func"].replace("d4", "d3")
                        )
                    else:
                        logger.warning(
                            "Due to a bug in TURBOMOLE it is currently not possible to use GCP "
                            "together with the D4 correction. Turning off GCP."
                        )
                        prepinfo[jt]["gcp"] = False

                break

        return prepinfo

    def _print_update(self) -> None:
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

    def _print_comparison(self) -> None:
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
