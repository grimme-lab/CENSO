"""
prescreening == part1, calculate free energy on GFNn-xTB input geometry
idea is to improve on E and (Gsolv)
"""
from censo.part import CensoPart
from censo.utilities import print, timeit, format_data


class Screening(CensoPart):
    def __init__(self, core: CensoCore, settings: CensoSettings):
        super().__init__(core, settings, "screening")


    @timeit
    def run(self) -> None:
        """
        Advanced screening of the ensemble by doing single-point calculations on the input geometries,
        but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

        Basically consists of two parts:
            - screening of the ensemble by doing single-point calculations on the input geometries,
            - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
        """
        # initialize process handler for current program with conformer geometries
        handler = ProcessHandler(self._instructions, [conf.geom for conf in self.core.conformers])

        # print instructions
        self.print_info()

        # PART (1) 
        # set jobtype to pass to handler
        jobtype: List[str] = []
        if self._instructions.get("gas-phase", None):
            if self._instructions.get("implicit", None):
                jobtype = ["sp", "gsolv"]
            else:
                jobtype = ["sp", "xtb_gsolv"]
        else:
            jobtype = ["sp"]
        
        # set folder to do the calculations in
        folder = os.path.join(self.core.workdir, 'screening')
        if os.path.isdir(folder):
            # TODO - warning? stderr?
            print(f"Folder {folder} already exists. Potentially overwriting files.")
        elif os.system(f"mkdir {folder}") != 0 and not os.path.isdir(folder):
            # TODO - stderr case?
            print(f"Could not create directory for {self.__class__.__name__.lower()}. Executing calculations in {self.core.workdir}.")
            folder = self.core.workdir
        
        # compute results
        # for structure of results from handler.execute look there
        results = handler.execute(jobtype, folder)

        # update results for each conformer
        for conf in self.core.conformers:
            # store results of single jobs for each conformer
            conf.results[self.__class__.__name__.lower()] = results[id(conf)]
            
        # sort conformers list
        self.core.conformers.sort()
        
        # update conformers with threshold
        threshold = self._instructions.get("threshold", None)
        if not threshold is None:
            # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
            limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
            
            # filter out all conformers below threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.key(x) > limit + threshold, 
                    self.core.conformers
                )
            ]
            
            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(filtered)  
        else:
            """
            TODO
            print warning that no threshold could be determined
            (should not happen but doesn't necessarily break the program)
            """
            print("...")

        # PART (2)
        if self._instructions["evaluate_rrho"]:
            jobtype = ["xtb_rrho"]

            # append results to previous results
            results = handler.execute(jobtype, folder)
            for conf in self.core.conformers:
                conf.results[self.__class__.__name__.lower()].update(results[id(conf)])

            # sort conformers list
            self.core.conformers.sort()

            # update conformers with threshold
            # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
            limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
            
            # filter out all conformers below threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.key(x) > limit + threshold, 
                    self.core.conformers
                )
            ]
            
            # update the conformer list in core (remove conf if below threshold)
            self.core.update_conformers(filtered)

            # calculate boltzmann weights from gtot values calculated here
            # trying to get temperature from instructions, set it to room temperature if that fails for some reason
            self.core.calc_boltzmannweights(
                self._instructions.get("temperature", 298.15),
                self.__class__.__name__.lower()
            )
            
            # if no conformers are filtered basically nothing happens
        
        # write results (analogous to deprecated print)
        self.write_results()
                
        # DONE


    def key(self, conf: CensoConformer) -> float:
        """
        Calculates the key value for the given CensoConformer object.

        Args:
            conf (CensoConformer): The CensoConformer object for which to calculate the key.

        Returns:
            float: The calculated key value.
        """
        pass


    def print_info(self) -> None:
        """
        A function that prints information.
        """
        pass


    def write_results(self) -> None:
        """
        Write the results to a file in formatted way.
        writes (1):
            E (xtb),
            δE (xtb),
            E (DFT),
            δGsolv (DFT),
            Gtot,
            δGtot
        
        writes (2):
            G (xtb),
            δG (xtb),
            E (DFT),
            δGsolv (DFT),
            Grrho,
            Gtot,
            δGtot
        
        Also writes them into an easily digestible format.
        """
        # PART (1) of writing
        # column headers
        headers = [
            "CONF#",
            "E (xtb)",
            "ΔE (xtb)",
            "E (DFT)",
            "ΔGsolv (DFT)",
            "Gtot",
            "ΔGtot",
        ]
        
        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[kcal/mol]",
        ]
        
        # variables for printmap
        # minimal xtb single-point energy
        xtbmin = min(
            conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas'] 
            for conf in self.core.conformers
        )
        
        # minimal dft single-point energy
        dftmin = min(
            conf.results[self.__class__.__name__.lower()]['sp']['energy'] 
            for conf in self.core.conformers
        )
        
        # minimal solvation free enthalpy
        gsolvmin = min(
            conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv'] 
            if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() 
            else 0.0 
            for conf in self.core.conformers
        ) 
        
        # minimal total free enthalpy
        gtotmin = min(self.key(conf) for conf in self.core.conformers)
        
        # determines what to print for each conformer in each column
        # TODO - remaining float accuracies
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xtb)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas']}",
            "ΔE (xtb)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}",
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']}",
            "δGsolv (xtb)": lambda conf: 
                f"{conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv']}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys()
                else "---",
            "Gtot": lambda conf: f"{self.key(conf)}",
            "ΔE (DFT)": lambda conf: f"{(conf.results[self.__class__.__name__.lower()]['sp']['energy'] - dftmin) * AU2KCAL:.2f}",
            "ΔδGsolv": lambda conf: 
                f"{(conf.results[self.__class__.__name__.lower()]['xtb_gsolv']['gsolv'] - gsolvmin) * AU2KCAL:.2f}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() 
                else "---",
            "ΔGtot": lambda conf: f"{(self.key(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows1 = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows1, units=units)

        # list all conformers still considered with their boltzmann weights
        # as well as the averaged free enthalpy of the ensemble
        lines.append(
            "\nBoltzmann averaged free energy/enthalpy of ensemble on input geometries (not DFT optimized):\n"
        )
        lines.append(f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} {'avG(T) /a.u.':>14}\n")
        print("".ljust(int(PLENGTH), "-") + "\n")

        # calculate averaged free enthalpy
        avG = sum([
            conf.results[self.__class__.__name__.lower()]["bmw"] 
            * conf.results[self.__class__.__name__.lower()]["gtot"] 
            for conf in self.core.conformers
        ])
        
        # calculate averaged free energy
        avE = sum([
            conf.results[self.__class__.__name__.lower()]["bmw"]
            * conf.results[self.__class__.__name__.lower()]["sp"]["energy"]
            for conf in self.core.conformers
        ])

        # append the lines for the free energy/enthalpy
        lines.append(f"{self._instructions.get('temperature', 298.15):^15} {avE:>14.7f}  {avG:>14.7f}     <<==part1==\n")
        lines.append("".ljust(int(PLENGTH), "-") + "\n\n")
        
        lines.append(">>> END of Screening <<<".center(PLENGTH, " ") + "\n")
        
        # PART (2) of writing


def part1(config, conformers, store_confs, ensembledata):
    # RANDOM STUFF - configuring jobs
    # GENERAL CONFIGURATION
    # IN SHORT:
    # do single point, if in gas-phase use no solvation
    # then there is the option to either use additive gsolv correction as in prescreening
    # or use implicit solvation directly in orca or tm
    # GETTING RESULTS
    # SORTING AND STUFF
    # ***************************************************************************
    # first sorting by E or Gsolv
    # (remove high lying conformers above part1_threshold)
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("Removing high lying conformers".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")

    for conf in calculate:
        rrho = None
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        e = "prescreening_sp_info"
        conf.calc_free_energy(
            e=e,
            solv=solv,
            rrho=rrho,
            t=config.fixed_temperature,
            consider_sym=config.consider_sym,
        )
    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
    except ValueError:
        raise
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
    try:
        maxreldft = max([i.rel_free_energy for i in calculate if i is not None])
    except ValueError:
        print_errors(
            f"{'ERROR:':{WARNLEN}}No conformer left or Error in maxreldft!", save_errors
        )
    # print sorting
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "xtb_energy"),
        lambda conf: getattr(conf, "rel_xtb_energy"),
        lambda conf: getattr(conf, "prescreening_sp_info")["energy"],
        lambda conf: getattr(conf, "prescreening_gsolv_info")
        .get("range", {})
        .get(config.fixed_temperature, 0.0),
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "rel_free_energy"),
    ]
    columnheader = [
        "CONF#",
        "E(GFNn-xTB)",
        "ΔE(GFNn-xTB)",
        "E [Eh]",
        "Gsolv [Eh]",
        "gtot",
        "Δgtot",
    ]
    columndescription = ["", "[a.u.]", "[kcal/mol]", "", "", "[Eh]", "[kcal/mol]"]
    columndescription2 = ["", "", "", "", "", "", "", ""]
    columnformat = ["", (12, 7), (5, 2), (12, 7), (12, 7), (12, 7), (5, 2)]

    if config.solvent == "gas":
        columnheader[5] = "Etot"
        columnheader[6] = "ΔEtot"
        columndescription[3] = instruction["method"]
        # ignore gsolv in printout
        columncall.pop(4)
        columnheader.pop(4)
        columndescription.pop(4)
        columnformat.pop(4)
    elif config.solvent != "gas":
        # energy
        columndescription[3] = instruction["method"]
        # gsolv
        columndescription[4] = instruction["method2"]

    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part1preG.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
        columndescription2=columndescription2,
    )

    if maxreldft > (config.part1_threshold):
        print("\n" + "".ljust(int(PLENGTH / 2), "-"))
        print("Conformers considered further".center(int(PLENGTH / 2), " "))
        print("".ljust(int(PLENGTH / 2), "-") + "\n")
        for conf in list(calculate):
            if conf.rel_free_energy > (config.part1_threshold):
                store_confs.append(calculate.pop(calculate.index(conf)))
        if calculate:
            print(
                f"Below the g_thr(1) threshold of {config.part1_threshold} kcal/mol.\n"
            )
            print_block(["CONF" + str(i.id) for i in calculate])
        else:
            print(f"{'ERROR:':{WARNLEN}}There are no more conformers left!")
    else:
        print(
            "\nAll relative (free) energies are below the g_thr(1) threshold "
            f"of {config.part1_threshold} kcal/mol.\nAll conformers are "
            "considered further."
        )
    ensembledata.nconfs_per_part["part1_firstsort"] = len(calculate)
    # reset
    for conf in calculate:
        conf.free_energy = 0.0
        conf.rel_free_energy = None
    print("".ljust(int(PLENGTH / 2), "-"))

    # DO RRHO (if wanted (general -> evaluate_rrho))

    # # printout for part1 -------------------------------------------------------
    print("\n" + "".ljust(int(PLENGTH / 2), "-"))
    print("* Gibbs free energies of part1 *".center(int(PLENGTH / 2), " "))
    print("".ljust(int(PLENGTH / 2), "-") + "\n")
    columncall = [
        lambda conf: "CONF" + str(getattr(conf, "id")),
        lambda conf: getattr(conf, "xtb_free_energy"),
        lambda conf: getattr(conf, "rel_xtb_free_energy"),
        lambda conf: getattr(conf, "prescreening_sp_info")["energy"],
        lambda conf: getattr(conf, "prescreening_gsolv_info")
        .get("range", {})
        .get(config.fixed_temperature, 0.0),
        lambda conf: conf.get_mrrho(
            config.fixed_temperature, "prescreening_grrho_info", config.consider_sym
        ),
        lambda conf: getattr(conf, "free_energy"),
        lambda conf: getattr(conf, "rel_free_energy"),
    ]
    columnheader = [
        "CONF#",
        "G(GFNn-xTB)",
        "ΔG(GFNn-xTB)",
        "E [Eh]",
        "Gsolv [Eh]",
        "GmRRHO [Eh]",
        "Gtot",
        "ΔGtot",
    ]
    columndescription = [
        "",  # CONFX
        "[a.u.]",  # xtb energy
        "[kcal/mol]",  # rel xtb_energy
        str(config.func),  # E
        "",  # GSolv
        "",
        "[Eh]",  # Gtot
        "[kcal/mol]",  # rel Gtot
    ]
    columndescription2 = ["", "", "", "", "", "", "", ""]
    columnformat = ["", (12, 7), (5, 2), (12, 7), (12, 7), (12, 7), (12, 7), (5, 2)]
    if config.solvent == "gas":
        # Energy
        columndescription[3] = instruction["method"]
    elif config.solvent != "gas":
        # Energy
        columndescription[3] = instruction["method"]
        # Gsolv
        columndescription[4] = instruction["method2"]
    if config.evaluate_rrho:
        columndescription[5] = instruction_prerrho["method"]  # Grrho

    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho:
            # ignore rrho in printout
            columncall.pop(5)
            columnheader.pop(5)
            columndescription.pop(5)
            columnformat.pop(5)
        if config.solvent == "gas":
            columncall.pop(4)
            columnheader.pop(4)
            columndescription.pop(4)
            columnformat.pop(4)

    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "prescreening_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        conf.calc_free_energy(
            e="prescreening_sp_info",
            solv=solv,
            rrho=rrho,
            t=config.fixed_temperature,
            consider_sym=config.consider_sym,
        )
        conf.xtb_free_energy = conf.calc_free_energy(
            e="xtb_energy",
            solv=None,
            rrho=rrho,
            out=True,
            t=config.fixed_temperature,
            consider_sym=config.consider_sym,
        )

    try:
        minfree = min([i.free_energy for i in calculate if i is not None])
        minfree_xtb = min([i.xtb_free_energy for i in calculate if i is not None])
    except ValueError:
        raise ValueError
    for conf in calculate:
        conf.rel_free_energy = (conf.free_energy - minfree) * AU2KCAL
        conf.rel_xtb_free_energy = (conf.xtb_free_energy - minfree_xtb) * AU2KCAL
    try:
        maxreldft = max([i.rel_free_energy for i in calculate if i is not None])
    except ValueError:
        print_errors(
            f"{'ERROR:':{WARNLEN}}No conformer left or Error in maxreldft!", save_errors
        )
    # print sorting
    calculate.sort(key=lambda x: int(x.id))
    printout(
        os.path.join(config.cwd, "part1.dat"),
        columncall,
        columnheader,
        columndescription,
        columnformat,
        calculate,
        minfree,
        columndescription2=columndescription2,
    )
    # --------------------------------------------------------------------------
    calculate = calc_boltzmannweights(
        calculate, "free_energy", config.fixed_temperature
    )
    conf_in_interval(calculate)
    # --------------------------------------------------------------------------
    for conf in calculate:
        if conf.free_energy == minfree:
            ensembledata.bestconf["part1"] = conf.id

    # write to enso.json
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    # ***************************************************************************
    # Fuzzy or smart sorting
    # increase the individual threshold for conformers with GRRHO differing from
    # the mean GmRRHO
    if len(calculate) == 1:
        std_dev = 0.0
    else:
        std_dev = calc_std_dev(
            [
                conf.get_mrrho(
                    config.fixed_temperature,
                    rrho="prescreening_grrho_info",
                    consider_sym=config.consider_sym,
                )
                * AU2KCAL
                for conf in calculate
                if conf.get_mrrho(
                    config.fixed_temperature,
                    rrho="prescreening_grrho_info",
                    consider_sym=config.consider_sym,
                )
                is not None
            ]
        )
    max_fuzzy = 1
    fuzzythr = max_fuzzy * (1 - math.exp(-1 * 5 * (std_dev ** 2)))
    print(
        "\nAdditional global 'fuzzy-threshold' based on the standard deviation of (G_mRRHO):"
    )
    print(f"Std_dev(G_mRRHO) = {std_dev:.3f} kcal/mol")
    print(f"Fuzzythreshold   = {fuzzythr:.3f} kcal/mol")
    print(
        f"Final sorting threshold G_thr(1) = {config.part1_threshold:.3f} + "
        f"{fuzzythr:.3f} = {config.part1_threshold + fuzzythr:.3f} kcal/mol"
    )
    for conf in calculate:
        conf.prescreening_grrho_info["fuzzythr"] = fuzzythr

    # spearman between DFT and DFT + RRHO
    if config.evaluate_rrho and len(calculate) > 1:
        for conf in calculate:
            rrho = None
            if config.solvent == "gas":
                solv = None
            else:
                solv = "prescreening_gsolv_info"
            e = "prescreening_sp_info"
            conf.calc_free_energy(
                e=e,
                solv=solv,
                rrho=rrho,
                t=config.fixed_temperature,
                consider_sym=config.consider_sym,
            )
        try:
            minfree = min([i.free_energy for i in calculate if i is not None])
        except ValueError:
            raise ValueError
        without_RRHO = []
        calculate.sort(key=lambda x: int(x.id))
        for conf in calculate:
            without_RRHO.append((conf.free_energy - minfree) * AU2KCAL)
        for conf in calculate:
            conf.free_energy = 0.0
        for conf in calculate:
            rrho = "prescreening_grrho_info"
            if config.solvent == "gas":
                solv = None
            else:
                solv = "prescreening_gsolv_info"
            e = "prescreening_sp_info"
            conf.calc_free_energy(
                e=e,
                solv=solv,
                rrho=rrho,
                t=config.fixed_temperature,
                consider_sym=config.consider_sym,
            )
        try:
            minfree = min([i.free_energy for i in calculate if i is not None])
        except ValueError:
            raise ValueError
        with_RRHO = []
        calculate.sort(key=lambda x: int(x.id))
        for conf in calculate:
            with_RRHO.append((conf.free_energy - minfree) * AU2KCAL)
        for conf in calculate:
            conf.free_energy = 0.0
        if config.solvent != "gas":
            print(
                f"Spearman correlation coefficient between (E + Solv) "
                f"and (E + Solv + mRRHO) = {spearman(without_RRHO, with_RRHO):.3f}"
            )
        else:
            print(
                f"Spearman correlation coefficient between (E) "
                f"and (E + mRRHO) = {spearman(without_RRHO, with_RRHO):.3f}"
            )

    # sorting
    if maxreldft > config.part1_threshold:
        print("\n" + "".ljust(int(PLENGTH / 2), "-"))
        print("Conformers considered further".center(int(PLENGTH / 2), " "))
        print("".ljust(int(PLENGTH / 2), "-") + "\n")
        for conf in list(calculate):
            if conf.rel_free_energy <= config.part1_threshold:
                conf.part_info["part1"] = "passed"
            elif conf.rel_free_energy <= (
                config.part1_threshold + conf.prescreening_grrho_info["fuzzythr"]
            ):
                print(f"Considered CONF{conf.id} because of increased fuzzythr.")
                conf.part_info["part1"] = "passed"
                continue
            else:
                conf.part_info["part1"] = "refused"
                store_confs.append(calculate.pop(calculate.index(conf)))
        if calculate:
            print(
                f"These conformers are below the {config.part1_threshold+fuzzythr:.3f} "
                f"kcal/mol G_thr(1) threshold.\n"
            )
            print_block(["CONF" + str(i.id) for i in calculate])
        else:
            print_errors(
                f"{'ERROR:':{WARNLEN}}There are no more conformers left!", save_errors
            )
    else:
        for conf in list(calculate):
            conf.part_info["part1"] = "passed"
        print(
            "\nAll relative (free) energies are below the initial G_thr(1) threshold "
            f"of {config.part1_threshold} kcal/mol.\nAll conformers are "
            "considered further."
        )
    ensembledata.nconfs_per_part["part1"] = len(calculate)
    # save current data to jsonfile
    config.write_json(
        config.cwd,
        [i.provide_runinfo() for i in calculate]
        + [i.provide_runinfo() for i in prev_calculated]
        + [i.provide_runinfo() for i in store_confs]
        + [ensembledata],
        config.provide_runinfo(),
    )

    # free energy:
    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "prescreening_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        conf.calc_free_energy(
            e="prescreening_sp_info",
            solv=solv,
            rrho=rrho,
            t=config.fixed_temperature,
            consider_sym=config.consider_sym,
        )

    # write coord.enso_best
    for conf in calculate:
        if conf.id == ensembledata.bestconf["part1"]:
            # copy the lowest optimized conformer to file coord.enso_best
            with open(
                os.path.join("CONF" + str(conf.id), config.func, "coord"),
                "r",
                encoding=CODING,
                newline=None,
            ) as f:
                coord = f.readlines()
            with open(
                os.path.join(config.cwd, "coord.enso_best"), "w", newline=None
            ) as best:
                best.write(
                    "$coord  # {}   {}   !CONF{} \n".format(
                        conf.free_energy,
                        conf.get_mrrho(
                            config.fixed_temperature,
                            rrho="prescreening_grrho_info",
                            consider_sym=config.consider_sym,
                        ),
                        conf.id,
                    )
                )
                for line in coord[1:]:
                    if "$" in line:  # stop at $end ...
                        break
                    best.write(line)
                best.write("$end \n")

    ################################################################################
    # calculate average G correction
    print(
        "\nCalculating Boltzmann averaged free energy of ensemble on "
        f"input geometries (not DFT optimized)!\n"
    )
    # calculate Boltzmannweights
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho and config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avG(T) /a.u.':>14} "
            )
        elif not config.evaluate_rrho:
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGsolv(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
            )
        elif config.solvent == "gas":
            line = (
                f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
                f"{'avGmRRHO(T) /a.u.':>16} {'avG(T) /a.u.':>14} "
            )
    else:
        line = (
            f"{'temperature /K:':<15} {'avE(T) /a.u.':>14} "
            f"{'avGmRRHO(T) /a.u.':>16} {'avGsolv(T) /a.u.':>16} "
            f"{'avG(T) /a.u.':>14}"
        )
    print(line)
    print("".ljust(int(PLENGTH), "-"))
    # get free energy at (T)
    for conf in calculate:
        if not config.evaluate_rrho:
            rrho = None
        else:
            rrho = "prescreening_grrho_info"
        if config.solvent == "gas":
            solv = None
        else:
            solv = "prescreening_gsolv_info"
        e = "prescreening_sp_info"
        conf.calc_free_energy(
            e=e,
            solv=solv,
            rrho=rrho,
            t=config.fixed_temperature,
            consider_sym=config.consider_sym,
        )

    calculate = calc_boltzmannweights(
        calculate, "free_energy", config.fixed_temperature
    )
    avG = 0.0
    avE = 0.0
    avGRRHO = 0.0
    avGsolv = 0.0
    for conf in calculate:
        avG += conf.bm_weight * conf.free_energy
        avE += conf.bm_weight * conf.prescreening_sp_info["energy"]
        # avGRRHO += conf.bm_weight * conf.prescreening_grrho_info["energy"]
        avGRRHO += conf.bm_weight * conf.get_mrrho(
            config.fixed_temperature, "prescreening_grrho_info", config.consider_sym
        )
        avGsolv += conf.bm_weight * conf.prescreening_gsolv_info["range"].get(
            config.fixed_temperature, 0.0
        )

    # printout:
    if not config.evaluate_rrho or config.solvent == "gas":
        if not config.evaluate_rrho and config.solvent == "gas":
            line = f"{config.fixed_temperature:^15} {avE:>14.7f}  {avG:>14.7f} "
        elif not config.evaluate_rrho:
            line = (
                f"{config.fixed_temperature:^15} {avE:>14.7f} {avGsolv:>16.7f} "
                f"{avG:>14.7f} "
            )
        elif config.solvent == "gas":
            line = (
                f"{config.fixed_temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
                f"{avG:>14.7f} "
            )
    else:
        line = (
            f"{config.fixed_temperature:^15} {avE:>14.7f} {avGRRHO:>16.7f} "
            f"{avGsolv:>16.7f} {avG:>14.7f} "
        )
    print(line, "    <<==part1==")
    print("".ljust(int(PLENGTH), "-"))
    print("")

    ############################################################################
    # Calculate unbiased GFNn-xTB energy:
    for conf in list(calculate):
        if conf.xtb_energy_unbiased is None:
            pass
        else:
            conf = calculate.pop(calculate.index(conf))
            conf.job["success"] = True
            prev_calculated.append(conf)

    if calculate:
        print("\nCalculating unbiased GFNn-xTB energy")
        instruction_gfn = {
            "jobtype": "xtb_sp",
            "func": getattr(config, "part1_gfnv"),
            "charge": config.charge,
            "unpaired": config.unpaired,
            "solvent": config.solvent,
            "sm": config.sm_rrho,
            "rmsdbias": config.rmsdbias,
            "temperature": config.fixed_temperature,
            "gfn_version": config.part1_gfnv,
            "energy": 0.0,
            "energy2": 0.0,
            "success": False,
            "onlyread": config.onlyread,
        }
        folder_gfn = "GFN_unbiased"
        save_errors, store_confs, calculate = new_folders(
            config.cwd, calculate, folder_gfn, save_errors, store_confs
        )
        # write coord to folder
        calculate, store_confs, save_errors = ensemble2coord(
            config, folder_gfn, calculate, store_confs, save_errors
        )
        calculate = run_in_parallel(
            config,
            q,
            resultq,
            job,
            config.maxthreads,
            config.omp,
            calculate,
            instruction_gfn,
            config.balance,
            folder_gfn,
        )
        for conf in list(calculate):
            if not conf.job["success"]:
                conf.xtb_energy_unbiased = conf.xtb_energy
            else:
                conf.xtb_energy_unbiased = conf.job["energy"]
    if prev_calculated:
        for conf in list(prev_calculated):
            calculate.append(prev_calculated.pop(prev_calculated.index(conf)))
    ############################################################################

    # write ensemble
    move_recursively(config.cwd, "enso_ensemble_part1.xyz")
    if config.evaluate_rrho:
        kwargs = {"energy": "xtb_energy_unbiased", "rrho": "prescreening_grrho_info"}
    else:
        kwargs = {"energy": "xtb_energy_unbiased"}
    write_trj(
        sorted(calculate, key=lambda x: float(x.free_energy)),
        config.cwd,
        "enso_ensemble_part1.xyz",
        config.func,
        config.nat,
        "free_energy",
        config.fixed_temperature,
        config.consider_sym,
        **kwargs,
    )

    # reset
    for conf in calculate:
        conf.free_energy = 0.0
        conf.rel_free_energy = None
        conf.bm_weight = 0.0
        conf.reset_job_info()
    if save_errors:
        print("\n***---------------------------------------------------------***")
        print("Printing most relevant errors again, just for user convenience:")
        for _ in list(save_errors):
            print(save_errors.pop())
        print("***---------------------------------------------------------***")

    tmp = int((PLENGTH - len("END of Part1")) / 2)
    print("\n" + "".ljust(tmp, ">") + "END of Part1" + "".rjust(tmp, "<"))
    if config.progress:
        print("#>>># CENSO: Finished part1", file=sys.stderr)
    return config, calculate, store_confs, ensembledata
