"""
Screening is basically just an extension of part0 (Prescreening).
Additionally to part0 it is also possible to calculate gsolv implicitly and include the RRHO contribution.
"""
from typing import List
import os
from math import isclose

from censo.prescreening import Prescreening
from censo.part import CensoPart
from censo.utilities import print, timeit, format_data
from censo.parallel import ProcessHandler
from censo.core import CensoCore
from censo.settings import CensoSettings
from censo.datastructure import MoleculeData
from censo.cfg import (
    PLENGTH,
    AU2KCAL,
)


class Screening(Prescreening):
    
    alt_name = "part1"
    
    def __init__(self, core: CensoCore, settings: CensoSettings):
        CensoPart.__init__(self, core, settings, "screening")


    @timeit
    def run(self) -> None:
        """
        Advanced screening of the ensemble by doing single-point calculations on the input geometries,
        but this time with the ability to additionally consider implicit solvation and finite temperature contributions.

        Basically consists of two parts:
            - screening of the ensemble by doing single-point calculations on the input geometries (just as prescreening),
            - conformers are sorted out using these values and RRHO contributions are calculated (if enabled), updating the ensemble a second time
        """
        # TODO - test this
        # TODO - maybe put folder/handler as instance variable such that it can be reused later instead of reinstantiating
        super().run()

        # PART (2)
        # TODO - overwrite 'gtot'?
        threshold = self._instructions.get("threshold", None)
        
        # folder should already exist if previous part didn't raise runtime error
        folder = os.path.join(self.core.workdir, self.__class__.__name__.lower())
        if self._instructions["evaluate_rrho"]:
            # initialize process handler for current program with conformer geometries
            handler = ProcessHandler(self._instructions, [conf.geom for conf in self.core.conformers])
            
            jobtype = ["xtb_rrho"]

            # append results to previous results
            results = handler.execute(jobtype, folder)
            for conf in self.core.conformers:
                # update results for each conformer
                conf.results[self.__class__.__name__.lower()].update(results[id(conf)])
                
                # calculate new gtot including RRHO contribution
                conf.results[self.__class__.__name__.lower()]["gtot"] = self.key2(conf)

            # sort conformers list
            self.core.conformers.sort(key=lambda conf: conf.results[self.__class__.__name__.lower()]["gtot"])

            # update conformers with threshold
            # pick the free enthalpy of the first conformer as limit, since the conformer list is sorted
            limit = self.core.conformers[0].results[self.__class__.__name__.lower()]["gtot"]
            
            # filter out all conformers below threshold
            # so that 'filtered' contains all conformers that should not be considered any further
            filtered = [
                conf for conf in filter(
                    lambda x: self.key2(x) > limit + threshold, 
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

            # TODO - write results for second part
            self.write_results2()
                
        # DONE


    def key2(self, conf: MoleculeData) -> float:
        """
        Calculate the total Gibbs free energy (Gtot) of a given molecule using DFT single-point and gsolv (if included) and RRHO contributions.
        
        Parameters:
            conf (MoleculeData): The MoleculeData object containing the necessary information for the calculation.
        
        Returns:
            float: The total Gibbs free energy (Gtot) of the molecule.
        """
        # Gtot = E(DFT) + Gsolv + Grrho
        # note: key2 should only be called if evaluate_rrho is True
        return self.key(conf) + conf.results[self.__class__.__name__.lower()]["xtb_rrho"]["gibbs"][self._instructions["temperature"]]


    def write_results(self) -> None:
        """
        Overrides the write_results function of Prescreening.
        Write the results to a file in formatted way.
        writes (1):
            E (xtb),
            δE (xtb),
            E (DFT),
            δGsolv (DFT),
            Gtot,
            δGtot
            
        Generates NO csv file. All info is included in the file written in write_results2.
        """
        # PART (1) of writing
        # column headers
        headers = [
            "CONF#",
            "E (xTB)",
            "ΔE (xTB)",
            "E (DFT)",
            "ΔGsolv",
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
        # minimal xtb single-point energy (taken from prescreening)
        # TODO - where do prescreening and screening get xtb single-point from?
        xtbmin = min(
            conf.results["prescreening"]['xtb_gsolv']['energy_xtb_gas'] 
            for conf in self.core.conformers
        )
        
        # minimal total free enthalpy (single-point and potentially gsolv)
        gtotmin = min(self.key(conf) for conf in self.core.conformers)
        
        # determines what to print for each conformer in each column
        # TODO - remaining float accuracies
        printmap = {
            "CONF#": lambda conf: conf.name,
            "E (xTB)": lambda conf: f"{conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas']:.6f}", # TODO
            "ΔE (xTB)": lambda conf: f"{(conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] - xtbmin) * AU2KCAL:.2f}", # TODO
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}",
            "ΔGsolv": lambda conf: 
                f"{self.key(conf) - conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}" 
                if "xtb_gsolv" in conf.results[self.__class__.__name__.lower()].keys() or "gsolv" in conf.results[self.__class__.__name__.lower()].keys()
                else "---", 
            "Gtot": lambda conf: f"{self.key(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.key(conf) - gtotmin) * AU2KCAL:.2f}",
        }

        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)

        # write everything to a file
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.out"), "w", newline=None) as outfile:
            outfile.writelines(lines)
        
        
    def write_results2(self) -> None:
        """
        Additional write function in case RRHO is used.
        Write the results to a file in formatted way. This is appended to the first file.
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
        # column headers
        headers = [
            "CONF#",
            "G (xTB)",
            "ΔG (xTB)",
            "E (DFT)",
            "ΔGsolv",
            "GmRRHO",
            "Gtot",
            "ΔGtot",
        ]
        
        # column units
        units = [
            "",
            "[Eh]",
            "[kcal/mol]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[Eh]",
            "[kcal/mol]",
        ]

        # minimal xtb energy from single-point (and mRRHO)
        gxtbmin = min(
            conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions["temperature"]] # TODO?
            if self._instructions["evaluate_rrho"] else conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas']
            for conf in self.core.conformers
        )

        # minimal gtot from E(DFT), Gsolv and GmRRHO
        gtotmin = min(
            self.key2(conf)
            for conf in self.core.conformers
        )

        printmap = {
            "CONF#": lambda conf: conf.name,
            "G (xTB)": lambda conf: f"{conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}", # TODO
            "ΔG (xTB)": lambda conf: f"{(conf.results['prescreening']['xtb_gsolv']['energy_xtb_gas'] + conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']] - gxtbmin) * AU2KCAL:.2f}", # TODO
            "E (DFT)": lambda conf: f"{conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}",
            "ΔGsolv": lambda conf: 
                f"{self.key(conf) - conf.results[self.__class__.__name__.lower()]['sp']['energy']:.6f}"
                if not isclose(self.key(conf), conf.results[self.__class__.__name__.lower()]['sp']['energy'])
                else "---", 
            "GmRRHO": lambda conf: 
                f"{conf.results[self.__class__.__name__.lower()]['xtb_rrho']['gibbs'][self._instructions['temperature']]:.6f}"
                if self._instructions["evaluate_rrho"]
                else "---", 
            "Gtot": lambda conf: f"{self.key2(conf):.6f}",
            "ΔGtot": lambda conf: f"{(self.key2(conf) - gtotmin) * AU2KCAL:.2f}",
        }
        
        rows = [[printmap[header](conf) for header in headers] for conf in self.core.conformers]

        lines = format_data(headers, rows, units=units)
        
        # append lines to already existing file
        with open(os.path.join(self.core.workdir, f"{self.__class__.__name__.lower()}.out"), "a", newline=None) as outfile:
            outfile.writelines(lines)


"""
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
"""