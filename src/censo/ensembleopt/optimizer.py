from ..params import AU2KCAL, DIGILEN
from ..ensembledata import EnsembleData
from ..part import CensoPart
from ..utilities import (
    timeit,
    format_data,
)
from ..logging import setup_logger

logger = setup_logger(__name__)


class EnsembleOptimizer(CensoPart):
    """
    Boilerplate class for all ensemble optimization steps.
    """

    def __init__(self, ensemble: EnsembleData):
        super().__init__(ensemble)

    @timeit
    @CensoPart._create_dir
    def run(self, cut: bool = True) -> None:
        """
        Boilerplate run logic for any ensemble optimization step. The 'optimize' method should be implemented for every
        class respectively.
        """
        # print instructions
        self.print_info()

        # Print information about ensemble before optimization
        self.print_update()

        # Perform the actual optimization logic
        self.optimize(cut=cut)

        # Print comparison with previous parts
        self.print_comparison()

        # Print information about ensemble after optimization
        self.print_update()

        # dump ensemble
        self.ensemble.dump_ensemble(self._name)

        # DONE

    def print_update(self) -> None:
        print("\n\n")

        print("Number of conformers:".ljust(DIGILEN // 2, " ") +
              f"{len(self.ensemble.conformers)}")
        print("Highest ranked conformer:".ljust(DIGILEN // 2, " ") +
              f"{self.ensemble.conformers[0].name}")

        print("\n\n")

    def print_comparison(self) -> None:
        headers = [
            "CONF#"
        ]

        parts = list(self.ensemble.conformers[0].results.keys())

        headers.extend(
            [f"ΔGtot ({part})" for part in parts]
        )

        # column units
        units = [
            "",
        ]

        units.extend(["[kcal/mol]" for _ in range(len(parts))])

        # variables for printmap
        gtotmin = {
            part: 0.0 for part in parts
        }
        for part in parts:
            gtotmin[part] = min(conf.results[part]["gtot"]
                                for conf in self.ensemble.conformers)

        # determines what to print for each conformer in each column
        printmap = {
            "CONF#": lambda conf: conf.name,
        }
        printmap.update({
            header: lambda conf: f"{(conf.results[part]['gtot'] - gtotmin[part]) * AU2KCAL:.2f}"
            for header, part in zip(headers[1:], parts)
        })

        rows = [
            [printmap[header](conf) for header in headers]
            for conf in self.ensemble.conformers
        ]

        lines = format_data(headers, rows, units=units)

        # Print everything
        for line in lines:
            print(line, flush=True, end="")
