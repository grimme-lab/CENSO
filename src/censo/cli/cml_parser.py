"""
defininition of internal defaults, checking of logic for parameter combinations,
cml parsing
"""

from ..params import START_DESCR
import argparse


def check_soft_requirements(args: argparse.Namespace) -> bool:
    """
    Checks for soft-required options (e.g. if you call --new-config you don't need to give -i, for a normal
    CENSO run you need to, though).
    """
    soft_required = [
        "inp",
        "maxcores",
        # "charge",
        # "unpaired",
    ]
    requirement_override = ["writeconfig", "cleanup", "cleanup_all", "version"]
    # If all settings, that override soft-requirement, are unused
    if all(getattr(args, s, None) is False for s in requirement_override):
        # Check, if all the soft-required settings are given
        if all(getattr(args, s, None) is not None for s in soft_required):
            return True
        # Else, the check fails
        else:
            return False
    # Else, the the soft requirement is overridden
    else:
        return True


def parse(argv=None) -> argparse.Namespace:
    """
    Process commandline arguments

    NOTE: on args with the action 'store_const' with const=True, this is on purpose so as long as the flag is not set,
    the arg Namespace evaluates to None.
    """

    parser = argparse.ArgumentParser(
        description=START_DESCR,
        prog="censo",
    )

    groups = []

    # RUN SETTINGS
    groups.append(parser.add_argument_group("RUN SETTINGS"))
    groups[0].add_argument(
        "-i",
        "--input",
        dest="inp",
        type=str,
        help="Relative path to ensemble file, e.g. crest_conformers.xyz. For a default run this is REQUIRED. ",
    )
    groups[0].add_argument(
        "-n",
        "--nconf",
        dest="nconf",
        type=int,
        help="The first 'nconf' conformers will be considered.",
    )
    groups[0].add_argument(
        "-c",
        "--charge",
        dest="charge",
        default=0,
        type=int,
        help="Integer charge of the investigated molecule.",
    )
    groups[0].add_argument(
        "-u",
        "--unpaired",
        dest="unpaired",
        default=0,
        type=int,
        help="Integer number of unpaired electrons of the investigated molecule.",
    )
    groups[0].add_argument(
        "-v",
        "--version",
        dest="version",
        action="store_true",
        help="Print CENSO version and exit.",
    )
    groups[0].add_argument(
        "--cleanup",
        dest="cleanup",
        action="store_true",
        help="Delete unneeded files from current working directory.",
    )
    groups[0].add_argument(
        "--cleanup_all",
        dest="cleanup_all",
        action="store_true",
        help="Delete all CENSO files from previous runs from current working directory. "
        "Stronger than -cleanup !",
    )
    groups[0].add_argument(
        "--new-config",
        dest="writeconfig",
        action="store_true",
        help="Write new configuration file, which is placed into the current "
        "directory.",
    )
    groups[0].add_argument(
        "--inprc",
        dest="inprcpath",
        help="Use to provide a path to the CENSO configuration file if you want to use a different one"
        " than the default (~/.censo2rc).",
    )
    groups[0].add_argument(
        "--maxcores",
        dest="maxcores",
        type=int,
        help="Number of cores that should be used for CENSO on the machine. If this is not provided CENSO will use "
        "the maximum number available. For a default run this is REQUIRED.",
    )
    groups[0].add_argument(
        "-O",
        "--omp",
        dest="omp",
        type=int,
        help="Number of OpenMP threads, e.g. 4. Effectively translates to the number of cores used per calculation "
        "if load balancing is disabled.",
    )
    groups[0].add_argument(
        "--loglevel",
        dest="loglevel",
        help="Set the loglevel for all modules to a specified level.",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    groups[0].add_argument(
        "--reload",
        dest="reload",
        nargs="+",
        help="Reload data from json output files. List all file names separated by spaces. "
        "Note that all conformers from the current ensemble need to be included in the output data keys.",
    )

    # GENERAL SETTINGS
    groups.append(parser.add_argument_group("GENERAL SETTINGS"))
    groups[1].add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        type=float,
        help="Temperature in Kelvin for thermostatistical evaluation.",
    )
    groups[1].add_argument(
        "--trange",
        dest="trange",
        nargs=3,
        metavar=("start", "end", "step"),
        type=float,
        help="specify a temperature range [start, end, step] e.g.: 250.0 300.0 10.0"
        "  resulting in the range [250.0, 260.0, 270.0, 280.0, 290.0, 300.0].",
    )
    groups[1].add_argument(
        "--bhess",
        dest="bhess",
        action="store_const",
        const=True,
        help="Uses SPH and applies structure constraint to input/DFT geometry "
        "for mRRHO calcuation. ",
    )
    groups[1].add_argument(
        "--consider-sym",
        dest="consider_sym",
        action="store_const",
        const=True,
        help="Consider symmetry in mRRHO calcuation (based on desy xtb threshold). ",
    )
    groups[1].add_argument(
        "--rmsdbias",
        dest="rmsdbias",
        action="store_const",
        const=True,
        help="Applies constraint to rmsdpot.xyz to be consistent to CREST. ",
    )
    groups[1].add_argument(
        "--sm-rrho",
        dest="sm_rrho",
        type=str,
        help="Solvation model used in xTB GmRRHO calculation. Applied if not in "
        "gas-phase. Options are 'gbsa' or 'alpb'.",
    )
    groups[1].add_argument(
        "--evaluate-rrho",
        dest="evaluate_rrho",
        action="store_const",
        const=True,
        help="Evaluate mRRHO contribution.",
    )
    groups[1].add_argument(
        "-s",
        "--solvent",
        dest="solvent",
        type=str,
        help="Solvent to be used for Gsolv calculation.",
    )
    groups[1].add_argument(
        "--gas-phase",
        dest="gas-phase",
        action="store_const",
        const=True,
        help="Run calculation in gas-phase, overriding all solvation settings.",
    )
    groups[1].add_argument(
        "--imagthr",
        dest="imagthr",
        type=float,
        help="threshold for inverting imaginary frequencies for thermo in cm-1,"
        " e.g. -30.0.",
    )
    groups[1].add_argument(
        "--sthr",
        dest="sthr",
        type=float,
        help="Rotor cut-off for thermo in cm-1, e.g. 50.0.",
    )
    groups[1].add_argument(
        "--scale",
        dest="scale",
        type=float,
        help="Scaling factor for frequencies, e.g. 1.0.",
    )
    """
    groups[1].add_argument(
        "--vapor_pressure",
        "-vp",
        dest="vapor_pressure",
        action="store_true",
        help="Gsolv is evaluated for the input molecule in its solution (same). "
        "Only possible with COSMO-RS.",
    )
    """

    # PRESCREENING SETTINGS
    groups.append(parser.add_argument_group("PRESCREENING SETTINGS"))

    # SCREENING SETTINGS
    groups.append(parser.add_argument_group("SCREENING SETTINGS"))

    # OPTIMIZATION SETTINGS
    groups.append(parser.add_argument_group("OPTIMIZATION SETTINGS"))

    # REFINEMENT SETTINGS
    groups.append(parser.add_argument_group("REFINEMENT SETTINGS"))

    # NMR SETTINGS
    groups.append(parser.add_argument_group("NMR SETTINGS"))

    # OPTROT SETTINGS
    groups.append(parser.add_argument_group("OPTROT SETTINGS"))

    # UVVIS SETTINGS
    groups.append(parser.add_argument_group("UVVIS SETTINGS"))

    # leave these options out for now, implementation for cml complicated
    """
    groups[7].add_argument(
        "-freqOR",
        "--freqOR",
        dest="freq_or",
        nargs="*",
        required=False,
        type=float,
        help="Frequencies to evaluate specific rotation at in nm, e.g. 589 "
        "or 589 700 to evaluate at 598 nm and 700 nm.",
    )
    groups[6].add_argument(
        "-couplings",
        "--couplings",
        dest="couplings",
        action="store_true",
        required=False,
        help="Option to run coupling constant calculations. Options are ???.",
    )
    groups[6].add_argument(
        "-shieldings",
        "--shieldings",
        dest="shieldings",
        action="store_true",
        required=False,
        help="Option to run shielding constant calculations. Options are ???.",
    )
    groups[6].add_argument(
        "-hactive",
        "--hactive",
        dest="h_active",
        action="store_true",
        required=False,
        help="Investigates hydrogen nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-cactive",
        "--cactive",
        dest="c_active",
        action="store_true",
        required=False,
        help="Investigates carbon nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-factive",
        "--factive",
        dest="f_active",
        action="store_true",
        required=False,
        help="Investigates fluorine nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-siactive",
        "--siactive",
        dest="si_active",
        action="store_true",
        required=False,
        help="Investigates silicon nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-pactive",
        "--pactive",
        dest="p_active",
        action="store_true",
        required=False,
        help="Investigates phosophorus nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[4].add_argument(
        "-crestcheck",
        "--crestcheck",
        dest="crestcheck",
        action="store_true",
        required=False,
        
        help="Option to sort out conformers after DFT ensembleopt which CREST "
        "identifies as identical or rotamers of each other. \nThe identification/"
        "analysis is always performed, but the removal of conformers has to "
        "be the choice of the user. Options are: [???]", # TODO
    )
    groups[4].add_argument(
        "-macro",
        dest="macrocycles",
        action="store_const",
        const=True,
        required=False,
        help="Option to use macrocycles for geometry optimization."
    )
    groups[4].add_argument(
        "-optlevel2",
        "--optlevel2",
        dest="optlevel2",
        default=None,
        required=False,
        help="Option to set the optlevel in part2, only if optimizing with the xTB-optimizer!"
        "Allowed values are ***", # TODO
    )
    groups[4].add_argument(
        "-optcycles",
        "--optcycles",
        dest="optcycles",
        required=False,
        type=int,
        help="number of cycles in ensemble optimizer.",
    )
    groups[4].add_argument(
        "-hlow",
        "--hlow",
        dest="hlow",
        required=False,
        type=float,
        help="Lowest force constant in ANC generation (real), used by xTB-optimizer.",
    )
    groups[4].add_argument(
        "-spearmanthr",
        "--spearmanthr",
        dest="spearmanthr",
        required=False,
        help="Value between -1 and 1 for the spearman correlation coeffient threshold, "
        "employed in the ensemlbe optimizer",
    )
    groups[4].add_argument(
        "-radsize",
        "--radsize",
        dest="radsize",
        required=False,
        type=int,
        help="Radsize used in the ensembleopt and only for r2scan-3c!",
    )
    group1.add_argument(
        "-func",
        "--functional",
        dest="func",
        choices=options.value_options["func"],
        action="store",
        required=False,
        
        help="Functional for geometry ensembleopt (used in part2) and "
        "single-points in part1",
    ) 
    group1.add_argument(
        "-basis",
        "--basis",
        dest="basis",
        action="store",
        required=False,
        
        help="Basis set employed together with the functional (func) for the "
        "low level single point in part1 und ensembleopt in part2.",
    ) 
    group1.add_argument(
        "-prog",
        "--prog",
        choices=options.value_options["prog"],
        dest="prog",
        required=False,
        
        help="QM-program used in part0, part1 and part2 either 'orca' or 'tm'.",
    ) 
    group10.add_argument(
        "-part0_gfnv",
        "--part0_gfnv",
        dest="part0_gfnv",
        choices=options.value_options["part0_gfnv"],
        
        action="store",
        required=False,
        help="GFNn-xTB version employed for calculating the GFNn-xTB "
        "single point in part0. "
        f"Allowed values are [{', '.join(options.value_options['part0_gfnv'])}]",
    ) 
    group3.add_argument(
        "-part1",
        "--part1",
        choices=["on", "off"],
        dest="part1",
        action="store",
        required=False,
        
        help="Option to turn the prescreening evaluation (part1) 'on' or 'off'.",
    )
    group3.add_argument(
        "-smgsolv1",
        "--smgsolv1",
        choices=options.value_options["smgsolv1"],
        dest="smgsolv1",
        action="store",
        required=False,
        
        help="Solvent model for the Gsolv evaluation in part1. This can either be"
        " an implicit solvation or an additive solvation model. "
        f"Allowed values are [{', '.join(options.value_options['smgsolv1'])}]",
    ) 
    group10.add_argument(
        "-prescreening_threshold",
        "-prethr",
        "--thresholdpre",
        dest="prescreening_threshold",
        
        action="store",
        type=float,
        required=False,
        help=(
            "Threshold in kcal/mol. All conformers in part0 (prescreening)"
            " with a relativ energy below the threshold are considered for part1."
        ),
    ) 
    group4.add_argument(
        "-sm2",
        "--solventmodel2",
        choices=options.value_options.get("sm2"),
        dest="sm2",
        action="store",
        required=False,
        
        help="Solvent model employed during the geometry ensembleopt in part2."
        "The solvent model sm2 is not used for Gsolv evaluation, but for the "
        "implicit effect on a property (e.g. the geometry in the ensembleopt).",
    ) 
    group4.add_argument(
        "-smgsolv2",
        "--smgsolv2",
        choices=options.value_options["smgsolv2"],
        dest="smgsolv2",
        action="store",
        required=False,
        
        help="Solvent model for the Gsolv (solvation contribution to free energy) "
        "calculation in part2. Either the solvent"
        " model of the ensembleopt (sm2) or an additive solvation model. "
        f"Allowed values are [{', '.join(options.value_options['smgsolv2'])}]",
    ) """

    # TODO - keep this?
    """ group1.add_argument(
        "-prog_rrho",
        "--prog_rrho",
        choices=options.value_options["prog_rrho"],
        dest="prog_rrho",
        required=False,
        
        help="QM-program for mRRHO contribution in part1 2 and 3, currently only 'xtb'.",
    ) """

    # TODO - keep?
    """ group4.add_argument(
        "-ancopt",
        choices=["on"],  # there is no other option right now!
        dest="ancopt",
        required=False,
        
        help="Option to use xtb as driver for the xTB-optimizer in part2. "
        "Which is currently not changeable!",
    ) """

    args = parser.parse_args(argv)
    if not check_soft_requirements(args):
        raise argparse.ArgumentError(
            None,
            "You must provide an input file via '-i' and provide number of cores via '--maxcores'.",
        )

    return args
