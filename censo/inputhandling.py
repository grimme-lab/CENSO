"""
defininition of internal defaults, checking of logic for parameter combinations,
cml parsing
"""

import argparse
import os
from censo.cfg import (
    __version__,
)

# TODO - lineup args with settings_options
# TODO - arg for assets_path?
# removed ALL choices, since these are checked in InternalSettings.check_logic anyways
def cml(startup_description, argv=None):
    """
    Process commandline arguments
    """
    
    parser = argparse.ArgumentParser(
        description=startup_description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    groups = []

    ### RUN SETTINGS
    groups.append(parser.add_argument_group("RUN SETTINGS"))
    # TODO - add option to run only certain parts via cml?
    groups[0].add_argument(
        "-inp",
        "--input",
        dest="inp",
        type=os.path.abspath,
        action="store",
        required=False,
        metavar="",
        help="Input name of ensemble file: e.g. crest_conformers.xyz ",
    )
    groups[0].add_argument(
        "-nc",
        "--nconf",
        dest="nconf",
        type=int,
        action="store",
        required=False,
        metavar="",
        help="Number of conformers which are going to be considered.",
    )
    groups[0].add_argument(
        "-chrg",
        "--charge",
        dest="charge",
        type=int,
        action="store",
        required=False,
        metavar="",
        help="Charge of the investigated molecule.",
    )
    groups[0].add_argument(
        "-u",
        "--unpaired",
        dest="unpaired",
        type=int,
        action="store",
        required=False,
        metavar="",
        help="Integer number of unpaired electrons of the investigated molecule.",
    )
    groups[0].add_argument(
        "-checkinput",
        "--checkinput",
        dest="checkinput",
        action="store_true", # ???
        type=bool,
        required=False,
        help="Option to check if all necessary information for the CENSO "
        "calculation are provided and check if certain setting combinations "
        "make sence. Option to choose from : [???]", # TODO
    )
    groups[0].add_argument(
        "-version",
        "--version",
        dest="version",
        action="store_true",
        required=False,
        help="Print CENSO version and exit.",
    )
    groups[0].add_argument(
        "-consider_unconverged",
        "--consider_unconverged",
        dest="consider_unconverged",
        type=bool,
        metavar="",
        required=False,
        help=("Expert user option for including conformers, which were removed "
             "in a previous run ('stopped_before_converged'), "
             "in the (current) optimization. Choices: ???." # TODO
        ),
    )
    groups[0].add_argument(
        "--debug",
        "-debug",
        dest="debug",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    groups[0].add_argument(
        "--restart",
        "-restart",
        dest="restart",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    groups[0].add_argument(
        "--cleanup",
        "-cleanup",
        dest="cleanup",
        action="store_true",
        default=False,
        help="Delete unneeded files from current working directory.",
    )
    groups[0].add_argument(
        "--cleanup_all",
        "-cleanup_all",
        dest="cleanup_all",
        action="store_true",
        default=False,
        help="Delete all unneeded files from current working directory. "
        "Stronger than -cleanup !",
    )
    groups[0].add_argument(
        "--readonly",
        "-readonly",
        dest="onlyread",
        action="store",
        help="Create new enso.json from the output of previous calculations. "
        "This option has to used with exactly the same input settings or results will be unusable! "
        "",
    )
    groups[0].add_argument(
        "-newconfig",
        "-write_censorc",
        "--write_censorc",
        dest="writeconfig",
        default=False,
        action="store_true",
        required=False,
        help="Write new configuration file, which is placed into the current "
        "directory.",
    )
    groups[0].add_argument(
        "-copyinput",
        "--copyinput",
        dest="copyinput",
        default=False,
        action="store_true",
        required=False,
        help="Write all currently selected settings to a censo.inp configuration "
        "file, which is placed into the current directory.",
    )
    groups[0].add_argument(
        "-progress",
        "--progress",
        dest="progress",
        required=False,
        default="off",
        help="Print progress to stderr when starting and finishing a sorting Part."
        "Choices are 'on' or 'off'.",
    )
    groups[0].add_argument(
        "-inprc",
        "--inprc",
        dest="inprcpath",
        required=False,
        help="Path to the configuration file .censorc "
        "The default is ~/.censorc",
    )
    groups[0].add_argument(
        "-tutorial",
        "--tutorial",
        dest="tutorial",
        required=False,
        action="store_true",
        help="Start the interactive CENSO documentation.",
    )
    groups[0].add_argument(
        "-create_SI",
        "--create_SI",
        dest="create_si",
        required=False,
        action="store_true",
        help="Start automatic SI generation after CENSO run. (Work in progress)",
    )
    
    ### GENERAL SETTINGS
    groups.append(parser.add_argument_group("GENERAL SETTINGS"))
    groups[1].add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        type=float,
        action="store",
        required=False,
        metavar="",
        help="Temperature in Kelvin for thermostatistical evaluation.",
    )
    groups[1].add_argument(
        "-multitemp",
        "--multitemp",
        dest="multitemp",
        type=bool,
        required=False,
        metavar="",
        help="Needs to be turned on if a temperature range should be evaluated"
        " (flag trange). Options for multitemp are: [???].", # TODO
    )
    groups[1].add_argument(
        "-trange",
        "--trange",
        dest="trange",
        nargs=3,
        required=False,
        metavar=("start", "end", "step"),
        type=float,
        help="specify a temperature range [start, end, step] e.g.: 250.0 300.0 10.0"
        "  resulting in [250.0, 260.0, 270.0, 280.0, 290.0].",
    )
    groups[1].add_argument(
        "-bhess",
        "--bhess",
        dest="bhess",
        action="store",
        type=bool,
        required=False,
        metavar="",
        help="Uses SPH and applies structure constraint to input/DFT geometry "
        "for mRRHO calcuation. "
        "Options are: [???].", # TODO
    )
    groups[1].add_argument(
        "-consider_sym",
        "--consider_sym",
        dest="consider_sym",
        action="store",
        type=bool,
        required=False,
        metavar="",
        help="Consider symmetry in mRRHO calcuation (based on desy xtb threshold). "
        "Options are: [???].", # TODO
    )
    groups[1].add_argument(
        "-rmsdbias",
        "--rmsdbias",
        dest="rmsdbias",
        action="store",
        type=bool,
        required=False,
        metavar="",
        help="Applies constraint to rmsdpot.xyz to be consistent to CREST. "
        "Options are: [???].", # TODO
    )
    groups[1].add_argument(
        "-sm_rrho",
        "--sm_rrho",
        dest="sm_rrho",
        action="store",
        required=False,
        metavar="",
        help="Solvation model used in xTB GmRRHO calculation. Applied if not in "
        "gas-phase. Options are 'gbsa' or 'alpb'.",
    )
    groups[1].add_argument(
        "-evaluate_rrho",
        "--evaluate_rrho",
        dest="evaluate_rrho",
        action="store",
        type=bool,
        required=False,
        metavar="",
        help="Evaluate mRRHO contribution. Options: ???.", # TODO
    )
    groups[1].add_argument(
        "-solvent",
        "--solvent",
        dest="solvent",
        metavar="",
        action="store",
        required=False,
        help="Solvent the molecule is solvated in {}, available solvents "
        "are: {}. They can be extended in the "
        "file ~/.assets/censo_solvents.json .", # TODO
    )
    groups[1].add_argument(
        "-check",
        "--check",
        dest="check",
        action="store",
        type=bool,
        required=False,
        help="Option to terminate the CENSO-run if too many calculations/preparation"
        " steps fail. Options are: [???].", # TODO
    )
    groups[1].add_argument(
        "-cosmorsparam",
        "--cosmorsparam",
        dest="cosmorsparam",
        required=False,
        action="store",
        metavar="",
        help="Choose a COSMO-RS parametrization for possible COSMO-RS G_solv "
        "calculations: e.g. 19-normal for 'BP_TZVP_19.ctd' or 16-fine for"
        " 'BP_TZVPD_FINE_C30_1601.ctd'.",
    )
    groups[1].add_argument(
        "-O",
        "--omp",
        dest="omp",
        type=int,
        action="store",
        metavar="",
        help="Number of cores each thread can use. E.g. (maxthreads) 5 threads "
        "with each (omp) 4 cores --> 20 cores need to be available on the machine.",
    )
    groups[1].add_argument(
        "-P",
        "--maxthreads",
        dest="maxthreads",
        type=int,
        action="store",
        metavar="",
        help="Number of independent calculations during the CENSO calculation. E.g."
        " (maxthreads) 5 independent calculation-"
        " threads with each (omp) 4 cores --> 20 cores need to be available on "
        "the machine.",
    )
    groups[1].add_argument(
        "-imagthr",
        "--imagthr",
        dest="imagthr",
        type=float,
        action="store",
        metavar="",
        help="threshold for inverting imaginary frequencies for thermo in cm-1."
        " (e.g. -30.0)",
    )
    groups[1].add_argument(
        "-sthr",
        "--sthr",
        dest="sthr",
        action="store",
        type=float,
        metavar="",
        help="Rotor cut-off for thermo in cm-1. (e.g. 50.0)",
    )
    groups[1].add_argument(
        "-scale",
        "--scale",
        dest="scale",
        action="store",
        type=float,
        metavar="",
        help="scaling factor for frequencies  (e.g. 1.0)",
    )
    groups[1].add_argument(
        "--vapor_pressure",
        "-vp",
        dest="vapor_pressure",
        action="store",
        type=bool,
        help="Gsolv is evaluated for the input molecule in its solution (same). "
        "Only possible with COSMO-RS.",
    )
    
    ### PRESCREENING SETTINGS
    groups.append(parser.add_argument_group("PRESCREENING SETTINGS"))
    
    ### SCREENING SETTINGS
    groups.append(parser.add_argument_group("SCREENING SETTINGS"))
    
    ### OPTIMIZATION SETTINGS
    groups.append(parser.add_argument_group("OPTIMIZATION SETTINGS"))
    groups[4].add_argument(
        "-crestcheck",
        "--crestcheck",
        dest="crestcheck",
        action="store",
        type=bool,
        required=False,
        metavar="",
        help="Option to sort out conformers after DFT optimization which CREST "
        "identifies as identical or rotamers of each other. \nThe identification/"
        "analysis is always performed, but the removal of conformers has to "
        "be the choice of the user. Options are: [???]", # TODO
    )
    groups[4].add_argument(
        "-opt_spearman",
        dest="opt_spearman",
        type=bool,
        required=False,
        metavar="",
        help="Option to use an optimizer which checks if the hypersurface of DFT and"
        "xTB is parallel and optimizes mainly low lying conformers",
    )
    groups[4].add_argument(
        "-optlevel2",
        "--optlevel2",
        dest="optlevel2",
        default=None,
        required=False,
        metavar="",
        help="Option to set the optlevel in part2, only if optimizing with the xTB-optimizer!"
        "Allowed values are ***", # TODO
    )
    groups[4].add_argument(
        "-optcycles",
        "--optcycles",
        dest="optcycles",
        action="store",
        required=False,
        type=int,
        metavar="",
        help="number of cycles in ensemble optimizer.",
    )
    groups[4].add_argument(
        "-hlow",
        "--hlow",
        dest="hlow",
        action="store",
        required=False,
        type=float,
        metavar="",
        help="Lowest force constant in ANC generation (real), used by xTB-optimizer.",
    )
    groups[4].add_argument(
        "-spearmanthr",
        "--spearmanthr",
        dest="spearmanthr",
        action="store",
        required=False,
        metavar="",
        help="Value between -1 and 1 for the spearman correlation coeffient threshold, "
        "employed in the ensemlbe optimizer",
    )
    groups[4].add_argument(
        "-radsize",
        "--radsize",
        dest="radsize",
        action="store",
        required=False,
        metavar="",
        type=int,
        help="Radsize used in the optimization and only for r2scan-3c!",
    )
    
    ### REFINEMENT SETTINGS
    groups.append(parser.add_argument_group("REFINEMENT SETTINGS"))
    
    ### NMR SETTINGS
    groups.append(parser.add_argument_group("NMR SETTINGS"))
    groups[6].add_argument(
        "-couplings",
        "--couplings",
        dest="couplings",
        type=bool,
        required=False,
        metavar="",
        help="Option to run coupling constant calculations. Options are ???.",
    )
    groups[6].add_argument(
        "-shieldings",
        "--shieldings",
        dest="shieldings",
        type=bool,
        required=False,
        metavar="",
        help="Option to run shielding constant calculations. Options are ???.",
    )
    groups[6].add_argument(
        "-hactive",
        "--hactive",
        dest="h_active",
        type=bool,
        required=False,
        metavar="",
        help="Investigates hydrogen nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-cactive",
        "--cactive",
        dest="c_active",
        type=bool,
        required=False,
        metavar="",
        help="Investigates carbon nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-factive",
        "--factive",
        dest="f_active",
        type=bool,
        required=False,
        metavar="",
        help="Investigates fluorine nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-siactive",
        "--siactive",
        dest="si_active",
        type=bool,
        required=False,
        metavar="",
        help="Investigates silicon nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    groups[6].add_argument(
        "-pactive",
        "--pactive",
        dest="p_active",
        type=bool,
        required=False,
        metavar="",
        help="Investigates phosophorus nuclei in coupling and shielding calculations."
        "choices=[???]",
    )
    
    ### OPTROT SETTINGS
    groups.append(parser.add_argument_group("OPTROT SETTINGS"))
    groups[7].add_argument(
        "-freqOR",
        "--freqOR",
        dest="freq_or",
        nargs="*",
        required=False,
        type=float,
        metavar="",
        help="Frequencies to evaluate specific rotation at in nm, e.g. 589 "
        "or 589 700 to evaluate at 598 nm and 700 nm.",
    )
    
    ### UVVIS SETTINGS
    groups.append(parser.add_argument_group("UVVIS SETTINGS"))
    
    # TODO - split this into func for all different parts
    # leave these options out for now, implementation for cml complicated
    """ group1.add_argument(
        "-func",
        "--functional",
        dest="func",
        choices=options.value_options["func"],
        action="store",
        required=False,
        metavar="",
        help="Functional for geometry optimization (used in part2) and "
        "single-points in part1",
    ) 
    group1.add_argument(
        "-basis",
        "--basis",
        dest="basis",
        action="store",
        required=False,
        metavar="",
        help="Basis set employed together with the functional (func) for the "
        "low level single point in part1 und optimization in part2.",
    ) 
    group1.add_argument(
        "-prog",
        "--prog",
        choices=options.value_options["prog"],
        dest="prog",
        required=False,
        metavar="",
        help="QM-program used in part0, part1 and part2 either 'orca' or 'tm'.",
    ) 
    group10.add_argument(
        "-part0_gfnv",
        "--part0_gfnv",
        dest="part0_gfnv",
        choices=options.value_options["part0_gfnv"],
        metavar="",
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
        metavar="",
        help="Option to turn the prescreening evaluation (part1) 'on' or 'off'.",
    )
    group3.add_argument(
        "-smgsolv1",
        "--smgsolv1",
        choices=options.value_options["smgsolv1"],
        dest="smgsolv1",
        action="store",
        required=False,
        metavar="",
        help="Solvent model for the Gsolv evaluation in part1. This can either be"
        " an implicit solvation or an additive solvation model. "
        f"Allowed values are [{', '.join(options.value_options['smgsolv1'])}]",
    ) 
    group10.add_argument(
        "-prescreening_threshold",
        "-prethr",
        "--thresholdpre",
        dest="prescreening_threshold",
        metavar="",
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
        metavar="",
        help="Solvent model employed during the geometry optimization in part2."
        "The solvent model sm2 is not used for Gsolv evaluation, but for the "
        "implicit effect on a property (e.g. the geometry in the optimization).",
    ) 
    group4.add_argument(
        "-smgsolv2",
        "--smgsolv2",
        choices=options.value_options["smgsolv2"],
        dest="smgsolv2",
        action="store",
        required=False,
        metavar="",
        help="Solvent model for the Gsolv (solvation contribution to free energy) "
        "calculation in part2. Either the solvent"
        " model of the optimization (sm2) or an additive solvation model. "
        f"Allowed values are [{', '.join(options.value_options['smgsolv2'])}]",
    ) """
    
    # TODO - keep this?
    """ group1.add_argument(
        "-prog_rrho",
        "--prog_rrho",
        choices=options.value_options["prog_rrho"],
        dest="prog_rrho",
        required=False,
        metavar="",
        help="QM-program for mRRHO contribution in part1 2 and 3, currently only 'xtb'.",
    ) """
    
    
    # TODO - keep?
    """ group4.add_argument(
        "-ancopt",
        choices=["on"],  # there is no other option right now!
        dest="ancopt",
        required=False,
        metavar="",
        help="Option to use xtb as driver for the xTB-optimizer in part2. "
        "Which is currently not changeable!",
    ) """
    
    args = parser.parse_args(argv)

    return args