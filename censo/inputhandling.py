"""
defininition of internal defaults, checking of logic for parameter combinations,
cml parsing
"""

# TODO - ensure backwards-compatibility, increase maintainability

import argparse
import shutil
import os
import json
import csv
import time
import math
import sys
import glob
from collections import OrderedDict
from censo.core import CensoCore
from censo.cfg import (
    PLENGTH,
    DIGILEN,
    ENVIRON,
    CODING,
    WARNLEN,
    censo_solvent_db,
    external_paths,
    __version__,
    cosmors_param,
    dfa_settings,
    knownbasissets,
)
from censo.utilities import frange, format_line, print, print_block

# TODO - fix choices
# TODO - lineup args with settings_options
# TODO - arg for assets_path?
def cml(startup_description, argv=None):
    """
    Process commandline arguments
    """

    parser = argparse.ArgumentParser(
        description=startup_description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )
    group1 = parser.add_argument_group("GENERAL SETTINGS")
    group1.add_argument(
        "-inp",
        "--input",
        dest="inp",
        type=os.path.abspath,
        action="store",
        required=False,
        metavar="",
        help="Input name of ensemble file: e.g. crest_conformers.xyz ",
    )
    group1.add_argument(
        "-nc",
        "--nconf",
        dest="nconf",
        type=int,
        action="store",
        required=False,
        metavar="",
        help="Number of conformers which are going to be considered (max number "
        "of conformers are all conformers from the input file).",
    )
    group1.add_argument(
        "-chrg",
        "--charge",
        dest="charge",
        type=int,
        action="store",
        required=False,
        metavar="",
        help="Charge of the investigated molecule.",
    )
    group1.add_argument(
        "-u",
        "--unpaired",
        dest="unpaired",
        type=int,
        action="store",
        required=False,
        metavar="",
        help="Integer number of unpaired electrons of the investigated molecule.",
    )
    group1.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        type=float,
        action="store",
        required=False,
        metavar="",
        help="Temperature in Kelvin for thermostatistical evaluation.",
    )
    group1.add_argument(
        "-multitemp",
        "--multitemp",
        dest="multitemp",
        choices=["on", "off"],
        required=False,
        metavar="",
        help="Needs to be turned on if a temperature range should be evaluated"
        " (flag trange). Options for multitemp are: ['on' or 'off'].",
    )
    group1.add_argument(
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
    group1.add_argument(
        "-bhess",
        "--bhess",
        dest="bhess",
        choices=["on", "off"],
        action="store",
        required=False,
        metavar="",
        help="Uses SPH and applies structure constraint to input/DFT geometry "
        "for mRRHO calcuation. "
        "Options are: ['on' or 'off'].",
    )
    group1.add_argument(
        "-consider_sym",
        "--consider_sym",
        dest="consider_sym",
        choices=["on", "off"],
        action="store",
        required=False,
        metavar="",
        help="Consider symmetry in mRRHO calcuation (based on desy xtb threshold). "
        "Options are: ['on' or 'off'].",
    )
    group1.add_argument(
        "-rmsdbias",
        "--rmsdbias",
        dest="rmsdbias",
        choices=["on", "off"],
        action="store",
        required=False,
        metavar="",
        help="Applies constraint to rmsdpot.xyz to be consistent to CREST. "
        "Options are: ['on' or 'off'].",
    )
    group1.add_argument(
        "-sm_rrho",
        "--sm_rrho",
        dest="sm_rrho",
        choices=["gbsa", "alpb"],
        action="store",
        required=False,
        metavar="",
        help="Solvation model used in xTB GmRRHO calculation. Applied if not in "
        "gas-phase. Options are 'gbsa' or 'alpb'.",
    )
    group1.add_argument(
        "-evaluate_rrho",
        "--evaluate_rrho",
        dest="evaluate_rrho",
        action="store",
        choices=["on", "off"],
        required=False,
        metavar="",
        help="Evaluate mRRHO contribution. Options: on or off.",
    )
    group1.add_argument(
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
        "-checkinput",
        "--checkinput",
        dest="checkinput",
        action="store_true",
        required=False,
        help="Option to check if all necessary information for the CENSO "
        "calculation are provided and check if certain setting combinations "
        "make sence. Option to choose from : ['on' or 'off']",
    )
    group1.add_argument(
        "-solvent",
        "--solvent",
        dest="solvent",
        #choices=options.value_options["solvent"],
        metavar="",
        action="store",
        required=False,
        help="Solvent the molecule is solvated in {}, available solvents "
        "are: {}. They can be extended in the "
        "file ~/.assets/censo_solvents.json .".format(
            options.value_options["solvent"]
        ),
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
    group1.add_argument(
        "-prog_rrho",
        "--prog_rrho",
        choices=options.value_options["prog_rrho"],
        dest="prog_rrho",
        required=False,
        metavar="",
        help="QM-program for mRRHO contribution in part1 2 and 3, currently only 'xtb'.",
    )
    group1.add_argument(
        "-crestcheck",
        "--crestcheck",
        dest="crestcheck",
        choices=["on", "off"],
        action="store",
        required=False,
        metavar="",
        help="Option to sort out conformers after DFT optimization which CREST "
        "identifies as identical or rotamers of each other. \nThe identification/"
        "analysis is always performed, but the removal of conformers has to "
        "be the choice of the user. Options are: ['on' or 'off']",
    )
    group1.add_argument(
        "-check",
        "--check",
        dest="check",
        choices=["on", "off"],
        action="store",
        required=False,
        help="Option to terminate the CENSO-run if too many calculations/preparation"
        " steps fail. Options are: ['on' or 'off'].",
    )
    group1.add_argument(
        "-version",
        "--version",
        dest="version",
        action="store_true",
        required=False,
        help="Print CENSO version and exit.",
    )
    group1.add_argument(
        "-part3only",
        "--part3only",
        dest="part3only",
        required=False,
        action="store_true",
        help="Option to turn off part1 and part2",
    )
    group1.add_argument(
        "-cosmorsparam",
        "--cosmorsparam",
        dest="cosmorsparam",
        required=False,
        action="store",
        choices=options.value_options["cosmorsparam"],
        metavar="",
        help="Choose a COSMO-RS parametrization for possible COSMO-RS G_solv "
        "calculations: e.g. 19-normal for 'BP_TZVP_19.ctd' or 16-fine for"
        " 'BP_TZVPD_FINE_C30_1601.ctd'.",
    )
    group10 = parser.add_argument_group("CRE CHEAP-PRESCREENING - PART0")
    group10.add_argument(
        "-part0",
        "--part0",
        choices=["on", "off"],
        dest="part0",
        action="store",
        required=False,
        metavar="",
        help="Option to turn the CHEAP prescreening evaluation (part0) which "
        "improves description of ΔE compared to the input from semi-empirical "
        "methods 'on' or 'off'.",
    )
    group10.add_argument(
        "-func0",
        "--func0",
        dest="func0",
        choices=options.value_options["func0"],
        action="store",
        required=False,
        metavar="",
        help="Functional for fast single-point (used in part0)",
    )
    group10.add_argument(
        "-basis0",
        "--basis0",
        dest="basis0",
        action="store",
        required=False,
        metavar="",
        help="Basis set employed together with the functional (func0) for the "
        "fast single point calculation in part0.",
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
    group10.add_argument(
        "-part0_threshold",
        "-thrpart0",
        "--thresholdpart0",
        dest="part0_threshold",
        metavar="",
        action="store",
        required=False,
        help=(
            "Threshold in kcal/mol. All conformers in part0 (cheap single-point)"
            " with a relativ energy below the threshold are considered for part1."
        ),
    )

    group3 = parser.add_argument_group("CRE PRESCREENING - PART1")
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
    group3.add_argument(
        "-part1_gfnv",
        "--part1_gfnv",
        dest="part1_gfnv",
        choices=options.value_options["part1_gfnv"],
        metavar="",
        action="store",
        required=False,
        help="GFNn-xTB version employed for calculating the "
        "mRRHO contribution in part1. "
        f"Allowed values are [{', '.join(options.value_options['part1_gfnv'])}]",
    )
    group3.add_argument(
        "-part1_threshold",
        "-thrpart1",
        "--thresholdpart1",
        dest="part1_threshold",
        metavar="",
        action="store",
        required=False,
        help=(
            "Threshold in kcal/mol. All conformers in part1 (lax_single-point)"
            " with a relativ energy below the threshold are considered for part2."
        ),
    )

    group4 = parser.add_argument_group("CRE OPTIMIZATION - PART2")
    group4.add_argument(
        "-part2",
        "--part2",
        choices=["on", "off"],
        dest="part2",
        action="store",
        required=False,
        metavar="",
        help="Option to turn the geometry optimization (part2) 'on' or 'off'.",
    )
    group4.add_argument(
        "-prog2opt",
        "--prog2opt",
        dest="prog2opt",
        choices=options.value_options["prog"],
        action="store",
        required=False,
        metavar="",
        help="QM program used only for the optimization in part2. Either 'tm' or 'orca'.",
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
    )
    group4.add_argument(
        "-part2_gfnv",
        "--part2_gfnv",
        dest="part2_gfnv",
        choices=options.value_options["part2_gfnv"],
        metavar="",
        action="store",
        required=False,
        help="GFNn-xTB version employed for calculating the "
        "GmRRHO contribution in part2. "
        f"Allowed values are [{', '.join(options.value_options['part2_gfnv'])}]",
    )
    group4.add_argument(
        "-ancopt",
        choices=["on"],  # there is no other option right now!
        dest="ancopt",
        required=False,
        metavar="",
        help="Option to use xtb as driver for the xTB-optimizer in part2. "
        "Which is currently not changeable!",
    )
    group4.add_argument(
        "-opt_spearman",
        choices=["on", "off"],
        dest="opt_spearman",
        required=False,
        metavar="",
        help="Option to use an optimizer which checks if the hypersurface of DFT and"
        "xTB is parallel and optimizes mainly low lying conformers",
    )
    group4.add_argument(
        "-optlevel2",
        "--optlevel2",
        choices=options.value_options["optlevel2"],
        dest="optlevel2",
        default=None,
        required=False,
        metavar="",
        help="Option to set the optlevel in part2, only if optimizing with the xTB-optimizer!"
        "Allowed values are " + ", ".join(options.value_options["optlevel2"]),
    )
    group4.add_argument(
        "-optcycles",
        "--optcycles",
        dest="optcycles",
        action="store",
        required=False,
        type=int,
        metavar="",
        help="number of cycles in ensemble optimizer.",
    )
    group4.add_argument(
        "-hlow",
        "--hlow",
        dest="hlow",
        action="store",
        required=False,
        type=float,
        metavar="",
        help="Lowest force constant in ANC generation (real), used by xTB-optimizer.",
    )
    group4.add_argument(
        "-spearmanthr",
        "--spearmanthr",
        dest="spearmanthr",
        action="store",
        required=False,
        metavar="",
        help="Value between -1 and 1 for the spearman correlation coeffient threshold, "
        "employed in the ensemlbe optimizer",
    )
    group4.add_argument(
        #"-opt_limit",
        #"--opt_limit",
        "-thrpart2",
        "--thresholdpart2",
        "-part2_threshold",
        dest="opt_limit",
        action="store",
        required=False,
        metavar="",
        help=(
            "Lower limit Threshold in kcal/mol. If the GFNn and DFT hypersurfaces are"
            "assumed parallel, the conformers above the threshold are removed and"
            " not optimized further."
            "The conformers in part2 with a relativ free energy below the "
            "threshold are fully optimized."
        ),
    )
    group4.add_argument(
        #"-thrpart2",
        #"--thresholdpart2",
        #"-part2_threshold",
        "-part2_Pthreshold",
        "--part2_Pthreshold",
        dest="part2_P_threshold",
        action="store",
        required=False,
        metavar="",
        help=(
            "Boltzmann population sum threshold for part2 in %%. The conformers with "
            "the highest Boltzmann weigths are summed up until the threshold is reached."
            "E.g. all conformers up to a Boltzmann population of 90 %% are considered."
            'Example usage: "-thrpart2 99"  --> considers a population of 99 %%'
        ),
    )
    group4.add_argument(
        "-radsize",
        "--radsize",
        dest="radsize",
        action="store",
        required=False,
        metavar="",
        type=int,
        help=("Radsize used in the optimization and only for r2scan-3c!"),
    )
    group4.add_argument(
        "-consider_unconverged",
        "--consider_unconverged",
        choices=["on", "off"],
        dest="consider_unconverged",
        metavar="",
        required=False,
        help=("Expert user option for including conformers, which were removed "
             "in a previous run ('stopped_before_converged'), "
             "in the (current) optimization. Choices: 'on' or 'off'."
        ),
    )
    group5 = parser.add_argument_group("CRE REFINEMENT - PART3")
    group5.add_argument(
        "-part3",
        "--part3",
        choices=["on", "off"],
        dest="part3",
        action="store",
        required=False,
        metavar="",
        help="Option to turn the high level free energy evaluation (part3) 'on' or 'off'. "
        "It can be calculated on DFT optimized geometries (after part2) or directly on the "
        "ensemble input GFNn-xTB gemoetries by turning part0 part1 and part2 off.",
    )
    group5.add_argument(
        "-prog3",
        "--prog3",
        choices=options.value_options["prog3"],
        dest="prog3",
        required=False,
        metavar="",
        help="QM-program used in part3 either 'orca' or 'tm'.",
    )
    group5.add_argument(
        "-func3",
        "--functionalpart3",
        dest="func3",
        # choices=func3,
        action="store",
        required=False,
        metavar="",
        help="Functional for the high level single-point in part3. When Gsolv is "
        "calculated with COSMO-RS the energy and density for Gsolv is taken from "
        "func3/basis3 .",
    )
    group5.add_argument(
        "-basis3",
        "--basis3",
        dest="basis3",
        action="store",
        required=False,
        metavar="",
        help="Basis set employed together with the functional (func3) for the "
        "high level single point in part3.",
    )
    group5.add_argument(
        "-smgsolv3",
        "--smgsolv3",
        choices=options.value_options["smgsolv3"],
        dest="smgsolv3",
        action="store",
        required=False,
        metavar="",
        help="Solvent model for the Gsolv calculation in part3. Either the solvent"
        " model of the optimization (sm2) or an additive solvation model.",
    )
    group5.add_argument(
        "-part3_gfnv",
        "--part3_gfnv",
        dest="part3_gfnv",
        choices=options.value_options["part3_gfnv"],
        metavar="",
        action="store",
        required=False,
        help="GFNn-xTB version employed for calculating the "
        "GmRRHO contribution in part3. "
        f"Allowed values are [{', '.join(options.value_options['part3_gfnv'])}]",
    )
    group5.add_argument(
        "-thrpart3",
        "--thresholdpart3",
        dest="part3_threshold",
        action="store",
        required=False,
        metavar="",
        help=(
            "Boltzmann population sum threshold for part3 in %%. The conformers with "
            "the highest Boltzmann weigths are summed up until the threshold is reached."
            "E.g. all conformers up to a Boltzmann population of 90 %% are considered"
            'Example usage: "-thrpart3 99"  --> considers a population of 99 %%'
        ),
    )
    group6 = parser.add_argument_group("NMR Mode")
    group6.add_argument(
        "-part4",
        "--part4",
        choices=["on", "off"],
        dest="part4",
        action="store",
        required=False,
        metavar="",
        help="Option to turn the NMR property calculation mode (part4) 'on' or 'off'.",
    )
    group6.add_argument(
        "-couplings",
        "--couplings",
        dest="couplings",
        required=False,
        choices=["on", "off"],
        metavar="",
        help="Option to run coupling constant calculations. Options are 'on' or 'off'.",
    )
    group6.add_argument(
        "-prog4J",
        "--prog4J",
        # choices=options.value_options["prog"],
        dest="prog4_j",
        required=False,
        metavar="",
        help="QM-program for the calculation of coupling constants.",
    )
    group6.add_argument(
        "-funcJ",
        "--funcJ",
        dest="func_j",
        action="store",
        required=False,
        metavar="",
        help="Functional for the coupling constant calculation.",
    )
    group6.add_argument(
        "-basisJ",
        "--basisJ",
        dest="basis_j",
        action="store",
        required=False,
        metavar="",
        help="Basis set for the calculation of coupling constants.",
    )
    group6.add_argument(
        "-sm4_j",
        "--sm4_j",
        dest="sm4_j",
        action="store",
        required=False,
        metavar="",
        help="Solvation model used in the coupling constant calculation.",
    )
    group6.add_argument(
        "-shieldings",
        "--shieldings",
        dest="shieldings",
        required=False,
        choices=["on", "off"],
        metavar="",
        help="Option to run shielding constant calculations. Options are 'on' or 'off'.",
    )
    group6.add_argument(
        "-prog4S",
        "--prog4S",
        dest="prog4_s",
        required=False,
        metavar="",
        help="QM-program for the calculation of shielding constants.",
    )
    group6.add_argument(
        "-funcS",
        "--funcS",
        dest="func_s",
        action="store",
        required=False,
        metavar="",
        help="Functional for shielding constant calculation.",
    )
    group6.add_argument(
        "-basisS",
        "--basisS",
        dest="basis_s",
        action="store",
        required=False,
        metavar="",
        help="Basis set for the calculation of shielding constants.",
    )
    group6.add_argument(
        "-sm4_s",
        "--sm4_s",
        dest="sm4_s",
        action="store",
        required=False,
        metavar="",
        help="Solvation model used in the shielding constant calculation.",
    )
    group6.add_argument(
        "-hactive",
        "--hactive",
        choices=["on", "off"],
        dest="h_active",
        required=False,
        metavar="",
        help="Investigates hydrogen nuclei in coupling and shielding calculations."
        "choices=['on', 'off']",
    )
    group6.add_argument(
        "-cactive",
        "--cactive",
        choices=["on", "off"],
        dest="c_active",
        required=False,
        metavar="",
        help="Investigates carbon nuclei in coupling and shielding calculations."
        "choices=['on', 'off']",
    )
    group6.add_argument(
        "-factive",
        "--factive",
        choices=["on", "off"],
        dest="f_active",
        required=False,
        metavar="",
        help="Investigates fluorine nuclei in coupling and shielding calculations."
        "choices=['on', 'off']",
    )
    group6.add_argument(
        "-siactive",
        "--siactive",
        choices=["on", "off"],
        dest="si_active",
        required=False,
        metavar="",
        help="Investigates silicon nuclei in coupling and shielding calculations."
        "choices=['on', 'off']",
    )
    group6.add_argument(
        "-pactive",
        "--pactive",
        choices=["on", "off"],
        dest="p_active",
        required=False,
        metavar="",
        help="Investigates phosophorus nuclei in coupling and shielding calculations."
        "choices=['on', 'off']",
    )
    group9 = parser.add_argument_group("OPTICAL ROTATION MODE")
    group9.add_argument(
        "-OR",
        "--OR",
        "-part5",
        metavar="",
        choices=["on", "off"],
        action="store",
        dest="optical_rotation",
        required=False,
        help="Do optical rotation calculation.",
    )
    group9.add_argument(
        "-funcOR",
        "--funcOR",
        dest="func_or",
        # choices=func_or,
        action="store",
        required=False,
        metavar="",
        help="Functional for optical rotation calculation.",
    )
    group9.add_argument(
        "-funcOR_SCF",
        "--funcOR_SCF",
        dest="func_or_scf",
        # choices=func_or,
        action="store",
        required=False,
        metavar="",
        help="Functional used in SCF for optical rotation calculation.",
    )
    group9.add_argument(
        "-basisOR",
        "--basisOR",
        dest="basis_or",
        # choices=func_or,
        action="store",
        required=False,
        metavar="",
        help="Basis set for optical rotation calculation.",
    )
    group9.add_argument(
        "-freqOR",
        "--freqOR",
        dest="freq_or",
        nargs="*",
        required=False,
        type=float,
        metavar="",
        help="Frequencies to evaluate specific rotation at in nm. E.g. 589 "
        "Or 589 700 to evaluate at 598 nm and 700 nm.",
    )

    group7 = parser.add_argument_group("OPTIONS FOR PARALLEL CALCULATIONS")
    group7.add_argument(
        "-O",
        "--omp",
        dest="omp",
        type=int,
        action="store",
        metavar="",
        help="Number of cores each thread can use. E.g. (maxthreads) 5 threads "
        "with each (omp) 4 cores --> 20 cores need to be available on the machine.",
    )
    group7.add_argument(
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
    group7.add_argument(
        "-balance",
        "--balance",
        dest="balance",
        choices=["on", "off"],
        action="store",
        metavar="",
        help="Automatically balance the number of threads and cores when a low number"
        "of conformers is left. (never exceed O*P cores). Options are: 'on' or 'off'.",
    )
    group11 = parser.add_argument_group("Concerning overall GmRRHO calculations")
    group11.add_argument(
        "-imagthr",
        "--imagthr",
        dest="imagthr",
        action="store",
        metavar="",
        help="threshold for inverting imaginary frequencies for thermo in cm-1."
        " (e.g. -30.0)",
    )
    group11.add_argument(
        "-sthr",
        "--sthr",
        dest="sthr",
        action="store",
        metavar="",
        help="Rotor cut-off for thermo in cm-1. (e.g. 50.0)",
    )
    group11.add_argument(
        "-scale",
        "--scale",
        dest="scale",
        action="store",
        metavar="",
        help="scaling factor for frequencies  (e.g. 1.0)",
    )
    group12 = parser.add_argument_group("SPECIAL PURPOSE")
    group12.add_argument(
        "--vapor_pressure",
        "-vp",
        dest="vapor_pressure",
        choices=["on", "off"],
        action="store",
        help="Gsolv is evaluated for the input molecule in its solution (same). "
        "Only possible with COSMO-RS.",
    )

    group8 = parser.add_argument_group("CREATION/DELETION OF FILES")
    group8.add_argument(
        "--debug",
        "-debug",
        dest="debug",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    group8.add_argument(
        "--restart",
        "-restart",
        dest="restart",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    group8.add_argument(
        "--cleanup",
        "-cleanup",
        dest="cleanup",
        action="store_true",
        default=False,
        help="Delete unneeded files from current working directory.",
    )
    group8.add_argument(
        "--cleanup_all",
        "-cleanup_all",
        dest="cleanup_all",
        action="store_true",
        default=False,
        help="Delete all unneeded files from current working directory. "
        "Stronger than -cleanup !",
    )
    group8.add_argument(
        "--readonly",
        "-readonly",
        dest="onlyread",
        choices=["on", "off"],
        action="store",
        help="Create new enso.json from the output of previous calculations. "
        "This option has to used with exactly the same input settings or results will be unusable! "
        "",
    )
    group8.add_argument(
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
    group8.add_argument(
        "-copyinput",
        "--copyinput",
        dest="copyinput",
        default=False,
        action="store_true",
        required=False,
        help="Write all currently selected settings to a censo.inp configuration "
        "file, which is placed into the current directory.",
    )
    group8.add_argument(
        "-progress",
        "--progress",
        dest="progress",
        required=False,
        choices=["on", "off"],
        default="off",
        help="Print progress to stderr when starting and finishing a sorting Part."
        "Choices are 'on' or 'off'.",
    )
    group8.add_argument(
        "-inprc",
        "--inprc",
        dest="inprcpath",
        required=False,
        help="Path to the destination of the configuration file .censorc "
        "The default is ~/.censorc",
    )
    group8.add_argument(
        "-tutorial",
        "--tutorial",
        dest="tutorial",
        required=False,
        action="store_true",
        help="Start the interactive CENSO documentation.",
    )
    group8.add_argument(
        "-create_SI",
        "--create_SI",
        dest="create_si",
        required=False,
        action="store_true",
        help="Start automatic SI generation after CENSO run. (Work in progress)",
    )
    # TODO - add args for uvvis
    args = parser.parse_args(argv)

    if args.part3only:
        setattr(args, "part0", "off")
        setattr(args, "part1", "off")
        setattr(args, "part2", "off")
    return args


class config_setup(internal_settings):
    """
    Read or write configuration or input files.
    """

    def __init__(self, path=os.getcwd(), *args, **kwargs):
        internal_settings.__init__(self, *args, **kwargs)
        # settings just to calm down pylint, real assignment is dynamically done
        # general settings
        self.nconf = None
        self.charge = 0
        self.unpaired = 0
        self.solvent = "gas"
        self.prog_rrho = "xtb"
        self.temperature = 298.15
        self.trange = [273.15, 378.15, 5]
        self.multitemp = False
        self.evaluate_rrho = True
        self.bhess = True
        self.consider_sym = False
        self.imagthr = -50
        self.sthr = 50
        self.scale = 1.0
        self.sm_rrho = "alpb"
        self.check = True
        self.crestcheck = False
        self.prog = "tm"
        self.func = "b97-3c"
        self.basis = "automatic"
        self.maxthreads = 1
        self.omp = 1
        self.balance = False
        self.cosmorsparam = "automatic"

        # part0
        self.part0 = False
        self.part0_gfnv = "gfnff"
        self.part0_threshold = 4.0
        self.func0 = "b97-d"
        self.basis0 = "def2-SV(P)"
        # part1
        self.part1 = True
        self.smgsolv1 = "sm2"
        self.part1_gfnv = "gfnff"
        self.part1_threshold = 1.0
        # part2
        self.part2 = True
        self.part2_P_threshold = 90
        self.sm2 = "default"
        self.smgsolv2 = "sm2"
        self.part2_gfnv = "gfnff"
        self.ancopt = True
        self.hlow = 0.01
        self.opt_spearman = False
        self.optcycles = 5
        self.optlevel2 = "automatic"
        self.spearmanthr = 0.9999
        self.radsize = 8
        # part3
        self.part3 = True
        self.prog3 = "prog"
        self.func3 = "b97-d"
        self.basis3 = "def2-TZVPD"
        self.smgsolv3 = "sm2"
        self.part3_gfnv = "gfn2"
        # part4
        self.part4 = False
        self.prog4_j = "tm"
        self.prog4_s = "tm"
        self.couplings = True
        self.func_j = "pbe0"
        self.basis_j = "def2-TZVP"
        self.shieldings = True
        self.func_s = "pbe0"
        self.basis_s = "def2-TZVP"
        self.sm4_j = "default"
        self.sm4_s = "default"
        self.h_ref = "TMS"
        self.c_ref = "TMS"
        self.f_ref = "CFCl3"
        self.si_ref = "TMS"
        self.p_ref = "TMP"
        self.resonance_frequency = 300.0
        # part5
        self.optical_rotation = False
        self.func_or = "pbe"
        self.func_or_scf = "r2scan-3c"
        self.basis_or = "def2-SVPD"
        self.freq_or = [589]
        # special
        self.vapor_pressure = False

        # settings the program operates with updated to the defaults
        for key in self.internal_defaults.keys():
            setattr(self, key, self.internal_defaults[key]["default"])       

        # workingdirectory
        self.cwd = path
        # path to .xyz file (input)
        self.ensemblepath = ""
        self.configpath = ""
        self.jsonpath = ""

        # formatting:
        self.lenconfx = 3

        self.save_errors = []
        self.save_infos = []



        
    def _set_fixed_temperature(self):
        """Initialize the temperature employed in part0 part1 and crude_rrho in the
        optimization (part2). This temperature is kept fixed for all restarts, to be able
        to compare the same (free) energies."""
        if self.fixed_temperature is None:
            setattr(self, "fixed_temperature", float(getattr(self, "temperature")))

    def cleanup_run(self, complete):
        """
        Delete all unneeded files.
        """
        files_in_cwd = [
            f for f in os.listdir(self.cwd) if os.path.isfile(os.path.join(self.cwd, f))
        ]
        for file in files_in_cwd:
            if (
                "enso.json." in file
                or "enso_ensemble_part1.xyz." in file
                or "enso_ensemble_part2.xyz." in file
                or "enso_ensemble_part3.xyz." in file
            ):
                if int(file.split(".")[2]) > 1:
                    print(f"Removing: {file}")
                    os.remove(os.path.join(self.cwd, file))
        # get files like amat.tmp from COSMO calculation (can be several Mb)
        files_in_sub_cwd = glob.glob(self.cwd + "/**/**/**/*mat.tmp", recursive=True)
        size = 0
        for tmpfile in files_in_sub_cwd:
            if any(x in tmpfile for x in ["a3mat.tmp", "a2mat.tmp", "amat.tmp"]):
                if os.path.isfile(tmpfile):
                    size += os.path.getsize(tmpfile)
                    # print(f"Removing {tmpfile} {os.path.getsize(tmpfile)} byte")
                    os.remove(tmpfile)
        print(f"Removed {size/(1024*1024): .2f} Mb")
        if complete:
            if "enso.json" in files_in_cwd:
                print(f"Removing: {'enso.json'}")
                os.remove(os.path.join(self.cwd, "enso.json"))
            if "enso.json.1" in files_in_cwd:
                print(f"Removing: {'enso.json.1'}")
                os.remove(os.path.join(self.cwd, "enso.json.1"))
            if os.path.isdir(os.path.join(self.cwd, "conformer_rotamer_check")):
                print("Removing conformer_rotamer_check")
                shutil.rmtree(os.path.join(self.cwd, "conformer_rotamer_check"))
            # ask if CONF folders should be removed

    def get_method_name(
        self,
        jobtype,
        func=None,
        basis=None,
        sm=None,
        gfn_version=None,
        bhess=None,
        solvent=None,
        prog=None,
        func2=None,
        disp=None,
    ):
        """
        Create method name for storing and retrieving data
        --> method energy
        --> method2 gsolv
        """
        if func is not None and basis is not None:
            if func in dfa_settings.composite_method_basis.keys():
                if basis == dfa_settings.composite_method_basis[func]:
                    # composite method (e.g. r2scan-3c)
                    tmp_func_basis = func
                else:
                    # FUNC/BASIS
                    tmp_func_basis = f"{func}/{basis}"
            elif disp is not None:
                # FUNC-DISP/BASIS
                tmp_func_basis = f"{func}-{disp}/{basis}"
            else:
                # FUNC/BASIS
                tmp_func_basis = f"{func}/{basis}"
        if jobtype in ("cosmors",):
            exc_name = {"cosmors": "COSMO-RS-normal", "cosmors-fine": "COSMO-RS-fine"}
            # energy  FUNC/BASIS
            method = tmp_func_basis
            # cosmors gsolv  COSMO-RS[FUNC/BASIS]
            method2 = f"{exc_name.get(sm)}[{tmp_func_basis}]"
        elif jobtype in ("gbsa_gsolv", "alpb_gsolv"):
            # energy  FUNC/BASIS
            method = tmp_func_basis
            # e.g. ALPB_Gsolv[GFN2]
            method2 = f"{sm}[{gfn_version}]"
        elif jobtype == "sp":
            # energy  FUNC/BASIS
            method = tmp_func_basis
        elif jobtype == "sp_implicit":
            # energy  FUNC/BASIS[DCOSMORS]
            method = f"{tmp_func_basis}[{str(sm).upper()}]"
            method2 = "incl. in E"
        elif jobtype == "smd_gsolv":
            # energy  FUNC/BASIS
            method = tmp_func_basis
            # SMD_gsolv  SMD_GSOLV[FUNC/BASIS]
            method2 = f"{sm}[{tmp_func_basis}]"
        elif jobtype == "rrhoxtb":
            # GFN2-bhess
            if bhess:
                if solvent != "gas":
                    method = f"{str(gfn_version).upper()}[{sm}]-bhess"
                else:
                    method = f"{str(gfn_version).upper()}-bhess"
            else:
                if solvent != "gas":
                    method = f"{str(gfn_version).upper()}[{sm}]"
                else:
                    method = f"{str(gfn_version).upper()}"
        elif jobtype in ("opt", "xtbopt"):
            if solvent == "gas":
                # energy  FUNC/BASIS
                method = tmp_func_basis
            else:
                # energy  FUNC/BASIS[DCOSMORS]
                method = f"{tmp_func_basis}[{str(sm).upper()}]"
        elif jobtype in ("couplings", "couplings_sp", "shieldings", "shieldings_sp"):
            if solvent == "gas":
                method = f"{tmp_func_basis}-{prog}"
            else:
                method = f"{tmp_func_basis}[{str(sm).upper()}]-{prog}"
        elif jobtype in ("opt-rot", "opt-rot_sp"):
            if func2 == "r2scan-3c" and basis != "def2-mTZVPP":
                func2 = f"r2scan/{basis}"
            if solvent == "gas":
                method = f"{tmp_func_basis}_[SCF={func2}]({prog})"
            else:
                method = f"{tmp_func_basis}[{str(sm).upper()}]_[SCF={func2}]({prog})"
        else:
            raise Exception(f"JOBTYPE {jobtype} not known in get_method_name")
        try:
            method2
        except NameError:
            method2 = ""
        return method, method2

    def provide_runinfo(self, extend=True):
        """
        Write dictionary structured like internal defaults.
        And extended with startup information.
        """
        runinfo = []
        keys = list(self.internal_defaults.keys())
        if extend:
            keys = keys + self.startupinfokeys
        for key in keys:
            runinfo.append((key, getattr(self, key)))
        return OrderedDict(runinfo)

    def _decomment(self, csvfile):
        """
        remove any comments from file before parsing with csv.DictReader
        comment symbols are # and $
        """
        for row in csvfile:
            raw = row.split("#")[0].strip()
            raw2 = raw.split("$")[0].strip()
            if raw2:
                yield raw2

    def _exchange_onoff(self, inp, reverse=False):
        """
        Exchange on --> True, off--> False, backward if reverse=True
        """
        exchange = {"on": True, "off": False}
        if reverse:
            if isinstance(inp, bool) and inp in {v: k for k, v in exchange.items()}:
                return {v: k for k, v in exchange.items()}[inp]
            else:
                return inp
        elif not reverse:
            if isinstance(inp, str) and inp in exchange.keys():
                return exchange[inp]
            else:
                return inp

    def read_config(self, path, startread, args):
        """
        Read from config data from file (here enso.inp or .censorc),
        cml > .censorc > internal defaults
        Creates settings for current run from censorc or similar
        => should be the primary way of inputting settings!!
        """
        # opens censorc with csvreader
        # iterates over k/v pairs given in censorc and adjusts config_setup attributes accordingly
        rcdata = {}
        with open(path, "r") as csvfile:
            # skip header:
            while True:
                line = csvfile.readline()
                if line.startswith(startread):
                    break
                elif line == "":
                    # EOF
                    break
                else:
                    pass
            reader = csv.DictReader(
                self._decomment(csvfile),
                fieldnames=("key", "value"),
                skipinitialspace=True,
                delimiter=":",
            )
            for row in reader:
                if "end" in row["key"]:
                    break
                else:
                    rcdata[row["key"]] = row["value"]
        if "end" in rcdata:
            del rcdata["end"]

        # read args from cml input and update censorc
        # FIXME - why change the rcfile and not just override settings from censorc temporarily?
        # if this stays like that this makes the C in CENSO pretty useless
        args_key = {v: k for k, v in self.key_args_dict.items()}
        cmlflags = vars(args)
        for key in cmlflags.keys():
            if key in args_key.keys():
                if cmlflags[key] is not None:
                    # print(f"SETTING cml: {key} to {cmlflags[key]}")
                    rcdata[key] = cmlflags[key]
                    # print(key, cmlflags[key])
        # end get commandline arguments
        # update censorc-key to internal key
        for key, value in list(rcdata.items()):
            if key in self.key_args_dict.keys():
                if key != self.key_args_dict[key]:
                    # print(f"updating: {key} to {self.key_args_dict[key]} "
                    #      "{value} {rcdata.get(self.key_args_dict[key])}")
                    rcdata[self.key_args_dict[key]] = rcdata.get(
                        self.key_args_dict[key], value
                    )
                    del rcdata[key]
        # end update censorc-key to internal key

        readinkeys = []
        for item in list(rcdata.keys()):
            if item not in self.internal_defaults.keys():
                self.save_errors.append(
                    f"{'WARNING:':{WARNLEN}}{item} is not a known "
                    f"keyword in {os.path.basename(path)}."
                )
                del rcdata[item]
            else:
                readinkeys.append(item)
        diff = list(set(self.internal_defaults.keys()) - set(readinkeys))
        if diff:
            self.save_errors.append(
                f"{'WARNING:':{WARNLEN}}These keywords were not found in the configuration "
                f"file {os.path.basename(path)}\n{'':{WARNLEN}}and therefore default "
                f"values are taken for:"
            )
            for item in diff:
                self.save_errors.append(f"{'':{WARNLEN}}{item}")
                rcdata[item] = self.internal_defaults[item]["default"]

        for key in rcdata:
            if rcdata[key] == "":
                rcdata[key] = None
        # -----> keys are checked, now check values!!!!
        for key in rcdata:
            # change on --> True , off --> False
            rcdata[key] = self._exchange_onoff(rcdata[key])
            if key != "nconf":
                if (
                    not isinstance(rcdata[key], self.internal_defaults[key]["type"])
                    and rcdata[key] is not None
                ):
                    try:
                        if self.internal_defaults[key]["type"] == list:
                            tmp = rcdata[key].strip("[")
                            tmp = tmp.strip("]")
                            tmp = tmp.split(",")
                            rcdata[key] = [float(i) for i in tmp]
                        else:
                            rcdata[key] = self.internal_defaults[key]["type"](
                                rcdata[key]
                            )
                    except (ValueError, TypeError):
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{key}= {rcdata[key]}"
                            " could not be"
                            " converted and default values are set to "
                            f"{self.internal_defaults[key]['default']}"
                        )
                        rcdata[key] = self.internal_defaults[key]["type"](
                            self.internal_defaults[key]["default"]
                        )
            if key == "nconf":
                if rcdata[key] == "all":
                    rcdata[key] = None
                else:
                    try:
                        rcdata[key] = int(rcdata[key])
                    except ValueError:
                        rcdata[key] = None
        for key, value in rcdata.items():
            if key in vars(self).keys():
                setattr(self, key, value)
            else:
                print(f"{'ERROR:':{WARNLEN}}", key)
                self.save_errors.append(f"{key} not known in config!")


    def print_parameters(self, cmlcall=None):
        """
        print settings at startup
        """

        # print parameter setting
        print("\n" + "".ljust(PLENGTH, "-"))
        print("PARAMETERS".center(PLENGTH, " "))
        print("".ljust(PLENGTH, "-") + "\n")

        if cmlcall is not None:
            print(f"program call: {' '.join(cmlcall)}")
        print(
            f"The configuration file {os.path.basename(self.configpath)} is read "
            f"from {self.configpath}."
        )
        print(f"Reading conformer rotamer ensemble from: {self.ensemblepath}.")
        if self.save_infos:
            for _ in list(self.save_infos):
                print(self.save_infos.pop(0))
        info = []
        info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
        info.append(["justprint", "CRE SORTING SETTINGS".center(int(PLENGTH / 2), " ")])
        info.append(["justprint", "".ljust(int(PLENGTH / 2), "-") + "\n"])
        info.append(["nat", "number of atoms in system"])
        info.append(["nconf", "number of considered conformers"])
        info.append(["maxconf", "number of all conformers from input"])
        info.append(["charge", "charge"])
        info.append(["unpaired", "unpaired"])
        info.append(["solvent", "solvent"])
        info.append(["temperature", "temperature"])
        if self.multitemp:
            info.append(["multitemp", "evaluate at different temperatures"])
            info.append(
                [
                    "printoption",
                    "temperature range",
                    [i for i in frange(self.trange[0], self.trange[1], self.trange[2])],
                ]
            )
        info.append(["evaluate_rrho", "calculate mRRHO contribution"])
        info.append(["consider_sym", "consider symmetry for mRRHO contribution"])
        info.append(["check", "cautious checking for error and failed calculations"])
        info.append(["crestcheck", "checking the DFT-ensemble using CREST"])
        info.append(["maxthreads", "maxthreads"])
        info.append(["omp", "omp"])
        info.append(["balance", "automatically balance maxthreads and omp"])
        if self.onlyread in ("on", True):
            info.append(
                [
                    "printoption",
                    "onlyread",
                    {True: "on"}.get(self.onlyread, self.onlyread),
                ]
            )
        if self.vapor_pressure in ("on", True):
            info.append(
                [
                    "printoption",
                    "vapor_pressure",
                    {True: "on"}.get(self.vapor_pressure, self.vapor_pressure),
                ]
            )

        if self.part0:
            # PART0:
            info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
            info.append(
                [
                    "justprint",
                    "CRE CHEAP-PRESCREENING - PART0".center(int(PLENGTH / 2), " "),
                ]
            )
            info.append(["justprint", "".ljust(int(PLENGTH / 2), "-")])
            info.append(["part0", "part0"])
            info.append(["nconf", "starting number of considered conformers"])
            info.append(["prog", "program for part0"])
            info.append(["func0", "functional for fast single-point"])
            info.append(["basis0", "basis set for fast single-point"])
            info.append(["part0_threshold", "threshold g_thr(0) for sorting in part0"])
            if self.solvent != "gas":
                info.append(["sm_rrho", "Solvent model used with xTB"])

            tmp_func_basis, _ = self.get_method_name(
                "sp", func=getattr(self, "func0"), basis=getattr(self, "basis0")
            )
            info.append(
                [
                    "justprint",
                    f"\nshort-notation:\n{tmp_func_basis} "
                    "// GFNn-xTB (Input geometry)",
                ]
            )
        if self.part1:
            # PART1:
            info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
            info.append(
                ["justprint", "CRE PRESCREENING - PART1".center(int(PLENGTH / 2), " ")]
            )
            info.append(["justprint", "".ljust(int(PLENGTH / 2), "-")])
            info.append(["part1", "part1"])
            if not self.part0:
                info.append(["nconf", "starting number of considered conformers"])
            info.append(["prog", "program for part1"])
            info.append(["func", "functional for initial evaluation"])
            info.append(["basis", "basis set for initial evaluation"])
            info.append(["evaluate_rrho", "calculate mRRHO contribution"])
            if self.evaluate_rrho:
                info.append(["prog_rrho", "program for mRRHO contribution"])
                if self.prog_rrho == "xtb" or self.smgsolv2 in (
                    "gbsa_gsolv",
                    "alpb_gsolv",
                ):
                    info.append(
                        ["part1_gfnv", "GFN version for mRRHO and/or GBSA_Gsolv"]
                    )
                    info.append(
                        [
                            "bhess",
                            "Apply constraint to input geometry during mRRHO calculation",
                        ]
                    )
                    if self.solvent != "gas":
                        info.append(["sm_rrho", "solvent model applied with xTB"])
            info.append(["printoption", "evaluate at different temperatures", "off"])
            info.append(
                [
                    "part1_threshold",
                    "threshold g_thr(1) and G_thr(1) for sorting in part1",
                ]
            )
            if self.solvent != "gas":
                info.append(
                    ["smgsolv1", "solvent model for Gsolv contribution of part1"]
                )
            # shortnotation:
            tmp_rrho_method, _ = self.get_method_name(
                "rrhoxtb",
                bhess=self.bhess,
                gfn_version=self.part1_gfnv,
                sm=self.sm_rrho,
                solvent=self.solvent,
            )
            tmp_func_basis, _ = self.get_method_name(
                "sp", func=getattr(self, "func"), basis=getattr(self, "basis")
            )
            if self.solvent != "gas":
                info.append(
                    [
                        "justprint",
                        f"\nshort-notation:\n{tmp_func_basis} + "
                        f"{str(getattr(self, 'smgsolv1')).upper()}[{self.solvent}] "
                        f"+ GmRRHO({tmp_rrho_method}) "
                        f"// GFNn-xTB (Input geometry)",
                    ]
                )
            else:
                info.append(
                    [
                        "justprint",
                        f"\nshort-notation:\n{tmp_func_basis} "
                        f"+ GmRRHO({str(getattr(self, 'part1_gfnv')).upper()}) "
                        "// GFNn-xTB (Input geometry)",
                    ]
                )
        if self.part2:
            # PART2:
            info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
            info.append(
                ["justprint", "CRE OPTIMIZATION - PART2".center(int(PLENGTH / 2), " ")]
            )
            info.append(["justprint", "".ljust(int(PLENGTH / 2), "-")])
            info.append(["part2", "part2"])
            info.append(["prog", "program"])
            info.append(["prog2opt", "QM program employed only in the geometry opt."])
            info.append(["func", "functional for part2"])
            info.append(["basis", "basis set for part2"])
            info.append(["ancopt", "using xTB-optimizer for optimization"])
            if self.opt_spearman:
                info.append(["opt_spearman", "using the new ensemble optimizer"])
                info.append(
                    [
                        "opt_limit",
                        "optimize all conformers below this G_thr(opt,2) threshold",
                    ]
                )
                info.append(["printoption", "spearmanthr", f"{self.spearmanthr:.3f}"])
            if self.ancopt and self.optlevel2 is not None:
                info.append(["optlevel2", "optimization level in part2"])
            if self.solvent != "gas":
                info.append(["sm2", "solvent model applied in the optimization"])
                if self.smgsolv2 not in (None, "sm2"):
                    info.append(["smgsolv2", "solvent model for Gsolv contribution"])
            info.append(["multitemp", "evaluate at different temperatures"])
            info.append(
                [
                    "part2_P_threshold",
                    "Boltzmann sum threshold G_thr(2) for sorting in part2",
                ]
            )
            info.append(["evaluate_rrho", "calculate mRRHO contribution"])
            if self.evaluate_rrho:
                info.append(["prog_rrho", "program for mRRHO contribution"])
                if self.prog_rrho == "xtb":
                    info.append(
                        ["part2_gfnv", "GFN version for mRRHO and/or GBSA_Gsolv"]
                    )
                    if self.bhess:
                        info.append(
                            [
                                "bhess",
                                "Apply constraint to input geometry "
                                + "during mRRHO calculation",
                            ]
                        )
                    if self.solvent != "gas":
                        info.append(["sm_rrho", "solvent model applied with xTB"])
            if self.consider_unconverged:
                info.append(["consider_unconverged", 
                "consider (unconverged) conformers from previous run"])            
            # shortnotation:
            tmp_rrho_method, _ = self.get_method_name(
                "rrhoxtb",
                bhess=self.bhess,
                gfn_version=self.part2_gfnv,
                sm=self.sm_rrho,
                solvent=self.solvent,
            )
            tmp_func_basis, _ = self.get_method_name(
                "sp", func=getattr(self, "func"), basis=getattr(self, "basis")
            )
            if self.solvent != "gas":
                info.append(
                    [
                        "justprint",
                        f"\nshort-notation:\n{tmp_func_basis} + "
                        f"{str(getattr(self, 'smgsolv2')).upper()}[{self.solvent}] "
                        f"+ GmRRHO({tmp_rrho_method}) // "
                        f"{tmp_func_basis}"
                        f"[{str(getattr(self, 'sm2')).upper()}] ",
                    ]
                )
            else:
                info.append(
                    [
                        "justprint",
                        f"\nshort-notation:\n{tmp_func_basis} "
                        f"+ GmRRHO({str(getattr(self, 'part2_gfnv')).upper()}) "
                        f"// {tmp_func_basis}",
                    ]
                )
        # PART3:
        if self.part3:
            info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
            info.append(
                ["justprint", "CRE REFINEMENT - PART3".center(int(PLENGTH / 2), " ")]
            )
            info.append(["justprint", "".ljust(int(PLENGTH / 2), "-")])
            info.append(["part3", "part3"])
            info.append(["part3_threshold", "Boltzmann sum threshold G_thr(3)"])
            info.append(["prog3", "program for part3"])
            info.append(["func3", "functional for part3"])
            info.append(["basis3", "basis set for part3"])
            if self.solvent != "gas":
                info.append(["smgsolv3", "solvent model"])
            info.append(["multitemp", "evaluate at different temperatures"])
            info.append(["prog_rrho", "program for mRRHO contribution"])
            if self.prog_rrho == "xtb":
                info.append(["part3_gfnv", "GFN version for mRRHO and/or GBSA_Gsolv"])
                if self.bhess:
                    info.append(
                        [
                            "bhess",
                            "Apply constraint to input geometry during mRRHO calculation",
                        ]
                    )
                if self.solvent != "gas":
                    info.append(["sm_rrho", "solvent model applied with xTB"])
            # shortnotation:
            tmp_rrho_method, _ = self.get_method_name(
                "rrhoxtb",
                bhess=self.bhess,
                gfn_version=self.part3_gfnv,
                sm=self.sm_rrho,
                solvent=self.solvent,
            )
            tmp_func3_basis3, _ = self.get_method_name(
                "sp", func=getattr(self, "func3"), basis=getattr(self, "basis3")
            )
            tmp_func_basis, _ = self.get_method_name(
                "sp", func=getattr(self, "func"), basis=getattr(self, "basis")
            )
            if self.part2:
                if self.solvent != "gas":
                    geometries = (
                        f"{tmp_func_basis}[{str(getattr(self, 'sm2')).upper()}] "
                    )
                else:
                    geometries = f"{tmp_func_basis}"
            else:
                # on unoptimized geometries:
                geometries = " GFNn-xTB (Input geometry)"
            if self.solvent != "gas":
                info.append(
                    [
                        "justprint",
                        f"\nshort-notation:\n{tmp_func3_basis3} + "
                        f"{str(getattr(self, 'smgsolv3')).upper()}[{self.solvent}] "
                        f"+ GmRRHO({tmp_rrho_method}) // {geometries}",
                    ]
                )
            else:
                info.append(
                    [
                        "justprint",
                        f"\nshort-notation:\n{tmp_func3_basis3}"
                        f" + GmRRHO({str(getattr(self, 'part3_gfnv')).upper()}) "
                        f"// {geometries}",
                    ]
                )
        # NMR MODE
        if self.nmrmode:
            info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
            info.append(
                ["justprint", " NMR MODE SETTINGS".center(int(PLENGTH / 2), " ")]
            )
            info.append(["justprint", "".ljust(int(PLENGTH / 2), "-")])
            info.append(["part4", "part4"])
            info.append(["couplings", "calculate couplings (J)"])
            if self.couplings:
                info.append(["prog4_j", "program for coupling calculations"])
                if self.solvent != "gas":
                    info.append(["sm4_j", "solvation model for coupling calculations"])
                info.append(["func_j", "functional for coupling calculation"])
                info.append(["basis_j", "basis set for coupling calculation"])
            info.append(["justprint", ""])
            info.append(["shieldings", "calculate shieldings (S)"])
            if self.shieldings:
                info.append(["prog4_s", "program for shielding calculations"])
                if self.solvent != "gas":
                    info.append(["sm4_s", "solvation model for shielding calculations"])
                info.append(["func_s", "functional for shielding calculation"])
                info.append(["basis_s", "basis set for shielding calculation"])
            info.append(["justprint", ""])
            if getattr(self, "h_active"):
                info.append(["h_active", "Calculating proton spectrum"])
                info.append(["h_ref", "reference for 1H"])
            if getattr(self, "c_active"):
                info.append(["c_active", "Calculating carbon spectrum"])
                info.append(["c_ref", "reference for 13C"])
            if getattr(self, "f_active"):
                info.append(["f_active", "Calculating fluorine spectrum"])
                info.append(["f_ref", "reference for 19F"])
            if getattr(self, "si_active"):
                info.append(["si_active", "Calculating silicon spectrum"])
                info.append(["si_ref", "reference for 29Si"])
            if getattr(self, "p_active"):
                info.append(["p_active", "Calculating phosphorus spectrum"])
                info.append(["p_ref", "reference for 31P"])
            if all(
                [
                    not getattr(self, active)
                    for active in (
                        "h_active",
                        "c_active",
                        "f_active",
                        "si_active",
                        "p_active",
                    )
                ]
            ):
                info.append(
                    [
                        "printoption",
                        "Calculating spectrum for all nuclei",
                        "H, C, F, Si, P",
                    ]
                )
            info.append(["resonance_frequency", "resonance frequency"])
            # short notation:

        if self.optical_rotation:
            info.append(["justprint", "\n" + "".ljust(int(PLENGTH / 2), "-")])
            info.append(
                [
                    "justprint",
                    "OPTICAL ROTATION MODE - PART5".center(int(PLENGTH / 2), " "),
                ]
            )
            info.append(["justprint", "".ljust(int(PLENGTH / 2), "-")])
            info.append(["optical_rotation", "part5"])
            info.append(["freq_or", "frequency in [nm]"])
            info.append(["func_or_scf", "functional for SCF"])
            info.append(["func_or", "functional for optical rotation"])
            info.append(["basis_or", "basis set for optical rotation"])
            if not self.part3:
                info.append(["part2_P_threshold", "Boltzmann sum threshold employed"])
            elif self.part3:
                info.append(["part3_threshold", "Boltzmann sum threshold employed"])

        # TODO - if-branch for uv/vis

        max_len_digilen = 0
        for item in info:
            if item[0] == "justprint":
                if "short-notation" in item[1]:
                    continue
                    # tmp = len(item[1]) -len('short-notation:')
                else:
                    tmp = len(item[1])
            else:
                tmp = len(item[1])
            if tmp > max_len_digilen:
                max_len_digilen = tmp
        max_len_digilen += 1
        if max_len_digilen < DIGILEN:
            max_len_digilen = DIGILEN

        optionsexchange = {True: "on", False: "off"}
        for item in info:
            if item[0] == "justprint":
                # print everything after justprint
                print(item[1:][0])
            else:
                if item[0] == "printoption":
                    option = item[2]
                else:
                    option = getattr(self, item[0])
                if option is True or option is False:
                    option = optionsexchange[option]
                elif isinstance(option, list):
                    option = [str(i) for i in option]
                    if len(str(option)) > 40:
                        length = 0
                        reduced = []
                        for i in option:
                            length += len(i) + 2
                            if length < 40:
                                reduced.append(i)
                        reduced.append("...")
                        option = reduced
                        length = 0
                    option = ", ".join(option)
                print(
                    "{}: {:{digits}} {}".format(
                        item[1], "", option, digits=max_len_digilen - len(item[1])
                    )
                )
        print("END of parameters\n")

    def create_SI(self, ensembledata):
        """Automatically create a supporting information block"""
        si_out = []
        l_length = 25
        si_out.append("\nAutomatic CENSO-SI creation:\n")
        si_out.append(
            "The preparation of this supporting information (SI) is by no means complete and "
            + "does not relieve the\nuser of the responsibility to verify/complete the SI."
        )
        si_out.append("In this CENSO run the following settings were employed:\n")
        # General information:
        tmp = {
            "CENSO Version": __version__,
            "Temperature": str(self.temperature)
            + "K used in sorting and Boltzmann factor evaluation",
            "Solvent": self.solvent,
        }
        si_out.append(f"{'':^{l_length}} | {'General information:'}")
        si_out.append("".ljust(int(PLENGTH), "-"))
        for key, value in tmp.items():
            si_out.append(f"{key:^{l_length}} | {value}")
        si_out.append("".ljust(int(PLENGTH), "-"))
        # ----------------------------PART0--------------------------------------
        if self.part0:
            # part0: E = fast single-point with normally b97-d/def2-SV(P)
            #        GRRHO is omitted for speed!
            #        Gsolv = ALPB_GSOLV
            si_out.append(f"{'':^{l_length}} | {'PART0 - cheap prescreening'}")
            si_out.append("".ljust(int(PLENGTH), "-"))
            for key, value in getattr(ensembledata, "si")["part0"].items():
                si_out.append(f"{key:^{l_length}} | {value}")
            si_out.append("".ljust(int(PLENGTH), "-"))
        # ==========================END PART0====================================
        # ----------------------------PART1--------------------------------------
        if self.part1:
            # part1: E normally with r2scan-3c
            # Gsolv can be anything
            # GmRRHO xtb with GFNn-xTB + ALPB GBSA gas-phase
            ###
            si_out.append(f"{'':^{l_length}} | {'PART1 - prescreening'}")
            si_out.append("".ljust(int(PLENGTH), "-"))
            for key, value in getattr(ensembledata, "si")["part1"].items():
                si_out.append(f"{key:^{l_length}} | {value}")
            si_out.append("".ljust(int(PLENGTH), "-"))
        # ==========================END PART1====================================
        # ----------------------------PART2--------------------------------------
        if self.part2:
            # part2: E normally with r2scan-3c
            # Gsolv can be anything
            # GmRRHO xtb with GFNn-xTB + ALPB GBSA gas-phase
            ###
            si_out.append(f"{'':^{l_length}} | {'PART2 - optimization'}")
            si_out.append("".ljust(int(PLENGTH), "-"))
            for key, value in getattr(ensembledata, "si")["part2"].items():
                si_out.append(f"{key:^{l_length}} | {value}")
            si_out.append("".ljust(int(PLENGTH), "-"))
        # ==========================END PART2====================================
        # ----------------------------PART3--------------------------------------
        if self.part3:
            # # part3: E normally hybrid dft
            # # Gsolv can be anything
            # # GmRRHO xtb with GFNn-xTB + ALPB GBSA gas-phase
            si_out.append(f"{'':^{l_length}} | {'PART3 - refinement'}")
            si_out.append("".ljust(int(PLENGTH), "-"))
            for key, value in getattr(ensembledata, "si")["part3"].items():
                si_out.append(f"{key:^{l_length}} | {value}")
            si_out.append("".ljust(int(PLENGTH), "-"))
        # ==========================END PART3====================================
        # ----------------------------PART4--------------------------------------
        if self.part4:
            ###
            si_out.append(f"{'':^{l_length}} | {'PART4 - NMR mode'}")
            si_out.append("".ljust(int(PLENGTH), "-"))
            for key, value in getattr(ensembledata, "si")["part4"].items():
                si_out.append(f"{key:^{l_length}} | {value}")
            si_out.append("".ljust(int(PLENGTH), "-"))
        # ==========================END PART4====================================
        # ----------------------------PART5--------------------------------------
        if self.optical_rotation:
            ###
            si_out.append(f"{'':^{l_length}} | {'PART5 - OR mode'}")
            si_out.append("".ljust(int(PLENGTH), "-"))
            for key, value in getattr(ensembledata, "si")["part5"].items():
                si_out.append(f"{key:^{l_length}} | {value}")
            si_out.append("".ljust(int(PLENGTH), "-"))
        # ==========================END PART5====================================
        # PRINT SI PART information
        for line in si_out:
            print(line)
        si_text = """
        Within CENSO free energies (G) are calculated from:

        G(T) = E + G_mRRHO(T) + G_solv(T)                                 (1)

        The input ensemble (originating from CREST or crest_combi) is sorted 
        in part0 and part1 at DFT level on GFNn-xTB input geometries. This provides an 
        improved description of the electronic energy and solvation contributions 
        compared to the input SQM level. Thermostatistical contributions, including 
        ZPVE are evaluated at the GFNn-xTB level in part1. This first filtering step 
        removes structures high lying in free energy and reduces the computational 
        cost by passing fewer conformers to the computational costly DFT optimizations.

        COSMO-RS calculations are performed with the energies and densities from 
        the functional and basis set combination of the 'energy' evaluation in the 
        respective part.

        """
        print(si_text)
        print("\nSome citations in bib style are provided below:")
        # PRINT bib of employed programs
        from .cfg import si_bib

        out_bib = []
        requirements = self.needed_external_programs()
        bib_need = {
            "needtm": "tm",
            "needxtb": "xtb",
            "needcosmors": "cosmors",
            "needorca": "orca",
        }
        out_bib.extend(si_bib.get("censo"))
        for item in ("needtm", "needxtb", "needcosmors", "needorca"):
            if requirements.get(item, False):
                out_bib.extend(si_bib.get(bib_need[item], []))
        for item in set([self.func0, self.func, self.func3, self.func_j, self.func_s]):
            if si_bib.get(item, None) is not None:
                out_bib.extend(si_bib.get(item, []))
        if self.evaluate_rrho and self.bhess:
            out_bib.extend(si_bib.get("sph", []))
        if (self.solvent != 'gas' and 
            any(
                [x in ('alpb', 'alpb_gsolv', 'gbsa', 'gbsa_gsolv') 
                for x in (self.sm_rrho, self.smgsolv1, self.smgsolv2, self.smgsolv3)]
                )
            ):
            out_bib.extend(si_bib.get("alpb", []))
        print("\nBib entries:")
        for line in out_bib:
            print(line)


    def read_rcfile():
        """read rcfile into buffer (supposed to be decorator)"""


    

    def needed_external_programs(self):
        """
        Automatically checks which external programs are required for the
        current run.
        """
        requirements = {}
        # xTB
        if (
            self.prog_rrho == "xtb"
            or self.part0
            or self.ancopt
            or self.smgsolv2 in ("gbsa_gsolv", "alpb_gsolv")
        ):
            requirements["needxtb"] = True
        # TM
        if (
            self.prog == "tm"
            or (self.prog2opt == "tm" and self.part2)
            or (self.prog3 == "tm" and self.part3)
            or (self.prog4_j == "tm" and self.part4 and self.couplings)
            or (self.prog4_s == "tm" and self.part4 and self.couplings)
            or (self.smgsolv1 in ("cosmors", "cosmors-fine") and self.part1)
            or (self.smgsolv2 in ("cosmors", "cosmors-fine") and self.part2)
            or (self.smgsolv3 in ("cosmors", "cosmors-fine") and self.part3)
        ):
            requirements["needtm"] = True
            requirements["needcefine"] = True
        if self.part4:
            if self.couplings and self.prog4_j == "tm":
                requirements["needescf"] = True
            if self.shieldings and self.prog4_s == "tm":
                requirements["needmpshift"] = True
        # COSMORS
        if (
            (self.part1 and self.smgsolv1 in  "cosmors")
            or (self.part2 and self.smgsolv2 == "cosmors")
            or (self.part3 and self.smgsolv3 == "cosmors")
        ):
            requirements["needcosmors"] = True
            requirements["needcosmors-normal"] = True
        if (
            (self.part1 and self.smgsolv1 == "cosmors-fine")
            or (self.part2 and self.smgsolv2 == "cosmors-fine")
            or (self.part3 and self.smgsolv3 == "cosmors-fine")
        ):
            requirements["needcosmors"] = True
            requirements["needcosmors-fine"] = True
        # ORCA
        if (
            self.prog == "orca"
            or (self.prog2opt == "orca" and self.part2)
            or (self.prog3 == "orca" and self.part3)
            or (self.prog4_j == "orca" and self.part4 and self.couplings)
            or (self.prog4_s == "orca" and self.part4 and self.shieldings)
            or (self.smgsolv1 == "smd_gsolv" and self.part1)
            or (self.smgsolv2 == "smd_gsolv" and self.part2)
            or (self.smgsolv3 == "smd_gsolv" and self.part3)
        ):
            requirements["needorca"] = True
        if self.run:
            requirements["startenso"] = True
        return requirements

    def _updateEnvironsettings(self, newsettings=None):
        """
        Update the environmentsettings which is needed for e.g. Turbomole
        calculations and is provided in each subroutine call.
        """
        if newsettings is not None:
            for key, value in newsettings.items():
                ENVIRON[key] = str(value)

    def processQMpaths(self, requirements, error_logical):
        """
        print path at startup and return error if programs don't exist
        """
        # print relevant Program paths:
        print("\n" + "".ljust(DIGILEN, "-"))
        print("PATHS of external QM programs".center(DIGILEN, " "))
        print("".ljust(DIGILEN, "-") + "\n")
        print("The following program paths are used:")
        if requirements.get("needorca", False):
            print("    ORCA:         {}".format(self.external_paths["orcapath"]))
            print("    ORCA Version: {}".format(self.external_paths["orcaversion"]))
        if requirements.get("needxtb", False):
            print("    xTB:          {}".format(self.external_paths["xtbpath"]))
        if requirements.get("needcrest", False):
            print("    CREST:        {}".format(self.external_paths["crestpath"]))
        if requirements.get("needtm", False):
            tmpath = shutil.which("ridft")
            if tmpath is not None:
                tmpath = os.path.dirname(tmpath)
            else:
                tmpath = "None"
            print("    TURBOMOLE:    {}".format(tmpath))
        if requirements.get("needescf", False):
            print("    escf:         {}".format(self.external_paths["escfpath"]))
        if requirements.get("needmpshift", False):
            print("    mpshift:      {}".format(self.external_paths["mpshiftpath"]))
        if requirements.get("needcosmors", False):
            try:
                tmp = self.external_paths["cosmorssetup"].split()
                if len(tmp) == 9:
                    print("    Setup of COSMO-RS:")
                    if self.cosmorsparam == "automatic":
                        print("        {}".format(" ".join(tmp[0:3])))  # ctd
                    else:
                        print(
                            "        ctd = {}".format(
                                cosmors_param.get(self.cosmorsparam, " ".join(tmp[0:3]))
                            )
                        )
                    print("        {}".format(" ".join(tmp[3:6])))  # cdir
                    print("        {}".format(" ".join(tmp[6:9])))  # ldir
                else:
                    print(
                        f"    Setup of COSMO-RS: {str(self.external_paths['cosmorssetup'])}"
                    )
            except Exception:
                print(
                    "    Setup of COSMO-RS: {}".format(
                        str(self.external_paths["cosmorssetup"])
                    )
                )
            if requirements.get("needcosmors-fine", False):
                # FINE
                db_path = os.path.join(
                    self.external_paths["dbpath"], "DATABASE-COSMO/BP-TZVPD-FINE"
                )
                print(
                    f"    Using {db_path}\n"
                    "    as path to the COSMO-RS FINE DATABASE."
                )
            if requirements.get("needcosmors-normal", False):
                # NORMAL
                db_path = os.path.join(
                    self.external_paths["dbpath"], "DATABASE-COSMO/BP-TZVP-COSMO"
                )
                print(
                    f"    Using {db_path}\n"
                    "    as path to the COSMO-RS NORMAL DATABASE."
                )
        print("")
        # Check if paths of needed programs exist:
        if requirements.get("needcrest", False):
            if (
                self.external_paths["crestpath"] is None
                or shutil.which(self.external_paths["crestpath"]) is None
            ):
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}path for CREST is not correct!")
                error_logical = True
        # xTB
        if requirements.get("needxtb", False):
            if (
                self.external_paths["xtbpath"] is None
                or shutil.which(self.external_paths["xtbpath"]) is None
            ):
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}path for xTB is not correct!")
                error_logical = True
            try:
                ENVIRON["OMP_NUM_THREADS"] = "{:d}".format(self.omp)
            except Exception:
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}can not set omp for xTB calculation!")
        # ORCA
        if requirements.get("needorca", False):
            if (
                self.external_paths["orcapath"] is None
                or shutil.which(os.path.join(self.external_paths["orcapath"], "orca"))
                is None
            ):
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}path for ORCA is not correct!")
                error_logical = True
        # cefine
        if requirements.get("needcefine", False):
            if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
                # print('running in a PyInstaller bundle')
                bundle_dir = getattr(
                    sys, "_MEIPASS", os.path.abspath(os.path.dirname(__file__))
                )
                path_to_cefine = os.path.abspath(os.path.join(bundle_dir, "cefine"))
                if not os.path.exists(path_to_cefine):
                    path_to_cefine = shutil.which("cefine")
            else:
                # print('running in a normal Python process')
                path_to_cefine = shutil.which("cefine")

            if path_to_cefine is not None and os.path.exists(path_to_cefine):
                print("    Using cefine from {}".format(path_to_cefine))
                self.external_paths["cefinepath"] = path_to_cefine
            else:
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}cefine (the commandline program for define) has not been found!"
                )
                self.save_errors.append(f"{'':{WARNLEN}}all programs needing TM can not start!")
                error_logical = True
        # TM
        if requirements.get("needtm", False):
            # preparation of parallel calculation with TM
            try:
                if ENVIRON.get("PARA_ARCH", None) == "SMP":
                    try:
                        ENVIRON["PARNODES"] = str(self.omp)
                        ENVIRON["OMP_NUM_THREADS"] = "{:d}".format(self.omp)
                        print(
                            "    PARNODES for TM or COSMO-RS calculation was set "
                            "to {}".format(ENVIRON["PARNODES"])
                        )
                    except Exception:
                        print(f"{'ERROR:':{WARNLEN}}PARNODES can not be changed!")
                        error_logical = True
                        raise
                else:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}PARA_ARCH has to be set to SMP for parallel TM "
                        "calculations!"
                    )
                    if self.run:
                        error_logical = True
            except Exception:
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}PARA_ARCH has to be set to SMP and PARNODES have to "
                    f"be set\n{'':{WARNLEN}}for parallel TM calculations!."
                )
                if requirements.get("startenso", False):
                    error_logical = True
                raise
        if requirements.get("needescf", False):
            if (
                self.external_paths["escfpath"] is None
                or shutil.which(self.external_paths["escfpath"]) is None
            ):
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}path for escf is not correct!")
                error_logical = True
        if requirements.get("needmpshift", False):
            if (
                self.external_paths["mpshiftpath"] is None
                or shutil.which(self.external_paths["mpshiftpath"]) is None
            ):
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}path for mpshift is not correct!")
                error_logical = True
        # COSMORS
        if requirements.get("needcosmors", False):
            if self.external_paths["cosmorssetup"] is None:
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Set up for COSMO-RS has to be written to .censorc!"
                )
                error_logical = True
            if self.external_paths["cosmothermversion"] is None:
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Version of COSMO-RS has to be written to .censorc!"
                )
                error_logical = True
            if shutil.which("cosmotherm") is not None:
                print("    Using COSMOtherm from {}".format(shutil.which("cosmotherm")))
            else:
                self.save_errors.append(f"{'ERROR:':{WARNLEN}}COSMOtherm has not been found!")
                error_logical = True
        # update cfg.external paths
        external_paths.update(self.external_paths)
        return error_logical

    def write_rcfile(self, pathtofile, usepaths=False, update=False):
        """
        write new global configruation file into the current directroy.
        """
        args_key = {v: k for k, v in self.key_args_dict.items()}
        if update:
            data = self.provide_runinfo(extend=False)
        else:
            data = {}
        with open(pathtofile, "w", newline=None) as outdata:
            outdata.write("$CENSO global configuration file: .censorc\n")
            outdata.write(f"$VERSION:{__version__} \n")
            outdata.write("\n")
            if usepaths:
                # write stored program paths to file
                outdata.write(f"ORCA: {self.external_paths['orcapath']}\n")
                outdata.write(f"ORCA version: {self.external_paths['orcaversion']}\n")
                outdata.write(f"GFN-xTB: {self.external_paths['xtbpath']}\n")
                outdata.write(f"CREST: {self.external_paths['crestpath']}\n")
                outdata.write(f"mpshift: {self.external_paths['mpshiftpath']}\n")
                outdata.write(f"escf: {self.external_paths['escfpath']}\n")
                outdata.write("\n")
                outdata.write("#COSMO-RS\n")
                outdata.write(f"{self.external_paths['cosmorssetup']}\n")
                # outdata.write("cosmothermversion: 16\n")
            else:
                outdata.write("ORCA: /path/excluding/binary/\n")
                outdata.write("ORCA version: 4.2.1\n")
                outdata.write("GFN-xTB: /path/including/binary/xtb-binary\n")
                outdata.write("CREST: /path/including/binary/crest-binary\n")
                outdata.write("mpshift: /path/including/binary/mpshift-binary\n")
                outdata.write("escf: /path/including/binary/escf-binary\n")
                outdata.write("\n")
                outdata.write("#COSMO-RS\n")
                outdata.write(
                    "ctd = BP_TZVP_C30_1601.ctd cdir = "
                    '"/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES" ldir = '
                    '"/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES"\n'
                )
                # outdata.write("cosmothermversion: 16\n")
            outdata.write("$ENDPROGRAMS\n\n")
            outdata.write("$CRE SORTING SETTINGS:\n")
            outdata.write("$GENERAL SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_general):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_general)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                if key == "nconf" and value is None:
                    value = "all"
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART0 - CHEAP-PRESCREENING - SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part0):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part0)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART1 - PRESCREENING - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part1):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part1)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART2 - OPTIMIZATION - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part2):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part2)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$PART3 - REFINEMENT - SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part3):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part3)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$NMR PROPERTY SETTINGS:\n")
            outdata.write("$PART4 SETTINGS:\n")
            for key in OrderedDict(self.defaults_nmrprop_part4):
                value = self._exchange_onoff(
                    data.get(
                        key, OrderedDict(self.defaults_nmrprop_part4)[key]["default"]
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("\n$OPTICAL ROTATION PROPERTY SETTINGS:\n")
            outdata.write("$PART5 SETTINGS:\n")
            for key in OrderedDict(self.defaults_optical_rotation_part5):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_optical_rotation_part5)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                options = self.value_options.get(key, "possibilities")
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, options))
            outdata.write("$END CENSORC\n")

    def write_censo_inp(self, path=None):
        """
         Write file "censo.inp" which stores current settings in the .censorc format
         """
        if path is None:
            path = self.cwd
        args_key = {v: k for k, v in self.key_args_dict.items()}
        data = self.provide_runinfo(extend=False)
        with open(os.path.join(path, "censo.inp"), "w", newline=None) as outdata:
            outdata.write("$File: censo.inp settings of current calculation\n")
            outdata.write(f"$VERSION:{__version__} \n")
            outdata.write("\n")
            # write stored program paths to file
            outdata.write(f"ORCA: {self.external_paths['orcapath']}\n")
            outdata.write(f"ORCA version: {self.external_paths['orcaversion']}\n")
            outdata.write(f"GFN-xTB: {self.external_paths['xtbpath']}\n")
            outdata.write(f"CREST: {self.external_paths['crestpath']}\n")
            outdata.write(f"mpshift: {self.external_paths['mpshiftpath']}\n")
            outdata.write(f"escf: {self.external_paths['escfpath']}\n")
            outdata.write("\n")
            outdata.write("#COSMO-RS\n")
            outdata.write(f"{self.external_paths['cosmorssetup']}\n")
            # outdata.write("cosmothermversion: 16\n")
            outdata.write("$ENDPROGRAMS\n\n")
            outdata.write("$CRE SORTING SETTINGS:\n")
            outdata.write("$GENERAL SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_general):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_general)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                if key == "nconf" and value is None:
                    value = "all"
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART0 - CHEAP-PRESCREENING - SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part0):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part0)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART1 - PRESCREENING - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part1):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part1)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART2 - OPTIMIZATION - SETTINGS:\n")
            outdata.write("# func and basis is set under GENERAL SETTINGS\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part2):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part2)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$PART3 - REFINEMENT - SETTINGS:\n")
            for key in OrderedDict(self.defaults_refine_ensemble_part3):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_refine_ensemble_part3)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$NMR PROPERTY SETTINGS:\n")
            outdata.write("$PART4 SETTINGS:\n")
            for key in OrderedDict(self.defaults_nmrprop_part4):
                value = self._exchange_onoff(
                    data.get(
                        key, OrderedDict(self.defaults_nmrprop_part4)[key]["default"]
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$OPTICAL ROTATION PROPERTY SETTINGS:\n")
            outdata.write("$PART5 SETTINGS:\n")
            for key in OrderedDict(self.defaults_optical_rotation_part5):
                value = self._exchange_onoff(
                    data.get(
                        key,
                        OrderedDict(self.defaults_optical_rotation_part5)[key][
                            "default"
                        ],
                    ),
                    reverse=True,
                )
                key = args_key.get(key, key)
                outdata.write(format_line(key, value, ""))
            outdata.write("\n$END censo.inp\n")

    ##########################

    