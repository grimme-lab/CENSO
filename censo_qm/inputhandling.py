"""
defininition of internal defaults, checking of logic for parameter combinations,
cml parsing
"""
import argparse
import shutil
import os
import json
import csv
import time
import math
import sys
import glob
from copy import deepcopy
from collections import OrderedDict
from typing import Text
from .cfg import (
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
from .utilities import frange, format_line, print, print_block


def cml(startup_description, options, argv=None):
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
        action="store",
        required=False,
        metavar="",
        help="Charge of the investigated molecule.",
    )
    group1.add_argument(
        "-u",
        "--unpaired",
        dest="unpaired",
        action="store",
        required=False,
        type=int,
        metavar="",
        help="Integer number of unpaired electrons of the investigated molecule.",
    )
    group1.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
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
        help="Solvent the molecule is solvated in, available solvents "
        "are: {}. They can be extended in the "
        "file ~/.censo_assets/censo_solvents.json .".format(
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
        "of conformers is left. (never exceed O*P cores).",
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
    args = parser.parse_args(argv)

    if args.part3only:
        setattr(args, "part0", "off")
        setattr(args, "part1", "off")
        setattr(args, "part2", "off")
    return args


class internal_settings:
    """
    All options are saved here.
    """

    # key in .censorc corresponds to name in cml
    key_args_dict = {
        "nconf": "nconf",
        "charge": "charge",
        "unpaired": "unpaired",
        "solvent": "solvent",
        "prog": "prog",
        "ancopt": "ancopt",
        "opt_spearman": "opt_spearman",
        "evaluate_rrho": "evaluate_rrho",
        "consider_sym": "consider_sym",
        "prog_rrho": "prog_rrho",
        "part0_gfnv": "part0_gfnv",
        "part1_gfnv": "part1_gfnv",
        "part2_gfnv": "part2_gfnv",
        "part3_gfnv": "part3_gfnv",
        "temperature": "temperature",
        "multitemp": "multitemp",
        "trange": "trange",
        "prog2opt": "prog2opt",
        "prog3": "prog3",
        "prog4_j": "prog4_j",
        "prog4_s": "prog4_s",
        "part0": "part0",
        "part1": "part1",
        "part2": "part2",
        "part3": "part3",
        "part4": "part4",
        "func0": "func0",
        "func": "func",
        "func3": "func3",
        "basis0": "basis0",
        "basis": "basis",
        "basis3": "basis3",
        "couplings": "couplings",
        "progJ": "prog4_j",
        "funcJ": "func_j",
        "basisJ": "basis_j",
        "shieldings": "shieldings",
        "progS": "prog4_s",
        "funcS": "func_s",
        "basisS": "basis_s",
        "part0_threshold": "part0_threshold",
        "part1_threshold": "part1_threshold",
        "part3_threshold": "part3_threshold",
        #"part2_threshold": "part2_threshold",
        #"opt_limit": "opt_limit",
        "part2_threshold": "opt_limit",            # updated to avoid confusion???
        "part2_P_threshold": "part2_P_threshold",     # updated to avoid confusion???
        "smgsolv1": "smgsolv1",
        "sm2": "sm2",
        "smgsolv2": "smgsolv2",
        "smgsolv3": "smgsolv3",
        "sm4J": "sm4_j",
        "sm4S": "sm4_s",
        "check": "check",
        "crestcheck": "crestcheck",
        "maxthreads": "maxthreads",
        "omp": "omp",
        "1H_active": "h_active",
        "13C_active": "c_active",
        "19F_active": "f_active",
        "31P_active": "p_active",
        "29Si_active": "si_active",
        "resonance_frequency": "resonance_frequency",
        "reference_1H": "h_ref",
        "reference_13C": "c_ref",
        "reference_31P": "p_ref",
        "reference_19F": "f_ref",
        "reference_29Si": "si_ref",
        "bhess": "bhess",
        "sm_rrho": "sm_rrho",
        "optcycles": "optcycles",
        "optlevel2": "optlevel2",
        "spearmanthr": "spearmanthr",
        "optical_rotation": "optical_rotation",
        "radsize": "radsize",
        "frequency_optical_rot": "freq_or",
        "funcOR": "func_or",
        "basisOR": "basis_or",
        "funcOR_SCF": "func_or_scf",
        "hlow": "hlow",
        "rmsdbias": "rmsdbias",
        "imagthr": "imagthr",
        "sthr": "sthr",
        "scale": "scale",
        "cosmorsparam": "cosmorsparam",
        "progress": "progress",
        "balance": "balance",
    }

    impgfnv = ["gfn1", "gfn2", "gfnff"]
    tmp_smd_solvents = [
        "1,1,1-TRICHLOROETHANE",
        "1,1,2-TRICHLOROETHANE",
        "1,2,4-TRIMETHYLBENZENE",
        "1,2-DIBROMOETHANE",
        "1,2-DICHLOROETHANE",
        "1,2-ETHANEDIOL",
        "1,4-DIOXANE",
        "1-BROMO-2-METHYLPROPANE",
        "1-BROMOOCTANE",
        "1-BROMOPENTANE",
        "1-BROMOPROPANE",
        "1-BUTANOL",
        "1-CHLOROHEXANE",
        "1-CHLOROPENTANE",
        "1-CHLOROPROPANE",
        "1-DECANOL",
        "1-FLUOROOCTANE",
        "1-HEPTANOL",
        "1-HEXANOL",
        "1-HEXENE",
        "1-HEXYNE",
        "1-IODOBUTANE",
        "1-IODOHEXADECANE",
        "1-IODOPENTANE",
        "1-IODOPROPANE",
        "1-NITROPROPANE",
        "1-NONANOL",
        "1-OCTANOL",
        "1-PENTANOL",
        "1-PENTENE",
        "1-PROPANOL",
        "2,2,2-TRIFLUOROETHANOL",
        "2,2,4-TRIMETHYLPENTANE",
        "2,4-DIMETHYLPENTANE",
        "2,4-DIMETHYLPYRIDINE",
        "2,6-DIMETHYLPYRIDINE",
        "2-BROMOPROPANE",
        "2-BUTANOL",
        "2-CHLOROBUTANE",
        "2-HEPTANONE",
        "2-HEXANONE",
        "2-METHOXYETHANOL",
        "2-METHYL-1-PROPANOL",
        "2-METHYL-2-PROPANOL",
        "2-METHYLPENTANE",
        "2-METHYLPYRIDINE",
        "2-NITROPROPANE",
        "2-OCTANONE",
        "2-PENTANONE",
        "2-PROPANOL",
        "2-PROPEN-1-OL",
        "E-2-PENTENE",
        "3-METHYLPYRIDINE",
        "3-PENTANONE",
        "4-HEPTANONE",
        "4-METHYL-2-PENTANONE",
        "4-METHYLPYRIDINE",
        "5-NONANONE",
        "ACETIC ACID",
        "ACETONE",
        "ACETONITRILE",
        "ACETOPHENONE",
        "ANILINE",
        "ANISOLE",
        "BENZALDEHYDE",
        "BENZENE",
        "BENZONITRILE",
        "BENZYL ALCOHOL",
        "BROMOBENZENE",
        "BROMOETHANE",
        "BROMOFORM",
        "BUTANAL",
        "BUTANOIC ACID",
        "BUTANONE",
        "BUTANONITRILE",
        "BUTYL ETHANOATE",
        "BUTYLAMINE",
        "N-BUTYLBENZENE",
        "SEC-BUTYLBENZENE",
        "TERT-BUTYLBENZENE",
        "CARBON DISULFIDE",
        "CARBON TETRACHLORIDE",
        "CHLOROBENZENE",
        "CHLOROFORM",
        "A-CHLOROTOLUENE",
        "O-CHLOROTOLUENE",
        "M-CRESOL",
        "O-CRESOL",
        "CYCLOHEXANE",
        "CYCLOHEXANONE",
        "MeCN",
        "CCl4",
        "CYCLOPENTANE",
        "CYCLOPENTANOL",
        "CYCLOPENTANONE",
        "DECALIN (CIS/TRANS MIXTURE)",
        "CIS-DECALIN",
        "N-DECANE",
        "DIBROMOMETHANE",
        "DIBUTYLETHER",
        "O-DICHLOROBENZENE",
        "E-1,2-DICHLOROETHENE",
        "Z-1,2-DICHLOROETHENE",
        "DICHLOROMETHANE",
        "DIETHYL ETHER",
        "DIETHYL SULFIDE",
        "DIETHYLAMINE",
        "DIIODOMETHANE",
        "DIISOPROPYL ETHER",
        "CIS-1,2-DIMETHYLCYCLOHEXANE",
        "DIMETHYL DISULFIDE",
        "N,N-DIMETHYLACETAMIDE",
        "N,N-DIMETHYLFORMAMIDE",
        "DIMETHYLSULFOXIDE",
        "DIPHENYLETHER",
        "DIPROPYLAMINE",
        "N-DODECANE",
        "ETHANETHIOL",
        "ETHANOL",
        "ETHYL ETHANOATE",
        "ETHYL METHANOATE",
        "ETHYL PHENYL ETHER",
        "ETHYLBENZENE",
        "FLUOROBENZENE",
        "FORMAMIDE",
        "FORMIC ACID",
        "N-HEPTANE",
        "N-HEXADECANE",
        "N-HEXANE",
        "HEXANOIC ACID",
        "IODOBENZENE",
        "IODOETHANE",
        "IODOMETHANE",
        "ISOPROPYLBENZENE",
        "P-ISOPROPYLTOLUENE",
        "MESITYLENE",
        "METHANOL",
        "METHYL BENZOATE",
        "METHYL BUTANOATE",
        "METHYL ETHANOATE",
        "METHYL METHANOATE",
        "METHYL PROPANOATE",
        "N-METHYLANILINE",
        "METHYLCYCLOHEXANE",
        "N-METHYLFORMAMIDE",
        "NITROBENZENE",
        "NITROETHANE",
        "NITROMETHANE",
        "O-NITROTOLUENE",
        "N-NONANE",
        "N-OCTANE",
        "N-PENTADECANE",
        "PENTANAL",
        "N-PENTANE",
        "PENTANOIC ACID",
        "PENTYL ETHANOATE",
        "PENTYLAMINE",
        "PERFLUOROBENZENE",
        "PROPANAL",
        "PROPANOIC ACID",
        "PROPANONITRILE",
        "PROPYL ETHANOATE",
        "PROPYLAMINE",
        "PYRIDINE",
        "TETRACHLOROETHENE",
        "TETRAHYDROFURAN",
        "TETRAHYDROTHIOPHENE-S,S-DIOXIDE",
        "TETRALIN",
        "THIOPHENE",
        "THIOPHENOL",
        "TOLUENE",
        "TRANS-DECALIN",
        "TRIBUTYLPHOSPHATE",
        "TRICHLOROETHENE",
        "TRIETHYLAMINE",
        "N-UNDECANE",
        "WATER",
        "XYLENE (MIXTURE)",
        "M-XYLENE",
        "O-XYLENE",
        "P-XYLENE",
        "DMF",
        "DMSO",
        "PhNO2",
        "MeNO2",
        "THF",
    ]
    solvents_smd = [i.lower() for i in tmp_smd_solvents]
    solvents_xtb = [
        "acetone",
        "acetonitrile",
        "aniline",
        "benzaldehyde",
        "benzene",
        "chcl3",
        "ch2cl2",
        "ccl4",
        "cs2",
        "dioxane",
        "dmf",
        "dmso",
        "ether",
        "ethanol",
        "ethylacetate",
        "furane",
        "hexadecane",
        "hexane",
        "h2o",
        "water",
        "methanol",
        "nitromethane",
        "thf",
        "toluene",
        "octanol",
        "woctanol",
        "phenol",
    ]
    solvents_cpcm = [
        "water",
        "acetone",
        "acetonitrile",
        "ammonia",
        "benzene",
        "chloroform",
        "ch2cl2",
        "ccl4",
        "cyclohexane",
        "dmf",
        "dmso",
        "ethanol",
        "hexane",
        "methanol",
        "octanol",
        "pyridine",
        "thf",
        "toluene",
    ]
    solvents_cosmors = [
        "propanone_c0",
        "chcl3_c0",
        "acetonitrile_c0",
        "ch2cl2_c0",
        "dimethylsulfoxide_c0",
        "h2o_c0",
        "methanol_c0",
        "thf_c0",
        "toluene_c0",
        "1-octanol_c0",
        "woctanol",  # this is a mixture and treated differently
        "n-hexadecane_c0",
        "dimethylformamide_c0",
        "aniline_c0",
        "cyclohexane_c0",
        "ccl4_c0",
        "diethylether_c0",
        "ethanol_c0",
        "hexane_c0",
        "nitromethane_c0",
        "benzaldehyde_c0",
        "benzene_c0",
        "cs2_c0",
        "dioxane_c0",
        "ethylacetate_c0",
        "furane_c0",
        "phenol_c0",
        "1,2-dichloroethane_c0",
    ]

    # only using the dielectric constant (DC) for cosmo

    # dcosmorsfile name = e.g. acetonitrile + '_25.pot'
    solvents_dcosmors = [
        "acetonitrile",
        "aniline",
        "benzene",
        "ccl4",
        "chcl3",
        "cyclohexane",
        "diethylether",
        "dimethylsulfoxide",
        "ethanol",
        "h2o",
        "hexadecane",
        "hexane",
        "methanol",
        "nitromethane",
        "octanol",
        "propanone",
        "thf",
        "toluene",
        "wet-octanol",
    ]

    smgsolv_1 = ["cosmors", "cosmors-fine", "gbsa_gsolv", "alpb_gsolv", "smd_gsolv"]
    sm2_tm = ["cosmo", "dcosmors"]
    sm2_orca = ["cpcm", "smd"]
    smgsolv_2 = ["cosmors", "cosmors-fine", "gbsa_gsolv", "alpb_gsolv", "smd_gsolv"]
    smgsolv3_tm = ["cosmo", "dcosmors"]
    smgsolv3_orca = ["cpcm", "smd"]
    smgsolv_3 = ["cosmors", "cosmors-fine", "gbsa_gsolv", "alpb_gsolv", "smd_gsolv"]
    sm4_j_tm = ["cosmo", "dcosmors"]
    sm4_s_tm = ["cosmo", "dcosmors"]
    sm4_j_orca = ["cpcm", "smd"]
    sm4_s_orca = ["cpcm", "smd"]

    imphref = ["TMS"]
    impcref = ["TMS"]
    impfref = ["CFCl3"]
    imppref = ["TMP", "PH3"]
    impsiref = ["TMS"]

    def __init__(self, **kwargs):
        self.func_info = dfa_settings()
        self.impsm2 = list(set(self.sm2_orca + self.sm2_tm + ["default"]))
        self.impsmgsolv1 = list(
            set(self.sm2_orca + self.sm2_tm + self.smgsolv_2 + ["sm2"])
        )
        self.impsmgsolv2 = list(
            set(self.sm2_orca + self.sm2_tm + self.smgsolv_2 + ["sm2"])
        )
        self.impsmgsolv3 = list(
            set(self.sm2_orca + self.sm2_tm + self.smgsolv_2 + ["sm2"])
        )
        self.impsm4_j = list(set(self.sm4_j_orca + self.sm4_j_tm))
        self.impsm4_s = list(set(self.sm4_s_orca + self.sm4_s_tm))

        self.defaults_refine_ensemble_general = [
            # general settings
            ("nconf", {"default": None, "type": int}),
            ("charge", {"default": 0, "type": int}),
            ("unpaired", {"default": 0, "type": int}),
            ("solvent", {"default": "gas", "type": str}),
            ("prog_rrho", {"default": "xtb", "type": str}),
            ("temperature", {"default": 298.15, "type": float}),
            ("trange", {"default": [273.15, 378.15, 5], "type": list}),
            ("multitemp", {"default": True, "type": bool}),
            ("evaluate_rrho", {"default": True, "type": bool}),
            ("consider_sym", {"default": True, "type": bool}),
            ("bhess", {"default": True, "type": bool}),
            ("imagthr", {"default": "automatic", "type": str}),
            ("sthr", {"default": "automatic", "type": str}),
            ("scale", {"default": "automatic", "type": str}),
            ("rmsdbias", {"default": False, "type": bool}),
            ("sm_rrho", {"default": "alpb", "type": str}),
            ("progress", {"default": False, "type": bool}),
            ("check", {"default": True, "type": bool}),
            ("prog", {"default": "tm", "type": str}),
            ("func", {"default": "r2scan-3c", "type": str}),
            ("basis", {"default": "automatic", "type": str}),
            ("maxthreads", {"default": 1, "type": int}),
            ("omp", {"default": 1, "type": int}),
            ("balance", {"default": False, "type": bool}),
            ("cosmorsparam", {"default": "automatic", "type": str}),
        ]
        self.defaults_refine_ensemble_part0 = [
            # part0
            ("part0", {"default": True, "type": bool}),
            ("func0", {"default": "b97-d", "type": str}),
            ("basis0", {"default": "def2-SV(P)", "type": str}),
            ("part0_gfnv", {"default": "gfn2", "type": str}),
            ("part0_threshold", {"default": 4.0, "type": float}),
        ]
        self.defaults_refine_ensemble_part1 = [
            # part1
            ("part1", {"default": True, "type": bool}),
            ("smgsolv1", {"default": "cosmors", "type": str}),
            ("part1_gfnv", {"default": "gfn2", "type": str}),
            ("part1_threshold", {"default": 3.5, "type": float}),
        ]
        self.defaults_refine_ensemble_part2 = [
            # part2
            ("part2", {"default": True, "type": bool}),
            ("prog2opt", {"default": "prog", "type": str}),
            ("opt_limit", {"default": 2.5, "type": float}),
            ("sm2", {"default": "default", "type": str}),
            ("smgsolv2", {"default": "cosmors", "type": str}),
            ("part2_gfnv", {"default": "gfn2", "type": str}),
            ("ancopt", {"default": True, "type": bool}),
            ("hlow", {"default": 0.01, "type": float}),
            ("opt_spearman", {"default": True, "type": bool}),
            ("part2_P_threshold", {"default": 99, "type": float}),
            ("optlevel2", {"default": "automatic", "type": str}),
            ("optcycles", {"default": 8, "type": int}),
            ("spearmanthr", {"default": -4.0, "type": float}),
            ("radsize", {"default": 10, "type": int}),
            ("crestcheck", {"default": False, "type": bool}),
        ]
        self.defaults_refine_ensemble_part3 = [
            # part3
            ("part3", {"default": False, "type": bool}),
            ("prog3", {"default": "prog", "type": str}),
            ("func3", {"default": "pw6b95", "type": str}),
            ("basis3", {"default": "def2-TZVPD", "type": str}),
            ("smgsolv3", {"default": "cosmors", "type": str}),
            ("part3_gfnv", {"default": "gfn2 ", "type": str}),
            ("part3_threshold", {"default": 99, "type": float}),
        ]
        self.defaults_nmrprop_part4 = [
            # part4
            ("part4", {"default": False, "type": bool}),
            ("couplings", {"default": True, "type": bool}),
            ("prog4_j", {"default": "prog", "type": str}),
            ("func_j", {"default": "pbe0", "type": str}),
            ("basis_j", {"default": "def2-TZVP", "type": str}),
            ("sm4_j", {"default": "default", "type": str}),
            ("shieldings", {"default": True, "type": bool}),
            ("prog4_s", {"default": "prog", "type": str}),
            ("func_s", {"default": "pbe0", "type": str}),
            ("basis_s", {"default": "def2-TZVP", "type": str}),
            ("sm4_s", {"default": "default", "type": str}),
            ("h_ref", {"default": "TMS", "type": str}),
            ("c_ref", {"default": "TMS", "type": str}),
            ("f_ref", {"default": "CFCl3", "type": str}),
            ("si_ref", {"default": "TMS", "type": str}),
            ("p_ref", {"default": "TMP", "type": str}),
            ("h_active", {"default": True, "type": bool}),
            ("c_active", {"default": True, "type": bool}),
            ("f_active", {"default": False, "type": bool}),
            ("si_active", {"default": False, "type": bool}),
            ("p_active", {"default": False, "type": bool}),
            ("resonance_frequency", {"default": 300.0, "type": float}),
        ]
        self.defaults_optical_rotation_part5 = [
            # part5
            ("optical_rotation", {"default": False, "type": bool}),
            ("func_or", {"default": "pbe", "type": str}),
            ("func_or_scf", {"default": "r2scan-3c", "type": str}),
            ("basis_or", {"default": "def2-SVPD", "type": str}),
            ("freq_or", {"default": [589.0], "type": list}),
        ]

        self.internal_defaults = OrderedDict(
            self.defaults_refine_ensemble_general
            + self.defaults_refine_ensemble_part0
            + self.defaults_refine_ensemble_part1
            + self.defaults_refine_ensemble_part2
            + self.defaults_refine_ensemble_part3
            + self.defaults_nmrprop_part4
            + self.defaults_optical_rotation_part5
        )

        # update internal defaults specific to QM package
        # orca
        self.internal_defaults_orca = deepcopy(self.internal_defaults)
        self.internal_defaults_orca["func0"]["default"] = "b97-d3"
        self.internal_defaults_orca["sm2"]["default"] = "smd"
        self.internal_defaults_orca["smgsolv1"]["default"] = "smd"
        self.internal_defaults_orca["smgsolv2"]["default"] = "smd"
        self.internal_defaults_orca["smgsolv3"]["default"] = "smd"
        self.internal_defaults_orca["sm4_j"]["default"] = "smd"
        self.internal_defaults_orca["sm4_s"]["default"] = "smd"
        # when basis == automatic and not "-3c"
        self.internal_defaults_orca["basis0"]["default"] = "def2-SV(P)"
        self.internal_defaults_orca["basis"]["default"] = "def2-TZVP(-f)"
        self.internal_defaults_orca["basis3"]["default"] = "def2-TZVP(-f)"
        self.internal_defaults_orca["basis_s"]["default"] = "def2-TZVP"
        self.internal_defaults_orca["basis_j"]["default"] = "def2-TZVP"

        # tm
        self.internal_defaults_tm = deepcopy(self.internal_defaults)
        self.internal_defaults_tm["sm2"]["default"] = "dcosmors"
        self.internal_defaults_tm["smgsolv1"]["default"] = "cosmors"
        self.internal_defaults_tm["smgsolv2"]["default"] = "cosmors"
        self.internal_defaults_tm["smgsolv3"]["default"] = "cosmors"
        self.internal_defaults_tm["sm4_j"]["default"] = "dcosmors"
        self.internal_defaults_tm["sm4_s"]["default"] = "dcosmors"
        # when basis == automatic and not "-3c"
        self.internal_defaults_tm["basis0"]["default"] = "def2-SV(P)"
        self.internal_defaults_tm["basis"]["default"] = "def2-TZVP"
        self.internal_defaults_tm["basis3"]["default"] = "def2-TZVP"
        self.internal_defaults_tm["basis_s"]["default"] = "def2-TZVP"
        self.internal_defaults_tm["basis_j"]["default"] = "def2-TZVP"
        self.internal_defaults_tm["basis_or"]["default"] = "def2-SVPD"

        for key, value in kwargs.items():
            if key == 'prog':
                update = list(self.internal_defaults.keys())
                for item in list(update):
                    if 'basis' in item:
                        # dont overwrite "automatic" 
                        update.remove(item)
                if value == 'tm':
                    tmp = {}
                    for key, value in self.internal_defaults_tm.items():
                        if key in update:
                            tmp[key] = value
                    self.internal_defaults.update(tmp)
                elif value == 'orca':
                    tmp = {}
                    for key, value in self.internal_defaults_orca.items():
                        if key in update:
                            tmp[key] = value
                    self.internal_defaults.update(tmp)

        self.value_options = {
            "nconf": ["all", "number e.g. 10 up to all conformers"],
            "charge": ["number e.g. 0"],
            "unpaired": ["number e.g. 0"],
            "solvent": ["gas"] + sorted([i for i in censo_solvent_db.keys()]),
            "prog": ["tm", "orca"],
            "prog2opt": ["tm", "orca", "prog", "automatic"],
            "part0": ["on", "off"],
            "part1": ["on", "off"],
            "part2": ["on", "off"],
            "part3": ["on", "off"],
            "part4": ["on", "off"],
            "optical_rotation": ["on", "off"],
            "prog3": ["tm", "orca", "prog"],
            "ancopt": ["on"],  # , "off"],
            "opt_spearman": ["on", "off"],
            "evaluate_rrho": ["on", "off"],
            "consider_sym": ["on", "off"],
            "prog_rrho": ["xtb",],
            "part0_gfnv": self.impgfnv,
            "part1_gfnv": self.impgfnv,
            "part2_gfnv": self.impgfnv,
            "part3_gfnv": self.impgfnv,
            "temperature": ["temperature in K e.g. 298.15"],
            "multitemp": ["on", "off"],
            "trange": ["temperature range [start, end, step]"],
            "func0": sorted(self.func_info.infos("func0", prog=None)),
            "basis0": ["automatic"]
            + sorted(
                list(self.func_info.composite_method_basis.values())
                + ["def2-SV(P)", "def2-TZVP"]
            ),
            "func": sorted(self.func_info.infos("func", prog=None)),
            "basis": ["automatic"]
            + sorted(
                list(self.func_info.composite_method_basis.values()) + ["def2-TZVP"]
            ),
            "func3": sorted(self.func_info.infos("func3", prog=None)),
            "basis3": sorted(knownbasissets),
            "part0_threshold": ["number e.g. 4.0"],
            "part1_threshold": ["number e.g. 5.0"],
            "opt_limit": ["number e.g. 4.0"],
            "part2_P_threshold": [
                "Boltzmann sum threshold in %. e.g. 95 (between 1 and 100)"
            ],
            "part3_threshold": [
                "Boltzmann sum threshold in %. e.g. 95 (between 1 and 100)"
            ],
            "sm2": sorted(self.impsm2),
            "smgsolv3": sorted(self.impsmgsolv3),
            "sm4_j": sorted(self.impsm4_j),
            "sm4_s": sorted(self.impsm4_s),
            "check": ["on", "off"],
            "crestcheck": ["on", "off"],
            "maxthreads": ["number of threads e.g. 2"],
            "omp": ["number cores per thread e.g. 4"],
            "smgsolv1": sorted(self.impsmgsolv1),
            "smgsolv2": sorted(self.impsmgsolv2),
            "bhess": ["on", "off"],
            "rmsdbias": ["on", "off"],
            "sm_rrho": ["alpb", "gbsa"],
            "optcycles": ["number e.g. 5 or 10"],
            "optlevel2": [
                "crude",
                "sloppy",
                "loose",
                "lax",
                "normal",
                "tight",
                "vtight",
                "extreme",
                "automatic",
            ],
            "spearmanthr": ["value between -1 and 1, if outside set automatically"],
            "couplings": ["on", "off"],
            "prog4_j": ["tm", "orca", "adf", "prog"],
            "prog4_s": ["tm", "orca", "adf", "prog"],
            "func_j": sorted(self.func_info.infos("func_j", prog=None)),
            "basis_j": sorted(knownbasissets),
            "func_s": sorted(self.func_info.infos("func_s", prog=None)),
            "basis_s": sorted(knownbasissets),
            "h_ref": sorted(self.imphref),
            "c_ref": self.impcref,
            "f_ref": self.impfref,
            "si_ref": self.impsiref,
            "p_ref": self.imppref,
            "h_active": ["on", "off"],
            "c_active": ["on", "off"],
            "f_active": ["on", "off"],
            "p_active": ["on", "off"],
            "si_active": ["on", "off"],
            "resonance_frequency": [
                "MHz number of your experimental spectrometer setup"
            ],
            "shieldings": ["on", "off"],
            "radsize": ["number e.g. 8 or 10"],
            "func_or": ["functional for opt_rot e.g. pbe"],
            "func_or_scf": ["functional for SCF in opt_rot e.g. r2scan-3c"],
            "basis_or": ["basis set for opt_rot e.g. def2-SVPD"],
            "freq_or": ["list of freq in nm to evaluate opt rot at e.g. [589, 700]"],
            "hlow": ["lowest force constant in ANC generation, e.g. 0.01"],
            "imagthr": ["automatic or e.g., -100    # in cm-1"],
            "sthr": ["automatic or e.g., 50     # in cm-1"],
            "scale": ["automatic or e.g., 1.0 "],
            "cosmorsparam": ["automatic"] + sorted(list(cosmors_param.keys())),
        }
        # must not be changed if restart(concerning optimization)
        self.restart_unchangeable = [
            "unpaired",
            "charge",
            "solvent",
            "prog",
            "prog2opt",
            "ancopt",
            "opt_spearman",
            "optlevel2",
            "func",
            "basis",
            "sm2",
            "nat",
            "radsize",
            "cosmorsparam",
        ]
        # might be changed, but data may be lost/overwritten
        self.restart_changeable = {
            "multitemp": False,
            # "temperature": False, # should not be changeable all solvent and
            # rrho values depend on this
            "trange": False,
            "bhess": False,
            "part1_gfnv": False,
            "part2_gfnv": False,
            "part3_gfnv": False,
            "smgsolv1": False,
            "smgsolv2": False,
            "smgsolv3": False,
            "func_or": False,
            "basis_or": False,
            "func_or_scf": False,
            "freq_or": False,
            "func3": False,
            "basis3": False,
            "func_j": False,
            "basis_j": False,
            "sm4_j": False,
            "func_s": False,
            "basis_s": False,
            "sm4_s": False,
            # "consider_sym": calculated on the fly
        }


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
        self.ensemblepath = ""
        self.configpath = ""
        self.jsonpath = ""

        # formatting:
        self.lenconfx = 3

        self.save_errors = []
        self.save_infos = []

        self.startupinfokeys = [
            "nat",
            "md5",
            "maxconf",
            "run",
            "configpath",
            "fixed_temperature",
            "vapor_pressure",
            "consider_unconverged",
        ]
        self.onlyread = False
        self.fixed_temperature = None
        self.nat = 0
        self.md5 = ""
        self.maxconf = 0
        self.run = True
        self.nmrmode = False
        self.consider_unconverged = False

        # pathsdefaults: --> read_program_paths
        self.external_paths = {}
        self.external_paths["orcapath"] = ""
        self.external_paths["orcaversion"] = ""
        self.external_paths["xtbpath"] = ""
        self.external_paths["crestpath"] = ""
        self.external_paths["cosmorssetup"] = ""
        self.external_paths["dbpath"] = ""
        self.external_paths["cosmothermversion"] = ""
        self.external_paths["mpshiftpath"] = ""
        self.external_paths["escfpath"] = ""

    def _set_fixed_temperature(self):
        """Initialize the temperature employed in part0 part1 and crude_rrho in the
        optimization (part2). This temperature is keept fixed for all restarts, to be able
        to compare the same (free) energies."""
        if self.fixed_temperature is None:
            setattr(self, "fixed_temperature", float(getattr(self, "temperature")))

    def cleanup_run(self, complete=False):
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
        """
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

    def check_logic(self, error_logical=False, silent=False):
        """
        Checks settings for impossible setting-comibinations, also checking
        if calculations are possible with the requested qm_codes.
        """
        if silent:
            store_errors = self.save_errors
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # check prog
        if self.prog not in self.value_options["prog"]:
            self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Prog can only be chosen from "
                    f"{self.value_options['prog']} and not --> {self.prog}!"
                )
            error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle prog2opt:
        if self.prog2opt in ("prog", "automatic"):
            if self.prog2opt in self.value_options["prog2opt"]:
                self.prog2opt = self.prog
            else:
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Prog2opt can only be chosen from "
                    "('tm', 'orca', 'prog')!"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle prog3:
        if self.prog3 == "prog" and self.prog in self.value_options["prog"]:
            self.prog3 = self.prog
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle prog4_j:
        if self.prog4_j == "prog" and self.prog in self.value_options["prog"]:
            self.prog4_j = self.prog
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle prog4:
        if self.prog4_s == "prog" and self.prog in self.value_options["prog"]:
            self.prog4_s = self.prog
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # set spearmanthr by number of atoms:
        if self.spearmanthr < -1 or self.spearmanthr > 1:
            self.spearmanthr = 1 / (math.exp(0.03 * (self.nat ** (1 / 4))))
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle prog_rrho
        if self.prog_rrho == "prog" and self.prog in self.value_options["prog"]:
            self.prog_rrho = self.prog
            if self.prog_rrho == "tm":
                if shutil.which("thermo") is not None:
                    # need thermo for reading thermostatistical contribution
                    self.prog_rrho = "tm"
                else:
                    self.prog_rrho = "xtb"
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Currently are only GFNn-xTB "
                        "hessians possible and no TM hessians"
                    )
        elif not self.prog_rrho:
            self.save_errors.append(
                f"{'WARNING:':{WARNLEN}}Thermostatistical contribution to "
                "free energy will not be calculated, since prog_rrho ist set to 'off'!"
            )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # COSMO-RS only works with TM
        if self.solvent != "gas":
            if (
                self.smgsolv1 in ("cosmors", "cosmors-fine")
                and self.prog == "orca"
                and self.part1
            ):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Part1 when used with COSMO-RS can only be prog = TM!\n"
                    f"{'':{WARNLEN}}Set prog to tm or change smgsolv1 to e.g. smd_gsolv!"
                )
                error_logical = True
            if (
                self.smgsolv2 in ("cosmors", "cosmors-fine")
                and self.prog == "orca"
                and self.part2
            ):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Part2 when used with COSMO-RS can only be prog = TM!\n"
                    f"{'':{WARNLEN}}Set prog to tm or change smgsolv2 to e.g. smd_gsolv!"
                )
                error_logical = True
            if (
                self.smgsolv3 in ("cosmors", "cosmors-fine")
                and self.prog3 == "orca"
                and self.part3
            ):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Part3 when used with COSMO-RS can only be prog3 = TM!\n"
                    f"{'':{WARNLEN}}Set prog3 to tm or change smgsolv3 to e.g. smd_gsolv!"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # adjust functional-names to new naming convention (cfg.functional)
        self.func0 = self.func_info.relay_functionals.get(self.func0, self.func0)
        self.func = self.func_info.relay_functionals.get(self.func, self.func)
        self.func3 = self.func_info.relay_functionals.get(self.func3, self.func3)
        self.func_j = self.func_info.relay_functionals.get(self.func_j, self.func_j)
        self.func_s = self.func_info.relay_functionals.get(self.func_s, self.func_s)
        self.func_or = self.func_info.relay_functionals.get(self.func_or, self.func_or)
        self.func_or_scf = self.func_info.relay_functionals.get(
            self.func_or_scf, self.func_or_scf
        )
        # extracheck for r2scan-3c and ORCA
        try:
            orcaversion = int(self.external_paths['orcaversion'].split('.')[0])
        except (ValueError, AttributeError):
            orcaversion = 2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func0
        # available with QM code
        if self.part0:
            if self.func0 not in self.func_info.infos("func0", prog=self.prog):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"(func0) {self.func0} is not implemented with the {self.prog} program package.\n"
                )
                error_logical = True
            if (self.func0 == "r2scan-3c" and
                self.prog == "orca" and
                orcaversion < 5):
                self.save_errors.append(f"{'WARNING:':{WARNLEN}}The composite "+
                "functional r2scan-3c is only available in ORCA since version 5.0.0 !")
            # check if available for part0
            if self.func0 not in self.func_info.infos("func0"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"{self.func0} is not available for part0.\n"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle basis0 for func0:
        if self.part0:
            if (
                self.basis0 == "None"
                or self.basis0 is None
                or self.basis0 == "automatic"
            ):
                if self.prog == "tm":
                    default = self.internal_defaults_tm.get(
                        "basis0", {"default": "def2-SV(P)"}
                    )["default"]
                elif self.prog == "orca":
                    default = self.internal_defaults_orca.get(
                        "basis0", {"default": "def2-SV(P)"}
                    )["default"]
                else:
                    default = "def2-SV(P)"
                self.basis0 = self.func_info.composite_method_basis.get(
                    self.func0, default
                )
            if self.basis0 not in knownbasissets:
                if self.basis0 == "def2-TZVP(-f)" and self.prog == "tm":
                    self.basis0 = "def2-TZVP"
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basis0: {self.basis0} "
                        f"is used instead of {'def2-TZVP(-f)'}"
                    )
                else:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basis0: {self.basis0} could not be "
                        + "checked, but is used anyway!"
                    )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func
        # now stuff gets fun:
        if self.prog != self.prog2opt:
            if self.part1 or self.part2:
                # check if func is available in both QM codes
                if self.func not in self.func_info.infos("func"):
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The functional "
                        f"{self.func} is not available for part1 / part2.\n"
                    )
                    error_logical = True
                if (self.func not in self.func_info.infos("func", prog=self.prog) or
                    self.func not in self.func_info.infos("func", prog=self.prog2opt)):
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The functional "
                        f"{self.func} is not available in both QM codes prog= {self.prog}"
                        f" and prog2opt= {self.prog2opt} !\n"
                        f"{'':{WARNLEN}}You can choose a functional for 'func' from:"
                    )
                    self.save_errors.append(print_block(sorted(list(set(self.func_info.infos("func", prog=self.prog)) &
                        set(self.func_info.infos("func", prog=self.prog2opt)))), redirect=True)
                    )
                    error_logical = True
                elif (self.func in self.func_info.infos("func", prog=self.prog) and
                    self.func in self.func_info.infos("func", prog=self.prog2opt)):
                    self.save_errors.append(
                        f"{'INFORMATION:':{WARNLEN}}The functional "
                        f"{self.func} is available in both QM codes prog= {self.prog}"
                        f" and prog2opt= {self.prog2opt} !\n"
                        f"{'':{WARNLEN}}Using func '{dfa_settings.functionals.get(self.func).get(self.prog)}' with prog {self.prog} and\n"
                        f"{'':{WARNLEN}}using func '{dfa_settings.functionals.get(self.func).get(self.prog2opt)}' with prog2opt {self.prog2opt} ."
                    )
                if self.func == "r2scan-3c" and orcaversion < 5:
                    self.save_errors.append(f"{'WARNING:':{WARNLEN}}The composite "+
                    "functional r2scan-3c is only available in ORCA since version 5.0.0 !"
                    )
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Handle basis for func:
            if self.part1 or self.part2:
                if self.basis == "None" or self.basis is None or self.basis == "automatic":
                    if self.prog == "tm":
                        default = self.internal_defaults_tm.get(
                            "basis", {"default": "def2-TZVP"}
                        )["default"]
                    elif self.prog == "orca":
                        default = self.internal_defaults_orca.get(
                            "basis", {"default": "def2-TZVP"}
                        )["default"]
                    else:
                        default = "def2-TZVP"
                    self.basis = self.func_info.composite_method_basis.get(
                        self.func, default
                    )
                    if self.basis == "def2-TZVP(-f)":
                        self.basis = "def2-TZVP"
                        self.save_errors.append(
                            f"{'INFORMATION:':{WARNLEN}}The basis set basis: {self.basis} "
                            f"is used instead of {'def2-TZVP(-f)'}"
                        )
                if self.basis not in knownbasissets:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basis: {self.basis} could not be "
                        + "checked, but is used anyway!"
                    )
        else:
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Handle func
            # available with QM code
            if self.part1 or self.part2:
                if self.func not in self.func_info.infos("func", prog=self.prog):
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The functional "
                        f"(func) {self.func} is not implemented with the {self.prog} program package.\n"
                    )
                    error_logical = True
                if (self.func == "r2scan-3c" and
                    self.prog == "orca" and
                    orcaversion < 5):
                    self.save_errors.append(f"{'WARNING:':{WARNLEN}}The composite "+
                    "functional r2scan-3c is only available in ORCA since version 5.0.0 !")
                # check if available for part1 / part2
                if self.func not in self.func_info.infos("func"):
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The functional "
                        f"{self.func} is not available for part1 / part2.\n"
                    )
                    error_logical = True
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Handle basis for func:
            if self.part1 or self.part2:
                if self.basis == "None" or self.basis is None or self.basis == "automatic":
                    if self.prog == "tm":
                        default = self.internal_defaults_tm.get(
                            "basis", {"default": "def2-TZVP"}
                        )["default"]
                    elif self.prog == "orca":
                        default = self.internal_defaults_orca.get(
                            "basis", {"default": "def2-TZVP"}
                        )["default"]
                    else:
                        default = "def2-TZVP"
                    self.basis = self.func_info.composite_method_basis.get(
                        self.func, default
                    )
                if self.basis not in knownbasissets:
                    if self.basis == "def2-TZVP(-f)" and self.prog == "tm":
                        self.basis = "def2-TZVP"
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}The basis set basis: {self.basis} "
                            f"is used instead of {'def2-TZVP(-f)'}"
                        )
                    else:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}The basis set basis: {self.basis} could not be "
                            + "checked, but is used anyway!"
                        )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func3 dsd-blyp with basis
        if (
            self.part3
            and self.func3 in ("dsd-blyp", "dsd-blyp-d3")
            and self.basis3 != "def2-TZVPP"
        ):
            self.save_errors.append(
                f"{'WARNING:':{WARNLEN}}DSD-BLYP is only available with the basis set def2-TZVPP!"
            )
            self.basis3 = "def2-TZVPP"
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle basis3:
        if self.part3:
            if (
                self.basis3 == "None"
                or self.basis3 is None
                or self.basis3 == "automatic"
            ):
                if self.prog3 == "tm":
                    default = self.internal_defaults_tm.get(
                        "basis3", {"default": "def2-TZVP"}
                    )["default"]
                elif self.prog3 == "orca":
                    default = self.internal_defaults_orca.get(
                        "basis3", {"default": "def2-TZVP"}
                    )["default"]
                else:
                    default = "def2-TZVP"
                self.basis3 = self.func_info.composite_method_basis.get(
                    self.func3, default
                )
            if self.basis3 not in knownbasissets:
                if self.basis3 == "def2-TZVP(-f)" and self.prog3 == "tm":
                    self.basis3 = "def2-TZVP"
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basis3: {self.basis3} "
                        f"is used instead of {'def2-TZVP(-f)'}"
                    )
                else:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basis3: {self.basis3} could not be "
                        + "checked, but is used anyway!"
                    )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func3
        if self.part3:
            if self.func3 in self.func_info.composite_method_basis.keys():
                if self.basis3 != self.func_info.composite_method_basis[self.func3]:
                    self.save_errors.append(
                        f"{'INFORMATION:':{WARNLEN}}You are using a basis "
                        f"set (basis3) {self.basis3} different to the original composite method"
                        f" basis set ({self.func_info.composite_method_basis[self.func3]})!"
                    )
            # available with QM code
            if self.func3 not in self.func_info.infos("func3", prog=self.prog3):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"(func3) {self.func3} is not implemented with the {self.prog3} program package.\n"
                )
                error_logical = True
            if (self.func3 == "r2scan-3c" and 
                self.prog3 == "orca" and 
                orcaversion< 5):
                self.save_errors.append(f"{'WARNING:':{WARNLEN}}The composite "+
                "functional r2scan-3c is only available in ORCA since version 5.0.0 !")
            # check if available for part3
            if self.func3 not in self.func_info.infos("func3"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"{self.func3} is not available for part3.\n"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.part4 and (self.couplings or self.shieldings):
            self.nmrmode = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func_j
        if self.part4:
            # available with QM code
            if self.func_j not in self.func_info.infos("func_j", prog=self.prog4_j):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"(funcJ) {self.func_j} is not implemented with the {self.prog4_j} "
                    "program package.\n"
                )
                error_logical = True
            if (self.func_j == "r2scan-3c" and 
                self.prog4_j == "orca" and 
                orcaversion < 5):
                self.save_errors.append(f"{'WARNING:':{WARNLEN}}The composite "+
                "functional r2scan-3c is only available in ORCA since version 5.0.0 !")
            # check if available for part4
            if self.func_j not in self.func_info.infos("func_j"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"{self.func_j} is not available for part4.\n"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle basis_j:
        if self.part4:
            if (
                getattr(self, "basis_j", None) is None
                or getattr(self, "basis_j", None) == "automatic"
            ):
                default_basis_j = self.func_info.composite_method_basis.get(
                    self.func_j, "def2-TZVP"
                )
                setattr(self, "basis_j", default_basis_j)
            if self.basis_j not in knownbasissets:
                if self.basis_j == "def2-TZVP(-f)" and self.prog4_j == "tm":
                    self.basis_j = "def2-TZVP"
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basisJ: {self.basis_j} "
                        f"is used instead of {'def2-TZVP(-f)'}"
                    )
                else:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basisJ: {self.basis_j} "
                        "could not be checked, but is used anyway!"
                    )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func_s
        if self.part4:
            # available with QM code
            if self.func_s not in self.func_info.infos("func_s", prog=self.prog4_s):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"(funcS) {self.func_s} is not implemented with the {self.prog4_s} "
                    "program package.\n"
                )
                error_logical = True
            if (self.func_s == "r2scan-3c" and 
                self.prog4_s == "orca" and 
                orcaversion < 5):
                self.save_errors.append(f"{'WARNING:':{WARNLEN}}The composite "+
                "functional r2scan-3c is only available in ORCA since version 5.0.0 !")
            # check if available for part4
            if self.func_s not in self.func_info.infos("func_s"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"{self.func_s} is not available for part4.\n"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle basis_s:
        if self.part4:
            if (
                getattr(self, "basis_s", None) is None
                or getattr(self, "basis_s", None) == "automatic"
            ):
                default_basis_s = self.func_info.composite_method_basis.get(
                    self.func_s, "def2-TZVP"
                )
                setattr(self, "basis_s", default_basis_s)
            if self.basis_s not in knownbasissets:
                if self.basis_s == "def2-TZVP(-f)" and self.prog4_s == "tm":
                    self.basis_s = "def2-TZVP"
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basisS: {self.basis_s} "
                        f"is used instead of {'def2-TZVP(-f)'}"
                    )
                else:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basisS: {self.basis_s} could not be "
                        + "checked, but is used anyway!"
                    )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # no unpaired electrons in coupling or shielding calculations!
        if self.unpaired > 0:
            if self.part4 and (self.couplings or self.shieldings):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}Coupling and shift calculations "
                    "(part4) are only available for closed-shell systems!"
                )
                error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Solvation:
        if self.solvent == "gas":
            self.smgsolv1 = "gas-phase"
            self.sm2 = "gas-phase"
            self.smgsolv2 = "gas-phase"
            self.smgsolv3 = "gas-phase"
            self.sm4_j = "gas-phase"
            self.sm4_s = "gas-phase"
        else:
            if self.vapor_pressure:
                self.save_errors.append(
                    f"{'INFORMATION:':{WARNLEN}}The vapor_pressure flag only affects "
                    f"settings for COSMO-RS!\n"
                    f"{'':{WARNLEN}}Information on solvents with properties similar "
                    f"to the input molecule must be provided for other solvent models!"
                )
            if self.part2:
                # Handle sm2 --> solvent model in optimization:
                exchange_sm = {
                    "cosmo": "cpcm",
                    "cpcm": "cosmo",
                    "dcosmors": "smd",
                    "smd": "dcosmors",
                }
                if self.sm2 not in self.impsm2:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The solvent model {self.sm2} is not implemented!"
                    )
                    error_logical = True
                if self.prog2opt == "orca":
                    if self.sm2 in self.sm2_tm:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm2} is not available with "
                            f"{self.prog2opt}! Therefore {exchange_sm[self.sm2]} is used!"
                        )
                        self.sm2 = exchange_sm[self.sm2]
                    elif self.sm2 == "default":
                        self.sm2 = self.internal_defaults_orca["sm2"]["default"]
                if self.prog2opt == "tm":
                    if self.sm2 in self.sm2_orca:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm2} is not available with "
                            f"{self.prog2opt}! Therefore { exchange_sm[self.sm2]} is used!"
                        )
                        self.sm2 = exchange_sm[self.sm2]
                    elif self.sm2 == "default":
                        self.sm2 = self.internal_defaults_tm["sm2"]["default"]
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if self.part1:
                # Handle smgsolv1
                exchange_sm = {
                    "cosmo": "cpcm",
                    "cpcm": "cosmo",
                    "dcosmors": "smd",
                    "smd": "dcosmors",
                }
                if self.smgsolv1 not in self.impsmgsolv1:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The solvent model {self.smgsolv1}"
                        " is not implemented for smgsolv1 !"
                    )
                    error_logical = True
                if self.smgsolv1 == "sm2":
                    self.smgsolv1 = self.sm2
                if self.prog == "tm" and self.smgsolv1 in self.sm2_orca:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}{self.smgsolv1} is not available with "
                        f"{self.prog}! Therefore {exchange_sm[self.smgsolv1]} is used!"
                    )
                    self.smgsolv1 = exchange_sm[self.smgsolv1]
                if self.prog == "orca" and self.smgsolv1 in self.sm2_tm:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}{self.smgsolv1} is not available with "
                        f"{self.prog}! Therefore {exchange_sm[self.smgsolv1]} is used!"
                    )
                    self.smgsolv1 = exchange_sm[self.smgsolv1]
                if (
                    self.smgsolv1 in ("alpb_gsolv", "gbsa_gsolv", "smd_gsolv")
                    and self.multitemp
                ):
                    self.save_errors.append(
                        f"{'INFORMATION:':{WARNLEN}}{self.smgsolv1} does not provide "
                        "information at different temperatures!"
                    )
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if self.part2:
                # Handle smgsolv2
                exchange_sm = {
                    "cosmo": "cpcm",
                    "cpcm": "cosmo",
                    "dcosmors": "smd",
                    "smd": "dcosmors",
                }
                if self.smgsolv2 not in self.impsmgsolv2:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The solvent model {self.smgsolv2}"
                        " is not implemented for smgsolv2 !"
                    )
                    error_logical = True
                if self.smgsolv2 == "sm2":
                    self.smgsolv2 = self.sm2
                if self.prog == "tm" and self.smgsolv2 in self.sm2_orca:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}{self.smgsolv2} is not available with "
                        f"{self.prog}! Therefore {exchange_sm[self.smgsolv2]} is used!"
                    )
                    self.smgsolv2 = exchange_sm[self.smgsolv2]
                if self.prog == "orca" and self.smgsolv2 in self.sm2_tm:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}{self.smgsolv2} is not available with "
                        f"{self.prog}! Therefore {exchange_sm[self.smgsolv2]} is used!"
                    )
                    self.smgsolv2 = exchange_sm[self.smgsolv2]
                if (
                    self.smgsolv2 in ("alpb_gsolv", "gbsa_gsolv", "smd_gsolv")
                    and self.multitemp
                ):
                    self.save_errors.append(
                        f"{'INFORMATION:':{WARNLEN}}{self.smgsolv2} does not provide "
                        "information at different temperatures!"
                    )
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if self.part3:
                # Handle smgsolv3
                exchange_sm = {
                    "cosmo": "cpcm",
                    "cpcm": "cosmo",
                    "dcosmors": "smd",
                    "smd": "dcosmors",
                }
                if self.smgsolv3 not in self.impsmgsolv3:
                    self.save_errors.append(
                        f"{'ERROR:':{WARNLEN}}The solvent model {self.smgsolv3}"
                        " is not implemented for smgsolv3 !"
                    )
                    error_logical = True
                if self.smgsolv3 == "sm2":
                    self.smgsolv3 = self.sm2
                if self.prog3 == "tm" and self.smgsolv3 in self.sm2_orca:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}{self.smgsolv3} is not available with "
                        f"{self.prog3}! Therefore {exchange_sm[self.smgsolv3]} is used!"
                    )
                    self.smgsolv3 = exchange_sm[self.smgsolv3]
                if self.prog3 == "orca" and self.smgsolv3 in self.sm2_tm:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}{self.smgsolv3} is not available with "
                        f"{self.prog3}! Therefore {exchange_sm[self.smgsolv3]} is used!"
                    )
                    self.smgsolv3 = exchange_sm[self.smgsolv3]
                if (
                    self.smgsolv3 in ("alpb_gsolv", "gbsa_gsolv", "smd_gsolv")
                    and self.multitemp
                ):
                    self.save_errors.append(
                        f"{'INFORMATION:':{WARNLEN}}{self.smgsolv3} does not provide "
                        "information at different temperatures!"
                    )
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if self.part4:
                # Handle sm4_j
                if self.prog4_j == "orca":
                    if self.sm4_j in self.sm4_j_tm:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_j} is not available with {self.prog4_j}!"
                            f" Therefore {exchange_sm[self.sm4_j]} is used!"
                        )
                        self.sm4_j = exchange_sm[self.sm4_j]
                    elif self.sm4_j == "default":
                        self.sm4_j = self.internal_defaults_orca["sm4_j"]["default"]
                if self.prog4_j == "tm":
                    if self.sm4_j in self.sm4_j_orca:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_j} is not available with {self.prog4_j}!"
                            f" Therefore {exchange_sm[self.sm4_j]} is used!"
                        )
                        self.sm4_j = exchange_sm[self.sm4_j]
                    elif self.sm4_j == "default":
                        self.sm4_j = self.internal_defaults_tm["sm4_j"]["default"]
                # Handle sm4_s
                if self.prog4_s == "orca":
                    if self.sm4_s in self.sm4_s_tm:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_s} is not available with {self.prog4_s}!"
                            f" Therefore { exchange_sm[self.sm4_s]} is used!"
                        )
                        self.sm4_s = exchange_sm[self.sm4_s]
                    elif self.sm4_s == "default":
                        self.sm4_s = self.internal_defaults_orca["sm4_s"]["default"]
                if self.prog4_s == "tm":
                    if self.sm4_s in self.sm4_s_orca:
                        self.save_errors.append(
                            f"{'WARNING:':{WARNLEN}}{self.sm4_s} is not available with {self.prog4_s}!"
                            f" Therefore {exchange_sm[self.sm4_s]} is used!"
                        )
                        self.sm4_s = exchange_sm[self.sm4_s]
                    elif self.sm4_s == "default":
                        self.sm4_s = self.internal_defaults_tm["sm4_s"]["default"]
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Check if solvent-information is available for solventmodel
            # Check which solvation models are applied:
            check_for = {
                "xtb": False,
                "cosmors": False,
                "dcosmors": False,
                "cpcm": False,
                "smd": False,
                "DC": False,
            }
            applied_solventmodels = []
            if self.evaluate_rrho:
                applied_solventmodels.append(self.sm_rrho)
            if self.part1:
                applied_solventmodels.append(self.smgsolv1)
            if self.part2:
                applied_solventmodels.append(self.sm2)
                applied_solventmodels.append(self.smgsolv2)
            if self.part3:
                applied_solventmodels.append(self.smgsolv3)
            if self.part4:
                applied_solventmodels.append(self.sm4_j)
                applied_solventmodels.append(self.sm4_s)
            if self.optical_rotation:
                applied_solventmodels.append("cosmo")

            for solventmodel in list(set(applied_solventmodels)):
                if solventmodel in ("alpb", "gbsa", "alpb_gsolv", "gbsa_gsolv"):
                    check_for["xtb"] = True
                elif solventmodel in ("cosmors", "cosmors-fine"):
                    check_for["cosmors"] = True
                elif solventmodel in ("dcosmors",):
                    check_for["dcosmors"] = True
                    check_for["DC"] = True
                elif solventmodel in ("cosmo",):
                    check_for["DC"] = True
                elif solventmodel in ("cpcm",):
                    check_for["cpcm"] = True
                elif solventmodel in ("smd", "smd_gsolv"):
                    check_for["smd"] = True
                else:
                    print("unexpected behaviour solvents")
            lookup = {
                "xtb": "solvents_xtb",
                "cosmors": "solvents_cosmors",
                "dcosmors": "solvents_dcosmors",
                "cpcm": "solvents_cpcm",
                "smd": "solvents_smd",
                "DC": "",
            }
            # check if solvent in censo_solvent_db
            if censo_solvent_db.get(self.solvent, "not_found") == "not_found":
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The solvent '{self.solvent}'' is not found in your file "
                    f"{os.path.expanduser(os.path.join('~/.censo_assets/', 'censo_solvents.json'))}!"
                    f"\n{'':{WARNLEN}}Check your input!"
                )
                error_logical = True
            else:
                for key, value in check_for.items():
                    if value:
                        if (
                            censo_solvent_db[self.solvent].get(key, "nothing_found")
                            == "nothing_found"
                        ):
                            self.save_errors.append(
                                f"{'ERROR:':{WARNLEN}}The solvent for solventmodel in "
                                "{key} is not found!"
                            )
                            error_logical = True
                        if key == "DC":
                            try:
                                if censo_solvent_db[self.solvent].get(key, None) is not None:
                                    _ = float(censo_solvent_db[self.solvent].get(key, None))
                                else:
                                    self.save_errors.append(
                                        f"{'ERROR:':{WARNLEN}}The dielectric constant for the solvent '{self.solvent}' "
                                        f"is not provided for the solventmodel {'cosmo / dcosmors'}!"
                                    )
                                    error_logical = True
                            except ValueError:
                                self.save_errors.append(
                                    f"{'ERROR:':{WARNLEN}}The dielectric constant can "
                                    "not be converted."
                                )
                                error_logical = True
                        elif key in ("smd", "cpcm"):
                            if (censo_solvent_db[self.solvent].get(key, ['', 'nothing_found'])[1].lower() 
                                not in getattr(self, lookup[key]) and
                                censo_solvent_db[self.solvent].get(key, ['', None])[1]):
                                self.save_errors.append(
                                    f"{'WARNING:':{WARNLEN}}The solvent "
                                    f"'{censo_solvent_db[self.solvent].get(key, 'nothing_found')[1]}'"
                                    f" for solventmodel/program {key} can not be checked "
                                    "but is used anyway."
                                )
                        else:
                            if (censo_solvent_db[self.solvent].get(key, ['', 'nothing_found'])[1]
                                not in getattr(self, lookup[key]) and
                                censo_solvent_db[self.solvent].get(key, ['', None])[1]):
                                self.save_errors.append(
                                    f"{'WARNING:':{WARNLEN}}The solvent "
                                    f"'{censo_solvent_db[self.solvent].get(key, 'nothing_found')[1]}' "
                                    f"for solventmodel/program {key} can not be checked "
                                    "but is used anyway."
                                )
                        if (key != "DC" and 
                            censo_solvent_db[self.solvent].get(key, ["",""])[0] is None and
                            censo_solvent_db[self.solvent].get(key, "nothing_found") != "nothing_found"
                            ):
                            if key == 'xtb':
                                tmp_sm = 'alpb'
                            else:
                                tmp_sm = key
                            if censo_solvent_db[self.solvent].get(key, ["",None])[1] is None:
                                self.save_errors.append(
                                    f"{'ERROR:':{WARNLEN}}The solvent '{self.solvent}' "
                                    f"is not parameterized for the solventmodel {tmp_sm}, and "
                                    f"no replacement is available!!!"
                                )
                                error_logical = True
                            else:
                                self.save_errors.append(
                                    f"{'WARNING:':{WARNLEN}}The solvent '{self.solvent}' "
                                    f"is not parameterized for the solventmodel {tmp_sm}, therefore"
                                    f" '{censo_solvent_db[self.solvent].get(key, ['',''])[1]}' is used!!!"
                                )

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle optlevel2:
        # sm2 needs to be set (not default!)
        if self.optlevel2 in ("None", None, "automatic"):
            if self.sm2 in ("smd", "dcosmors") and self.solvent != "gas":
                self.optlevel2 = "lax"
            else:
                # gas phase
                self.optlevel2 = "normal"
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.part4 and not self.couplings and not self.shieldings:
            self.part4 = False
            self.save_errors.append(
                f"{'INFORMATION:':{WARNLEN}}Neither coupling nor "
                "shielding constants are activated! Part4 is not executed."
            )
        elif self.part4 and not any(
            [
                getattr(self, flag)
                for flag in (
                    "h_active",
                    "c_active",
                    "f_active",
                    "si_active",
                    "p_active",
                )
            ]
        ):
            self.save_errors.append(
                f"{'INFORMATION:':{WARNLEN}}No active element for the calculation of NMR is "
                "activated in the .censorc! Therefore all nuclei are calculated!"
            )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # optical rotation
        if self.optical_rotation and self.prog == "orca":
            self.save_errors.append(
                f"{'ERROR:':{WARNLEN}}Optical rotation calculations are only possible with TM."
            )
            error_logical = True
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle func_or
        if self.optical_rotation:
            # available with QM code
            if self.func_or not in self.func_info.infos("func_or", prog="tm"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"(funcOR) {self.func_or} is not implemented with the "
                    "{'tm'} program package.\n"
                )
                error_logical = True
            # check if available for part5
            if self.func_or not in self.func_info.infos("func_or"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"{self.func_or} is not available for part5.\n"
                )
                error_logical = True
            # available with QM code
            if self.func_or_scf not in self.func_info.infos("func_or_scf", prog="tm"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"(funcOR) {self.func_or_scf} is not implemented with the "
                    "{'tm'} program package.\n"
                )
                error_logical = True
            # check if available for part5
            if self.func_or_scf not in self.func_info.infos("func_or_scf"):
                self.save_errors.append(
                    f"{'ERROR:':{WARNLEN}}The functional "
                    f"{self.func_or_scf} is not available for part5.\n"
                )
                error_logical = True
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Handle basis_or:
            if (
                getattr(self, "basis_or", None) is None
                or getattr(self, "basis_or", None) == "automatic"
            ):
                default_basis_or = self.func_info.composite_method_basis.get(
                    self.func_or_scf, "def2-SVPD"
                )
                setattr(self, "basis_or", default_basis_or)
            if self.basis_or not in knownbasissets:
                if self.basis_or == "def2-TZVP(-f)" and self.prog == "tm":
                    self.basis_or = "def2-TZVP"
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basisOR: {self.basis_or} "
                        f"is used instead of {'def2-TZVP(-f)'}"
                    )
                else:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}The basis set basisOR: {self.basis_or} "
                        "could not be checked, but is used anyway!"
                    )

        if silent:
            self.save_errors = store_errors
        return error_logical

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

    def read_program_paths(self, configpath, silent=False):
        """
        Get absolute paths of external programs employed in censo
        Read from the configuration file .censorc
        """
        if silent:
            keep = self.save_errors

        with open(configpath, "r") as inp:
            stor = inp.readlines()
        for line in stor:
            if "ctd =" in line:
                try:
                    self.external_paths["cosmorssetup"] = str(line.rstrip(os.linesep))
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read settings for COSMO-RS from .censorc!"
                    )
                try:
                    normal = "DATABASE-COSMO/BP-TZVP-COSMO"
                    fine = "DATABASE-COSMO/BP-TZVPD-FINE"
                    tmp_path = self.external_paths["cosmorssetup"].split()[5].strip('"')
                    if "OLDPARAM" in tmp_path:
                        tmp_path = os.path.split(tmp_path)[0]
                    tmp_path = os.path.split(tmp_path)[0]
                    self.external_paths["dbpath"] = tmp_path
                    self.external_paths["dbpath_fine"] = os.path.join(tmp_path, fine)
                    self.external_paths["dbpath_normal"] = os.path.join(
                        tmp_path, normal
                    )
                except Exception as e:
                    self.save_errors.append(e)
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read settings for COSMO-RS from "
                        f".censorc!\n{'':{WARNLEN}}Most probably there is a user "
                        "input error."
                    )
            if "ORCA:" in line:
                try:
                    self.external_paths["orcapath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for ORCA from .censorc!."
                    )
            if "ORCA version:" in line:
                try:
                    tmp = line.split()[2]
                    tmp = tmp.split(".")
                    tmp.insert(1, ".")
                    tmp = "".join(tmp)
                    self.external_paths["orcaversion"] = tmp
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read ORCA version from .censorc!"
                    )
            if "GFN-xTB:" in line:
                try:
                    self.external_paths["xtbpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for GFNn-xTB from .censorc!"
                    )
                    if shutil.which("xtb") is not None:
                        self.external_paths["xtbpath"] = shutil.which("xtb")
                        self.save_errors.append(
                            f"{'':{WARNLEN}}Going to use {self.external_paths['xtbpath']} instead."
                        )
            if "CREST:" in line:
                try:
                    self.external_paths["crestpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for CREST from .censorc!"
                    )
                    if shutil.which("crest") is not None:
                        self.external_paths["crestpath"] = shutil.which("crest")
                        self.save_errors.append(
                            f"{'':{WARNLEN}}Going to use {self.external_paths['crestpath']} instead."
                        )
            if "mpshift:" in line:
                try:
                    self.external_paths["mpshiftpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for mpshift from .censorc!"
                    )
            if "escf:" in line:
                try:
                    self.external_paths["escfpath"] = str(line.split()[1])
                except Exception:
                    self.save_errors.append(
                        f"{'WARNING:':{WARNLEN}}Could not read path for escf from .censorc!"
                    )
            if "$ENDPROGRAMS" in line:
                break

            if silent:
                self.save_errors = keep

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

            if os.path.exists(path_to_cefine):
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

    def read_json(self, path, silent=False):
        """
        Reading stored data on conformers and information on settings of
        previous run.
        """
        if os.path.isfile(path):
            if not silent:
                print("Reading file: {}\n".format(os.path.basename(path)))
            try:
                with open(path, "r", encoding=CODING, newline=None) as inp:
                    save_data = json.load(inp, object_pairs_hook=OrderedDict)
            except (ValueError, TypeError, FileNotFoundError):
                print(f"{'ERROR:':{WARNLEN}}Your Jsonfile (enso.json) is corrupted!\n")
                time.sleep(0.02)
                raise
        return save_data

    def write_json(self, path, conformers, settings, outfile="enso.json"):
        """
        Dump conformer data and settings information of current run to json file
        """
        data = {}
        if not isinstance(settings, OrderedDict):
            data["settings"] = vars(settings)
        else:
            data["settings"] = settings
        try:
            conformers.sort(key=lambda x: int(x["id"]))
        except Exception:
            pass
        for conf in conformers:
            if not isinstance(conf, OrderedDict):
                data[conf.id] = vars(conf)
            else:
                data[conf["id"]] = conf
        with open(os.path.join(path, outfile), "w") as out:
            json.dump(data, out, indent=4, sort_keys=False)
