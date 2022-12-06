"""
Guide for setting up CENSO and performing a CENSO calculation
"""
from .cfg import DESCR, censo_solvent_db, dfa_settings
from .inputhandling import internal_settings

def make_block(strlist, width=80):
    """Print all elements of strlist in block mode
    e.g. within 80 characters then newline
    - width [int] width of block
    """
    length = 4
    stor = []
    try:
        maxlen = max([len(str(x)) for x in strlist])
    except (ValueError, TypeError):
        maxlen = 12
    tmp = "    "
    for item in strlist:
        length += maxlen + 2
        if length <= width:
            tmp += f"{item}, "
        else:
            tmp += f"{item},\n"
            stor.append(tmp)
            tmp="    "
            length=4
    if tmp not in stor:
        stor.append(tmp)
    return stor



def interactiv_doc():
    """Guide for interactive explaination of CENSO"""

    default_settings = internal_settings()
    func_info = dfa_settings()

    general = """
    Commandline Energetic Sorting (CENSO) is a sorting algorithm for efficient
    evaluation of Structure Ensembles (SE). The input ensemble (or single structure)
    originating from a CREST[SQM/FF] run can be ranked by free energy at DFT level
    and structures can be optimized at DFT level. CENSO has a modular structure, 
    enabling efficient sorting at different levels of theory and targeting individual 
    error sources, e.g., energy or solvation. Sorting is based on (free) energy windows
    or thresholds, within which all conformers/structures are considered. By choosing
    appropriate thresholds, the following "properties" can be obtained:

    * lowest lying conformer and corresponding free energy
    * fully optimized SE and averaged ensemble free energy
    * Boltzmann populated SE at a given temperature

    Apart from the refinement/sorting of SE, CENSO offers an automated approach to 
    calculate NMR properties (shielding and coupling constants) or optical rotatory 
    (OR) dispersion.

    To perform fast calculations CENSO interfaces with QM codes like xTB, TURBOMOLE 
    and ORCA. The results are read by CENSO and evaluated in the sorting algorithm.

    CENSO is structured into several parts, which can all be turned on or off by the 
    user. These parts are presented in their intended order: 

    1) Part0: Cheap prescreening
    2) Part1: Prescreening
    3) part2: Optimization (and free energy calculation)
    4) part3: Refinement   (optional at higher hybrid-DFA level)
    ------------------------------------------------------------
    The following parts are "property" parts for NMR or OR and not necessary if only 
    a SE refinement is needed.
    5) part4: NMR-Mode
    6) part5: OR-Mode


    Nomenclature: Free energies are not available in each sorting step, due to the 
    computational cost involved. Sorting-thresholds based on incomplete free energies
    e.g., without thermostatistical contribution (G_mRRRHO) are denoted by lower case
    g_thr and full free energies by upper case G_thr.

    Parts 0-3 are concerned with efficient SE sorting, optimization and calculation
    of Boltzmann weigths for populated structures. The parts are described in the 
    following:

    Part0 - Cheap prescreening:

    Flexible and/or large molecules can have many conformers (i.e., several hundred)
    and sorting out truly high lying conformers fast is crucial for efficiency. 
    This is the goal of part0. Here the electronic energy description is improved upon
    the initial SQM/FF energy by performing very fast B97-D3(0)/def2-SV(P)+gcp 
    single-point calculations. If the molecule is in solution phase, solvation is treated
    at GFN2-xTB[ALPB] level. Sorting is based on g_thr(0) which has to be rather large, 
    e.g. 4 kcal/mol or above. 

    Part1 - Prescreening:

    Accurate electronic and solvation energies are calculated in part1. To be efficient
    COSMO-RS solvation contributions are calculated at r2SCAN-3c level, whereby
    both the gas phase electronic energy and solvation contribution is obtained and 
    no additional calculations are necessary. Conformers above the threshold g_thr(1)
    are discarded. After sorting, thermostatistical contributions (G_mRRHO) are 
    calculated at GFN2-xTB[ALPB] level using single-point hessian (SPH) calculations. 
    Full free energies are calculated and sorting is performed based on G_thr(1). 
    The threshold can be automatically increased (fuzzy-threshold) as a function of
    the standard deviation of G_mRRHO, which is the case for flexible or diverse SE.
    Up to now all calculations have been performed on the SQM/FF input geometries!


    Part2 - Optimization:

    The relevant conformers/structures have to be optimized at DFT (r2SCAN-3c) level 
    in implicit solvation (DCOSMO-RS). An efficient ensemble optimizer has been 
    implemented, where all conformers are optimized for 8 iterations. Then a spearman 
    correlation coefficient is calculated to check for parallel potential energy surfaces.
    If parallellity can be assumed, the sorting threshold G_thr(2) is decreased
    and conformers above the threshold are discarded, if their gradient norm is below
    a predefined threshold. For large ensembles, this decreases the number of high 
    lying conformers fast. The batch wise optimization is repeated until all confomers 
    within the energy window G_thr(2) are fully optimized. On the DFT optimized geometries
    free energies are calculated like in part1, with COSMO-RS(r2SCAN-3c) for E and
    dG_solv and GFN2-xTB[ALPB]-SPH for G_mRRHO. Boltzmann weights are calculated and
    an ensemble averaged free energy can be obtained. The conformers for use in 
    parts 3-5 are considered up to a Boltzmann sum threshold, e.g., all conformers 
    up to 90 % population are further considered.


    Part3 - Refinement:

    If it is necessary to refine the calculated Boltzmann weights at a higher (hybrid)
    DFT level, this can be performed in part3. Also it has to be noted, that Boltzmann
    weights calculated from rather accurate r2SCAN-3c energies are reliable enough 
    for most applications. Like in part2 free energy contributions are calculated 
    on DFT optimized geometries, although using a higher (hybrid) DFT level. Conformers
    below a Boltzmann sum threshold are considered further (e.g. in part4 or part5).


    Optional property related parts (part4 and part5) are described below:

    Part4 - NMR-Mode:
    
    In part4 NMR properties can be calculated for the populated conformers. The 
    Boltzmann weights are taken from either part1, part2 or part3 if they are available.
    All populated conformers up to the population part2_P_threshold or part3_threshold
    are considered. Coupling and shielding constants are calculated separately 
    and can be calculated for the elements H, C, F, Si, P, or all elements.
    Files for further processing with ANMR are created.
    After the CENSO run, NMR spectra can be calculated using the ANMR code.

    Part5 - OR-Mode:

    In part5 optical rotatory (OR) dispersion of the populated structure ensemble
    can be calculated. Boltzmann weights can be taken from part1, part2 or part3.

    """

    censorc = """

    Set up CENSO:

    When running CENSO for the first time you are asked to creat a global 
    configuration file .censorc. This file contains all settings CENSO operates 
    with and you should adjust the settings to your needs. In the first part of 
    the .censorc the paths of the employed programs have to be provided. The .censorc
    can be read by CENSO, if the file is located in your home folder e.g. ~/.censorc
    or it is located in your current working directory. The current directory is 
    always preferred to the global censorc in your home directory.

    Definition of keywords in the censorc:

    $GENERAL SETTINGS:

    nconf           = how many conformers should be considered. Either a number  
                      or the flag 'all'.
    charge          = molecular charge of the molecule under investigation.
    unpaired        = number of unpaired electrons in the molecule under 
                      investigation.
    solvent         = Solvent if the molecule is in solution phase, else 'gas'.
    prog_rrho       = QM-code used for the calculation of thermostatistical 
                      contributions. This is only feasible with xtb, since normally 
                      a large number of hessian calculations have to be performed.
    temperature     = Temperature (in Kelvin) used for the Boltzmann evaluations.
    trange          = temperature range which is used to calculate free energies 
                      at different temperatures (considered in G_mRRHO and 
                      dG_solv[only COSMO-RS]). The temperature range will only be
                      evaluated if multitemp is set to 'on'.
    multitemp       = Evaluate free energies at different temperatures defined 
                      in trange.
    evaluate_rrho   = Option to consider /not consider thermostatistical contributions.
    consider_sym    = Option to consider symmetry in the thermostatistical contribution
                      (only xtb).
    bhess           = Calculate single point hessians on the input geometry,
                      instead of "ohess" (optimization + hessian calculation).
    imagthr         = threshold for inverting imaginary frequencies for thermostatistical
                      contributions (in cm-1). Internal defaults are applied if set to 
                      'automatic'.
    sthr            = rotor cut-off (cm-1) used for the thermostatical contributions.
                      Internal defaults are applied if set to 'automatic'.
    scale:          = scaling factor for frequencies in vibrational partition function.
                      Internal defaults are applied if set to 'automatic'.
    rmsdbias        = gESC related, using rmsdpot.xyz to be consistent to CREST
    sm_rrho         = solvent model applied in the GFNn-xTB thermostatistical 
                      contribution calculation.
    check           = Terminate the CENSO run if too many calculations crash.
    prog            = QM code used for part0, part1 and part2, this can be 
                      TURBOMOLE or ORCA.
    func            = functional used in part1 (prescreening) and part2 (optimization)
    basis           = basis set used in combination with func in part1 
                      (prescreening) and part2 (optimization). If basis is set to
                      'automatic' the basis set is chosen internally.
    maxthreads      = Used for parallel calculation. Maxthreads determines the 
                      number of independent calculations running in parallel.
                      E.g. resulting in 4 independent single-point /optimization
                      calculations.
    omp             = Used for parallel calculation. Omp determines the number of 
                      cores each independent calculation can use. Eg. maxthreads = 4
                      and omp = 5 resulting in 4 independent calculations and each
                      independent calculation uses 5 cores.
    cosmorsparam    = Flag for choosing COSMO-RS parameterizations. If set to
                      'automatic' the input from the COSMO-RS input line is choosen.

    $PART0 - CHEAP-PRESCREENING - SETTINGS:  

    part0           = Option to turn the "cheap prescreening part" on or off.
    func0           = functional used in part0
    basis0          = Basis set used in combination with func0. If basis0 is set to
                      'automatic' the basis set is choosen internally.
    part0_gfnv      = GFN version employed in the thermostatistical contribution 
                      in part0. 
    part0_threshold = Threshold/Energy-window (kcal/mol) within which all conformers
                      are considered.

    $PART1 - PRESCREENING - SETTINGS:

    part1:          = Option to turn the "prescreening part" on or off.
    smgsolv1:       = Additive solvation contribution employed in part1. 
    part1_gfnv:     = GFN version employed in the thermostatistical contribution 
                      in part1.
    part1_threshold = Threshold/Energy-window (kcal/mol) within which all conformers
                      are considered further.

    $PART2 - OPTIMIZATION - SETTINGS:

    part2           = Option to turn the "optimization part" on or off.
    prog2opt        = QM program used only for the optimization in part2.
                      Either 'tm' or 'orca'. It is implemented mainly for the following use case
                      optimization of polar molecules with ORCA/SMD instead of TM/DCOSMO-RS and 
                      still be able to calculated GSolv and the Electronic energy with TM.
    part2_threshold = Threshold/Energy-window (kcal/mol) within which all conformers
                      are fully optimized.
    sm2             = Implicit solvation model used in the optimization.
    smgsolv2        = Additive solvation model used for calculation of dG_solv
                      in part2.
    part2_gfnv      = GFN version employed in the thermostatistical contribution 
                      in part2.
    ancopt          = Using ANCoptimizer implemented in xTB for geometry optimization.
    hlow            = Lowest force constant in ANC generation, used with ancopt.
    opt_spearman    = Using the new "ensemble-optimizer".
    part2_P_threshold = Boltzmann threshold in '%' within which all conformers
                      are considered further. E.g. 90 --> all conformers up to 
                      a sum of 90 '%' are considered.
    optlevel2       = Optimization threshold in the geometry optimization. If set to
                      "automatic", internal default will be used.
    optcycles       = Number of optimization iterations performed in the ensemble
                      optimizer.
    spearmanthr     = Spearman rank correlation coeff. used to determine if PES 
                      during geometry optimization can be assumed parallel. 
    radsize         = Setting of the radial grid size for func used in part2.
    crestcheck      = Automatically sort out conformers which might have become
                      identical or rotamers during geometry optimization, using 
                      CREST (this is threshold based, so use with care).

    $PART3 - REFINEMENT - SETTINGS:
    part3           = Option to turn the "refinement part" on or off.
    prog3           = QM code used for part3 this can be TURBOMOLE or ORCA.
    func3           = functional used in part3 (refinement)
    basis3          = basis set employed in combination with func3.
    smgsolv3        = Additive solvation model used for calculation of dG_solv in 
                      part3.
    part3_gfnv      = GFN version employed in the thermostatistical contribution 
                      in part3.
    part3_threshold = Boltzmann threshold in '%' within which all conformers
                      are considered further. E.g. 90 --> all conformers up to 
                      a sum of 90 '%' are considered.

    $PART4 SETTINGS:
    part4           = Option to turn the "NMR property part" on or off.
    couplings       = Perform coupling constant calculations [options are on or off].
    progJ           = QM code (TM, ORCA) used for coupling constant calculations.
    funcJ           = Density functional employed for the coupling constant calculation.
    basisJ          = basis set employed with the DFA (funcJ) for coupling 
                      constant calculations.
    sm4J            = implicit solvent model employed in the coupling constant calculation.
    shieldings      = Perform shielding constant calculations [options are on or off].
    progS           = QM code (TM, ORCA) used for shielding constant calculations.
    funcS           = Density functional employed for the shielding constant calculation.
    basisS          = basis set employed with the DFA (funcS) for shielding 
                      constant calculations.
    sm4S            = implicit solvent model employed in the shielding constant calculation.
    reference_1H    = Reference molecule to convert 1H shieldings to shifts e.g. TMS.
    reference_13C   = Reference molecule to convert 13C shieldings to shifts e.g. TMS.
    reference_19F   = Reference molecule to convert 19F shieldings to shifts e.g. CFCl3.
    reference_29Si  = Reference molecule to convert 29Si shieldings to shifts e.g. TMS.
    reference_31P   = Reference molecule to convert 31P shieldings to shifts e.g. TMP.
    1H_active       = Calculate 1H NMR properties [options are on or off].
    13C_active      = Calculate 13C NMR properties [options are on or off].
    19F_active      = Calculate 19F NMR properties [options are on or off].
    29Si_active     = Calculate 29Si NMR properties [options are on or off].
    31P_active      = Calculate 31P NMR properties [options are on or off].
    resonance_frequency = Resonance frequency of the experimental spectrometer.

    $OPTICAL ROTATION PROPERTY SETTINGS:
    $PART5 SETTINGS:
    optical_rotation        = Option to turn the "OR property part" on or off.
    funcOR:                 = Functional employed to calculate the optical rotatory (OR) dispersion.
    funcOR_SCF: r2scan-3c   = Functional to generate converged MOs.
    basisOR: def2-SVPD      = Basis set employed for the OR calculation.
    frequency_optical_rotv  = List of frequencies in nm to evaluate OR at e.g. [589.0].

    """
    setup = """

    How to perform a calculation:

    If not done before, create the global configuration file .censorc.

    $ censo -newconfig

    Edit the configuration file, e.g. adjust thresholds and provide paths to 
    employed programs.

    Check if all input combinations go along by running:

    $ censo -inp <inputfile.xyz> -part0 on -part1 on -part2 off -checkinput

    And then run the calculation:

    $ censo -inp <inputfile.xyz> -part0 on -part1 on -part2 off > censo.out &

    When the calculation succeded and additional information is wanted, 
    e.g. increase the number of considered conformers or perform DFT 
    geometry optimizations (part2) then this can be accieved with the keyword:
    '-restart' and additionally the new request:

    $ censo -restart -part2 on > censo.out &


"""
    files = """
    # Files written by censo

    CENSO generates several files. Most important is the file enso.json which 
    contains all conformer/ensemble and setting information and is required for 
    restarting calculations.

    The file coord.enso_best contains the geometry of the lowest lying conformer
    determined by the sorting procedure. The sorted ensemble geometries are written to
    the files enso_ensemble_part*.xyz.

    The file coord.enso_best is (over-) written in parts 0-3 and contains the geometry of
    the conformer with the lowest free energy. 

    The files enso_ensemble_part1.xyz and enso_ensemble_part2.xyz contain the sorted 
    ensembles after part1 or part2. The files enso_ensemble_part2_p_XX.xyz contain
    the sorted ensemble up to a Boltzmann population of XX %.

    The files partX.dat are just for convenience and contain the same printout 
    of free energies as in the complete censo output.

    The files .anmrrc and anmr_enso are written in part4 (NMR) and are needed for
    further processing with the ANMR code.

    """
    jobscript = """
    Example jobscript:

    It is important to source the correct QM programs e.g. Turbomole or ORCA otherwise
    CENSO might not find them or simply use any version which is found in your PATH.

    #!/bin/bash
    # PBS Job
    #PBS -V
    #PBS -N JOB_NAME
    #PBS -m ae
    #PBS -q batch
    #PBS -l nodes=1:ppn=28
    # 
    cd $PBS_O_WORKDIR

    ### setup programs
    ## XTB
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    ulimit -s unlimited
    export OMP_STACKSIZE=1000m

    ## TM
    export PARA_ARCH=SMP
    source /home/$USER/bin/.turbo751
    export PARNODES=4  ## omp 
    export TM_PAR_FORK=1

    ### ORCA4.2.1
    ORCAPATH="/tmp1/orca_4_2_1_linux_x86-64_openmpi216";
    MPIPATH="/software/openmpi-2.1.5/bin";
    MPILIB="/software/openmpi-2.1.5/lib64";
    PATH=${ORCAPATH}:${MPIPATH}:$PATH 
    LD_LIBRARY_PATH=${ORCAPATH}:${MPILIB}:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=/software/intel/parallel_studio_xe_2017.1/parallel_studio_xe_2017.4.056/compilers_and_libraries_2017/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=/software/intel/parallel_studio_xe_2017/mkl/lib/intel64:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH

    ## PATH
    PATH=/home/abt-grimme/AK-bin:$PATH
    PATH=/home/$USER/bin:$PATH
    export PATH
    ### end programs + PATH

    export HOSTS_FILE=$PBS_NODEFILE
    cat $HOSTS_FILE>hosts_file

    TMP_DIR=/tmp1/$USER
    DIR1=$PWD

    mkdir -p $TMP_DIR/$PBS_JOBID

    #check file system access
    if [ ! -d $TMP_DIR/$PBS_JOBID ]; then
     echo "Unable to create $TMP_DIR/$PBS_JOBID  on $HOSTNAME. Must stop."
     exit
    fi

    #check current location
    if [ "$PWD" == "$HOME" ]; then
     echo "Cowardly refusing to copy the whole home directory"
     exit
    fi

    #copy everything to node (will NOT copy directories for safety reasons.
    #Add an 'r' only if absolutely sure what you are doing)
    #bwlimit limits bandwidth to 5000 kbytes/sec

     rsync -q --bwlimit=5000 $DIR1/* $TMP_DIR/$PBS_JOBID/
     rsync -rq --ignore-missing-args --bwlimit=5000 $DIR1/CONF* $TMP_DIR/$PBS_JOBID/
     rsync -rq --ignore-missing-args --bwlimit=5000 $DIR1/NMR* $TMP_DIR/$PBS_JOBID/
     rsync -q --bwlimit=5000 $DIR1/.* $TMP_DIR/$PBS_JOBID/
     cd $TMP_DIR/$PBS_JOBID

    ####################################################################################
    #Gettimings
    start=$(date +%s)
    #####################################################################################
    #jobs start here (if you have no idea what this script does, only edit this part...)

    echo "Calculation from $(date)" >> RUNTIME
    export PYTHONUNBUFFERED=1

    censo -inp inputfile.xyz -P 7 -O 4 > censo.out

    #end of job      (....and stop editing here.)
    #####################################################################################
    #Print timings to file
    end=$(date +%s)
    secs=$(expr $end - $start)
    printf '%dh:%dm:%02ds\\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60)) >> RUNTIME
    #####################################################################################
    #copy everything back that is smaller than 5 gbytes

     rsync -rq --bwlimit=5000 --max-size=5G $TMP_DIR/$PBS_JOBID/* $DIR1/
     rsync -q --bwlimit=5000 --max-size=5G $TMP_DIR/$PBS_JOBID/.* $DIR1/

    #to be safe, get mos alpha and beta seperately. 
    #Note that the rsync syntax is strange; you need to first include everything, 
    #then exclude the rest ("*\" includes subdirectories)

     rsync -rq --bwlimit=5000 --include="*/" --include="mos" --include="alpha" --include="beta" --exclude=* $TMP_DIR/$PBS_JOBID/* $DIR1/

    #if you want the large files as well, comment in the following

    #rsync -r --bwlimit=1000 --min-size=5G $TMP_DIR/$PBS_JOBID/* $DIR1/

     cd $DIR1
     rm -r $TMP_DIR/$PBS_JOBID

    """

    thresholds = f"""
    CENSO is a threshold based sorting algorithm. The thresholds have been determined
    for typical drug-like organic molecules up to 100 atoms.

    For cases with high flexibility or very different conformer ordering between SQM/FF and
    DFT it can be necessary to increase the thresholds.

    If your target quantity is an averaged ensemble free energy, be sure not to 
    use small energy windows in part0 and part1, as this can falsely reduce the 
    ensemble size.

    In the following the employed thresholds are listed and their use is explained:

    Part0 - Cheap Prescreening:
            
            part0_threshold = g_thr(0)
            cml:      -part0_threshold
            censorc:  part0_threshold

            part0_threshold g_thr(0) is an energy window /threshold in kcal/mol 
            within which all conformers are considered. The CENSO internal default 
            for g_thr(0) = {default_settings.internal_defaults.get("part0_threshold",{})["default"]} kcal/mol. The sorting threshold is designed to 
            remove conformers very high lying in electronic energy (E).


    Part1 - Prescreening:

            In part1 two sorting steps are applied. One based on g_thr(1) with 
            improved energy and solvation description and one on G_thr(1) where 
            thermostatistical contributions have been calculated additionally. 
            

            part1_threshold = g_thr(1) and G_thr(1)
            cml:      -part1_threshold 
            censorc:  part1_threshold

            part1_threshold g_thr(1) is an energy window /threshold in kcal/mol 
            within which all conformers are considered. The CENSO internal default 
            for g_thr(1) = {default_settings.internal_defaults.get("part1_threshold",{})["default"]} kcal/mol. The sorting threshold g_thr(1) is designed to 
            remove high lying conformers based on improved electronic energy (E)
            and solvation contributions (dG_solv). G_thr(1) sorting is based on
            full free energy including thermostatistical contributions (G_mRRHO).
            Both, g_thr(1) and G_thr(1) use the same base threshold e.g. 3.5 kcal/mol
            and G_thr(1) can be automatically increased as a function of the standard
            deviation of G_mRRHO in the structure ensemble, indicating high structural
            diversity and possibly larger errors (fuzzy sorting).


    Part2 - Optimization:

            In part2 two sorting thresholds are applied. One which is applied during
            geometry optimization and one Boltzmann sum threshold.

            threshold applied during optimization G_thr(opt,2):
            cml:      -thrpart2
            censorc:  part2_threshold

            Spearman-threshold for testing for parallel PES:
            cml:      -spearmanthr
            censorc:  spearmanthr

            Boltzmann sum threshold G_thr(2):
            cml:      -part2_Pthreshold
            censorc:  part2_P_threshold

            The internal default for G_thr(opt,2) is set to {default_settings.internal_defaults.get("opt_limit",{})["default"]} kcal/mol. During 
            the geometry optimization the initial threshold G_thr(opt,2) is
            increased by 60 % (or at least 1.5 kcal/mol). If the PES can be assumed
            parallel (tested by a spearman rank coefficient) the threshold is 
            decresed until it reaches the initial value. All conformers with (free)
            energies above the G_thr(opt,2) threshold are discarded and all conformers
            below G_thr(opt,2) will by fully DFT optimized. The last threshold
            employed in part2 G_thr(2) is a Boltzmann sum threshold and all populated
            conformers of the ensemble, up to the Boltzmann population (in %)
            are considered further. 

    Part3 - Refinement:

            In part3 a Boltzmann sum threshold is employed.
            
            Boltzmann sum threshold G_thr(3):
            cml:      -thrpart3
            censorc:  part3_threshold
            
            Based on high level free energies Boltzmann weights are calculated 
            and all conformers up to the Boltzmann sum threshold (in %) are 
            considered further.

    Part4 - NMR properties:

            In part4 only Boltzmann sum thresholds of part2 G_thr(2) or part3 G_thr(3) 
            are applied. It is possible to calculate part4 using Boltzmann weights
            from part1, but since there is no Boltzmann sum threshold in part1, in this 
            case, the threshold G_thr(2) is used. The NMR properties are calculated 
            on populated conformers up to the Boltzmann threshold in %. 

    Part5 - Optical rotatory dispersion:

            In part4 only Boltzmann sum thresholds of part2 G_thr(2) or part3 G_thr(3) 
            are applied. It is possible to calculate part5 using Boltzmann weights
            from part1, but since there is no Boltzmann sum threshold in part1, in this 
            case, the threshold G_thr(2) is used. The OR property is calculated 
            on populated conformers up to the Boltzmann threshold in %. For optical
            rotation it is necessary to include almost the entire ensemble e.g. 99 %
    """

    solvation = """
    CENSO uses several QM-packages and not all solvents are available for all 
    solvation models throughout the QM-packages.
    For this reason a user editable file is created in the folder:

        $  ~/.censo_assets/censo_solvents.json

    which contains a dictionary of all available solvent models and solvents.
    If a solvent is not available with a certain solvent model, the user can then
    choose a replacement solvent. E.g. if benzene is not available choose toluene. 

    Example for the solvent name in CENSO, the solvent models (xtb represents ALPB or GBSA).
    The list [None, "toluene"] means that the solvent is not found in the solvent
    model and that the solvent toluene in the same solvent model is used instead.
    DC = Dieclectric constant (used in COSMO).

    "acetone": {
        "cosmors": ["propanone_c0", "propanone_c0"], 
        "dcosmors": ["propanone", "propanone"],
        "xtb": ["acetone", "acetone"],
        "cpcm": ["acetone", "acetone"],
        "smd": ["ACETONE", "ACETONE"],
        "DC": 20.7,
    },
    "benzene": {
        "cosmors": ["benzene_c0", "benzene_c0"],
        "dcosmors": [None, "toluene"],
        "xtb": ["benzene", "benzene"],
        "cpcm": ["Benzene", "Benzene"],
        "smd": ["BENZENE", "BENZENE"],
        "DC": 2.3,
    },

    The solvent file is directly used in `CENSO` and typos will cause the 
    calculations to crash! Adding a new solvent is as easy as adding a new 
    dictionary entry to the file.


    In CENSO several solvent models can be applied. The intent is either a good 
    description of the free energy (keyword: smgsolv*) or an implicit effect on
    a property or geometry (keyword: sm*).

    (SM) implicit solvation for properties:
    * COSMO         [TM]
    * CPCM          [ORCA]
    * DCOSMO-RS     [TM]
    * ALPB          [xtb]
    * GBSA          [xtb]
    * SMD           [ORCA]

    (SMGSOLV) implicit solvation for free energies:

    * COSMO-RS      [COSMO-RS]
    * SMD_Gsolv     [ORCA]
    * ALPB_Gsolv    [xtb]
    * GBSA_Gsolv    [xtb]

    Available solvents and naming convention employed in CENSO:
    
"""
    tmp = "    "
    for count, solvent in enumerate(list(censo_solvent_db.keys())):
        tmp += f"{solvent}, "
        if (count + 1) % 5 == 0:
            solvation += tmp + "\n"
            tmp = "    "
    solvation += tmp

    example_applications = """
    CENSO can be used for several applications / target quantities. Some are 
    listed below:

    For the demonstration purpose it is assumed that all parts are turned 'off'
    in the global configuration file of the user!


    * Calculate fast DFT (B97-D3(0)/def2-SV(P)+gcp) single-point energies on 
      GFNn-xTB geometries:

      $ censo -inp ensemble.xyz -part0 on  > censo.out &

    * Calculate part1 (cheap prescreening) in parallel, with 7 independent threads
      using each 4 cores:

      $ censo -inp ensemble.xyz -part0 on -P 7 -O 4 > censo.out &

    * Calculate free energies in solution phase (CHCl3) on GFNn-xTB geometries:

      $ censo -inp ensemble.xyz -part1 on -solvent chcl3 -func r2scan-3c -basis automatic > censo.out &

    * Reduce the ensemble to relevant populated conformers and optimize the 
      SQM/FF ensemble at r2scan-3c level:

      $ censo -inp ensemble.xyz -part0 on -part1 on -part2 on -func r2scan-3c -basis automatic > censo.out &

    * Restart on a previous calculation and calculate optical rotation:

      $ censo -restart -OR on > censo.out &

    """

    functionals = f"""
    Functionals that can be employed by TM for func0:\n
{''.join(make_block(func_info.infos('func0', prog='tm'), width=120))}

    Functionals that can be employed by TM for func:\n
{''.join(make_block(func_info.infos('func', prog='tm'), width=120))}

    Functionals that can be employed by TM for func3:\n
{''.join(make_block(func_info.infos('func3', prog='tm'), width=120))}

    Functionals that can be employed by TM for funcJ:\n
{''.join(make_block(func_info.infos('func_j', prog='tm'), width=120))}

    Functionals that can be employed by TM for funcS:\n
{''.join(make_block(func_info.infos('func_s', prog='tm'), width=120))}

    Functionals that can be employed by TM for funcOR:\n
{''.join(make_block(func_info.infos('func_or', prog='tm'), width=120))}

    Functionals that can be employed by TM for funcOR_SCF:\n
{''.join(make_block(func_info.infos('func_or_scf', prog='tm'), width=120))}

    Functionals that can be employed by ORCA for func0:\n
{''.join(make_block(func_info.infos('func0', prog='orca'), width=120))}

    Functionals that can be employed by ORCA for func:\n
{''.join(make_block(func_info.infos('func', prog='orca'), width=120))}

    Functionals that can be employed by ORCA for func3:\n
{''.join(make_block(func_info.infos('func3', prog='orca'), width=120))}

    Functionals that can be employed by ORCA for funcJ:\n
{''.join(make_block(func_info.infos('func_j', prog='orca'), width=120))}

    Functionals that can be employed by ORCA for funcS:\n
{''.join(make_block(func_info.infos('func_s', prog='orca'), width=120))}
    """

    # CONFORMER numbering kept from crest input

    # folders which are created
    # part0_sp
    # GFN_unbiased/
    # rrho_part1/
    # b97-3c/ folder of func name

    # which functionals are available

    tutorial_data = {
        "general": general,
        "censorc": censorc,
        "setup": setup,
        "thresholds": thresholds,
        "solvation": solvation,
        "examples": example_applications,
        "files": files,
        "functionals": functionals,
        "jobscript": jobscript,
    }
    tutorial_data["everything"] = "\n\n".join(tutorial_data.values())

    options = list(tutorial_data.keys())
    print(DESCR)
    print("This is the CENSO tutorial / interactive documentation:\n")
    print("Topic options are:")
    for option in options:
        print(f"    {option}")
    print("\nTo exit please type one of the following: exit or q")
    print("\nPlease input your information request:")
    while True:
        user_input = input()
        if user_input.strip() in ("q", "exit"):
            break
        if user_input.strip() not in options:
            print(f"Options are: {options}")
        else:
            print(f"You requested information on {user_input}:\n")
            print(tutorial_data.get(user_input.strip(), "Not found."))
            # if user_input == 'everything':
            #     print(everything)
            # if user_input == 'jobscript':
            #     print(jobscript)
            # break
            print("\n\nDo you want information on any of the other topics?")
            print("\nIf you want to exit please type one of the following: exit or q")
    print("\n****CENSO TUTORIAL END****")
