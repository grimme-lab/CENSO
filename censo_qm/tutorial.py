"""
Guide for setting up CENSO and performing a CENSO calculation
"""
from .cfg import DESCR
def interactiv_doc():
    """Guide for interactive explaination of CENSO"""
    
    everything = """
    Commandline Energetic Sorting (CENSO) is a sorting algorithm for efficient
    evaluation of Structure Ensembles (SE). The input ensemble (or single structure)
    originating from a CREST[SQM/FF] run can be ranked by free energy at DFT level
    and structures can be optimized at DFT level. CENSO has a modular structure, 
    enabling efficient sorting at different levels of theory and targeting individual 
    error sources, e.g., energy or solvation. Sorting is based on (free) windows or
    thresholds, within which all conformers/structures are considered. By choosing
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

    Parts 0-3 are concerned with efficent SE sorting, optimization and calculation
    of Boltzmann weigths for populated structures. The parts are described in the 
    following:

    Part0 - Cheap prescreening:

    Flexible and/or large molecules can have many conformers (i.e., several hundred)
    and sorting out truely high lying conformers fast is crucial for efficiency. 
    This is the goal of part0. Here the electronic energy description is improved upon
    the initial SQM/FF energy by performing very fast B97-D3(0)/def2-SV(P)+gcp 
    single-point calculations. If the molecule is in solution phase, solvation is treated
    at GFN2-xTB[ALPB] level. Sorting is based on g_thr,0 wich has to be rather large, 
    e.g. 4 kcal/mol or above. 

    Part1 - Prescreening:

    Accurate electronic and solvation energies are calculated in part1. To be efficient
    COSMO-RS solvation contributions are calculated at r2SCAN-3c level, whereby
    both the gas phase electronic energy and solvation contribution is obtained and 
    no additional calculations are necessary. Conformers above the threshold g_thr,1
    are discarded. After sorting, thermostatistical contributions (G_mRRHO) are 
    calculated at GFN2-xTB[ALPB] level using single-point hessian (SPH) calculations. 
    Full free energies are calculated and sorting is performed based on G_thr,1. 
    The threshold can be automatically increased (fuzzy-threshold) as a function of
    the standard deviation of G_mRRHO, which is the case for flexible or diverse SE.
    Up to now all calculations have been performed on the SQM/FF input geometries!


    Part2 - Optimiation:

    The relevant conformers/structures have to be optimized at DFT (r2SCAN-3c) level 
    in implicit solvation (DCOSMO-RS). An efficient ensemble optimizer has been 
    implemented, where all conformers are optimized for 8 iterations. Then a spearman 
    correlation coefficient is calculated to check for parallel potential energy surface.
    In case of paralellity can be assumed, the soring threshold G_thr,2 is decreased
    and conformers above the threshold are discarded, if their gradient norm is below
    a predefined threshold. For large ensembles, this decreases the not populated 
    conformers fast. The batch wise optimization is repeated until all confomers within
    the energy window G_thr,2 are fully optimized. On the DFT optimized geometries
    free energies are calculated like in part1, with COSMO-RS(r2SCAN-3c) for E and
    dG_solv and GFN2-xTB[ALPB]-SPH for G_mRRHO. Boltzmann weights are calculated and
    an ensemble averaged free energy can be obtained. The conformers which are further
    considered for use in parts 3-5 are considered up to Boltzmann sum threshold, e.g.,
    all conformers up to 90 % population are further considered.


    Part3 - Refinement:

    If it is necessary to refine the calculated Boltzmann weights at a higher (hybrid)
    DFT level, this can be performed in part3. Also it has to be noted, that Boltzmann
    weights calculated from rather accurate r2SCAN-3c energies are reliable enough 
    for most applications. Like in part2 free energy contributions are calculated 
    on DFT optimized geometries, although on a higher (hybrid) DFT level. Conformers
    below a Boltzmann sum threshold are considered further (e.g. in part4 or part5).


    Part4 - NMR-Mode:


    work in progress ...



    """
    # options = [
    #             "everything",
    #             "setting_up",
    #             "part0",
    #             "part1",
    #             "part2",
    #             "part3",
    #             "part4",
    #             "thresholds"
    # ]
    options = ["everything",]
    print(DESCR)
    print("This is the CENSO tutorial / interactive documentation:\n")
    print("Topic options are:")
    for option in options:
        print(f"  {option}")
    print("\nTo exit please type one of the following: exit or q")
    print("\nPlease input your information request:")
    while True:
        user_input = input()
        if user_input.strip() in ('q', 'exit'):
            break
        if user_input.strip() not in options:
            print(f"Options are: {options}")
        else:
            print(f"You requested information on {user_input}")
            if user_input == 'everything':
                print(everything)
            break

    # setup cosmors spaces are necessary 