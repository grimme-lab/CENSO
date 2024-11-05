"""
Storing constants for the use in all CENSO modules.
Storing program paths --> still in transition
Storing censo_solvent_db solvent database across all solvation models (as fallback)
"""

import os
import sys

from .__version__ import __version__


class Config:
    ENVIRON = os.environ.copy()
    if getattr(sys, "frozen", False):  # if bundled by pyinstaller ...
        # workaround for LD_LIBRARY_PATH pyinstaller and suse
        LP_KEY = "LD_LIBRARY_PATH"  # for GNU/Linux and *BSD.
        lp_orig = ENVIRON.get(LP_KEY + "_ORIG")
        if lp_orig is not None:
            ENVIRON[LP_KEY] = lp_orig  # restore the original, unmodified value
        else:
            # This happens when LD_LIBRARY_PATH was not set.
            # Remove the env var as a last resort:
            ENVIRON.pop(LP_KEY, None)
        # end workaround

    CODING = "ISO-8859-1"

    ASSETS_PATH = __file__.replace("params.py", "assets")

    USER_ASSETS_PATH = os.path.join(os.path.expanduser("~"), ".censo2_assets")

    PROGS = ("orca", "tm")

    SOLV_MODS: dict[str, tuple] = {
        "orca": ("cpcm", "smd"),
        "tm": ("cosmo", "dcosmors", "cosmors", "cosmors-fine"),
        "xtb": ("alpb", "gbsa"),
    }

    GRIDOPTIONS = (
        "low",
        "low+",
        "high",
        "high+",
    )

    GFNOPTIONS = (
        "gfnff",
        "gfn1",
        "gfn2",
    )

    CENSORCNAME = ".censo2rc"

    OMPMIN = 4

    OMPMAX = 32

    OMP = OMPMIN

    NCORES = os.cpu_count()

    COSMORS_PARAM = {
        "12-normal": "BP_TZVP_C30_1201.ctd",
        "13-normal": "BP_TZVP_C30_1301.ctd",
        "14-normal": "BP_TZVP_C30_1401.ctd",
        "15-normal": "BP_TZVP_C30_1501.ctd",
        "16-normal": "BP_TZVP_C30_1601.ctd",
        "17-normal": "BP_TZVP_C30_1701.ctd",
        "18-normal": "BP_TZVP_18.ctd",
        "19-normal": "BP_TZVP_19.ctd",
        "12-fine": "BP_TZVPD_FINE_HB2012_C30_1201.ctd",
        "13-fine": "BP_TZVPD_FINE_HB2012_C30_1301.ctd",
        "14-fine": "BP_TZVPD_FINE_C30_1401.ctd",
        "15-fine": "BP_TZVPD_FINE_C30_1501.ctd",
        "16-fine": "BP_TZVPD_FINE_C30_1601.ctd",
        "17-fine": "BP_TZVPD_FINE_C30_1701.ctd",
        "18-fine": "BP_TZVPD_FINE_18.ctd",
        "19-fine": "BP_TZVPD_FINE_19.ctd",
    }


DESCR = f"""
         ______________________________________________________________
        |                                                              |
        |                                                              |
        |                   CENSO - Commandline ENSO                   |
        |{'v ' + __version__:^{62}}|
        |    energetic sorting of CREST Conformer Rotamer Ensembles    |
        |                    University of Bonn, MCTC                  |
        |                           Oct 2024                           |
        |                 based on ENSO version 2.0.1                  |
        |             L. M. Seidler, F. Bohle and S. Grimme            |
        |                                                              |
        |______________________________________________________________|

        Please cite:
        (TBA)
        S. Grimme, F. Bohle, A. Hansen, P. Pracht, S. Spicher, and M. Stahn
        J. Phys. Chem. A 2021, 125, 19, 4039-4054.
        DOI: https://doi.org/10.1021/acs.jpca.1c00971

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

START_DESCR = "Energetic sorting of Conformer Rotamer Ensembles (command line version)."

DIGILEN = 60

PLENGTH = 100

AU2J = 4.3597482e-18  # a.u.(hartree/mol) to J
KB = 1.3806485279e-23  # J/K
R = 1.987203585e-03  # kcal/(mol*K)
AU2KCAL = 627.50947428
BOHR2ANG = 0.52917721067
PLANCK = 6.62607015e-34
C = 2.998e8
WARNLEN = max(len(i) for i in ["WARNING:", "ERROR:", "INFORMATION:"]) + 1

si_bib = {
    "tm": [
        r"@misc{TURBOMOLE,",
        r"  title = {{TURBOMOLE V7.5 2020}, a development of {University of Karlsruhe} and",
        r"  {Forschungszentrum Karlsruhe GmbH}, 1989-2007,",
        r"  {TURBOMOLE GmbH}, since 2007; available from \\",
        r"  {\tt https://www.turbomole.org}.}",
        r"}",
        r"@Article{TURBOMOLE.2020,",
        r"  author = {Balasubramani, Sree Ganesh  and Chen, Guo P. ",
        r"      and Coriani, Sonia and Diedenhofen, Michael and ",
        r"      Frank, Marius S. and Franzke, Yannick J. and ",
        r"      Furche, Filipp and Grotjahn, Robin and Harding, Michael E. ",
        r"      and H\"attig, Christof and Hellweg, Arnim and ",
        r"      Helmich-Paris, Benjamin and Holzer, Christof and Huniar, Uwe",
        r"      and Kaupp, Martin and Marefat Khah, Alireza ",
        r"      and Karbalaei Khani, Sarah and M\"uller, Thomas and Mack, Fabian",
        r"      and Nguyen, Brian D. and Parker, Shane M. and Perlt, Eva ",
        r"      and Rappoport, Dmitrij and Reiter, Kevin and Roy, Saswata and",
        r"      R\"uckert, Matthias and Schmitz, Gunnar and Sierka, Marek",
        r"      and Tapavicza, Enrico and Tew, David P. and van W\"ullen, Christoph",
        r"      and Voora, Vamsee K. and Weigend, Florian and",
        r"      Wody{\’n}ski, Artur and Yu, Jason M.},",
        r"  title = {TURBOMOLE: Modular program suite for \textit{ab initio}",
        r"            quantum-chemical and condensed-matter simulations},",
        r"  journal   = {J. Chem. Phys.},",
        r"  volume    = {152},",
        r"  issue     = {18},",
        r"  pages     = {184107},",
        r"  year      = {2020},",
        r"  url       = { https://doi.org/10.1063/5.0004635},",
        r"  DOI       = {10.1063/5.0004635}",
        r"}",
    ],
    "orca": [
        r"@article{ORCA_generic,",
        r"  author = {Neese, Frank},",
        r"  title = {The ORCA program system},",
        r"  journal = {WIREs Computational Molecular Science},",
        r"  volume = {2},",
        r"  number = {1},",
        r"  pages = {73-78},",
        r"  doi = {https://doi.org/10.1002/wcms.81},",
        r"  url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.81},",
        r"  eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/wcms.81},",
        r"  year = {2012}",
        r"}",
        r"@article{ORCA4.0,",
        r"  author = {Neese, Frank},",
        r"  title = {Software update: the ORCA program system, version 4.0},",
        r"  journal = {WIREs Computational Molecular Science},",
        r"  volume = {8},",
        r"  number = {1},",
        r"  pages = {e1327},",
        r"  doi = {https://doi.org/10.1002/wcms.1327},",
        r"  url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1327},",
        r"  eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/wcms.1327},",
        r"  year = {2018}",
        r"}",
    ],
    "cosmors": [
        r"@article{Klamt1995,",
        r"  author = {Klamt, Andreas},",
        r"  doi = {10.1021/j100007a062},",
        r"  journal = {The Journal of Physical Chemistry},",
        r"  number = {7},",
        r"  pages = {2224--2235},",
        r"  title = {{Conductor-like Screening Model for Real Solvents: A New "
        + "Approach to the Quantitative Calculation of Solvation Phenomena}},",
        r"  url = {https://pubs.acs.org/doi/abs/10.1021/j100007a062},",
        r"  volume = {99},",
        r"  year = {1995}",
        r"}",
        r"@article{Klamt1998,",
        r"  author = {Klamt, Andreas and Jonas, Volker and B{\"{u}}rger, Thorsten"
        + " and Lohrenz, John C. W.},",
        r"  doi = {10.1021/jp980017s},",
        r"  journal = {The Journal of Physical Chemistry A},",
        r"  number = {26},",
        r"  pages = {5074--5085},",
        r"  title = {{Refinement and Parametrization of COSMO-RS}},",
        r"  url = {https://pubs.acs.org/doi/10.1021/jp980017s},",
        r"  volume = {102},",
        r"  year = {1998}",
        r"}",
        r"@article{Eckert2002,",
        r"  author = {Eckert, Frank and Klamt, Andreas},",
        r"  doi = {10.1002/aic.690480220},",
        r"  journal = {AIChE Journal},",
        r"  number = {2},",
        r"  pages = {369--385},",
        r"  title = {{Fast solvent screening via quantum chemistry: COSMO-RS approach}},",
        r"  url = {http://doi.wiley.com/10.1002/aic.690480220},",
        r"  volume = {48},",
        r"  year = {2002}",
        r"}",
        r"@misc{COSMOtherm,",
        r"   title = {COSMOtherm, Release 19; COSMOlogic GmbH & Co. KG, "
        + "{\tt http://www.cosmologic.de}.}",
        r"}",
    ],
    "xtb": [
        r"@article{xtb_generic,",
        r"  title={Extended tight‐binding quantum chemistry methods},",
        r"  author={Bannwarth, Christoph and Caldeweyher, Eike and Ehlert,"
        + " Sebastian and Hansen, Andreas and Pracht, Philipp and Seibert, Jakob "
        + "and Spicher, Spicher and Grimme, Stefan},",
        r"  journal={{WIREs} Comput{.} Mol{.} Sci{.}},",
        r"  volume = {11},",
        r"  year={2020},",
        r"  pages={e01493},",
        r"  doi={10.1002/wcms.1493},",
        r"  url={https://dx.doi.org/10.1002/wcms.1493}",
        r"}",
        r"@article{GFN2,",
        r"  title={GFN2-xTB—An accurate and broadly parametrized self-consistent"
        + " tight-binding quantum chemical method with multipole electrostatics and"
        + " density-dependent dispersion contributions},",
        r"  author={Bannwarth, Christoph and Ehlert, Sebastian and Grimme, Stefan},",
        r"  journal={J{.} Chem{.} Theory Comput{.}},",
        r"  volume={15},",
        r"  number={3},",
        r"  pages={1652--1671},",
        r"  year={2019},",
        r"  doi={10.1021/acs.jctc.8b01176},",
        r"  url={https://dx.doi.org/10.1021/acs.jctc.8b01176},",
        r"}",
        r"@article{GFN1,",
        r"  title={A robust and accurate tight-binding quantum chemical method "
        + "for structures, vibrational frequencies, and noncovalent interactions of"
        + " large molecular systems parametrized for all spd-block elements (Z=1--86)},",
        r"  author={Grimme, Stefan and Bannwarth, Christoph and Shushkov, Philip},",
        r"  journal={J{.} Chem{.} Theory Comput{.}},",
        r"  volume={13},",
        r"  number={5},",
        r"  pages={1989--2009},",
        r"  year={2017},",
        r"  doi={10.1021/acs.jctc.7b00118},",
        r"  url={https://dx.doi.org/10.1021/acs.jctc.7b00118},",
        r"}",
    ],
    "censo": [
        r"@article{CENSO",
        r"  author = {Grimme, Stefan and Bohle, Fabian and Hansen, Andreas and "
        + "Pracht, Philipp and Spicher, Sebastian and Stahn, Marcel},",
        r"  title = {Efficient Quantum Chemical Calculation of Structure Ensembles"
        + " and Free Energies for Nonrigid Molecules},",
        r"  journal = {The Journal of Physical Chemistry A},",
        r"  volume = {125},",
        r"  number = {19},",
        r"  pages = {4039--4054},",
        r"  year = {2021},",
        r"  doi = {10.1021/acs.jpca.1c00971},",
        r"  note ={PMID: 33688730},",
        r"  URL = {https://doi.org/10.1021/acs.jpca.1c00971},",
        r"  eprint = {https://doi.org/10.1021/acs.jpca.1c00971}",
        r"}",
    ],
    # selected functionals:
    "r2scan-3c": [
        r"@article{r2scan-3c,",
        r"  author = {Grimme, Stefan and Hansen, Andreas and Ehlert, Sebastian "
        + "and Mewes, Jan-Michael},",
        r"  doi = {10.1063/5.0040021},",
        r"  journal = {The Journal of Chemical Physics},",
        r"  number = {6},",
        r"  pages = {064103},",
        r"  title = {{r 2 SCAN-3c: A “Swiss army knife” composite electronic-structure method}},",
        r"  url = {https://doi.org/10.1063/5.0040021},",
        r"  volume = {154},",
        r"  year = {2021}",
        r"}",
    ],
    "pbeh-3c": [
        r"@article{Grimme2015,",
        r"  author = {Grimme, Stefan and Brandenburg, Jan Gerit and Bannwarth, "
        + "Christoph and Hansen, Andreas},",
        r"  doi = {10.1063/1.4927476},",
        r"  journal = {The Journal of Chemical Physics},",
        r"  number = {5},",
        r"  pages = {054107},",
        r"  pmid = {26254642},",
        r"  title = {{Consistent structures and interactions by density functional "
        + "theory with small atomic orbital basis sets}},",
        r"  url = {http://scitation.aip.org/content/aip/journal/jcp/143/5/10.1063/1.4927476},",
        r"  volume = {143},",
        r"  year = {2015}",
        r"}",
    ],
    "b97-3c": [
        r"@article{b97-3c,",
        r"  author = {Brandenburg, Jan Gerit and Bannwarth, Christoph and Hansen, "
        + "Andreas and Grimme, Stefan},",
        r"  doi = {10.1063/1.5012601},",
        r"  journal = {The Journal of Chemical Physics},",
        r"  number = {6},",
        r"  pages = {064104},",
        r"  title = {{B97-3c: A revised low-cost variant of the B97-D density functional method}},",
        r"  url = {http://aip.scitation.org/doi/10.1063/1.5012601},",
        r"  volume = {148},",
        r"  year = {2018}",
        r"}",
    ],
    "pbe0": [
        r"@article{pbe0_one,",
        r"  author = {Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias},",
        r"  doi = {10.1103/PhysRevLett.77.3865},",
        r"  journal = {Physical Review Letters},",
        r"  number = {18},",
        r"  pages = {3865--3868},",
        r"  title = {{Generalized Gradient Approximation Made Simple}},",
        r"  url = {https://link.aps.org/doi/10.1103/PhysRevLett.77.3865},",
        r"  volume = {77},",
        r"  year = {1996}",
        r"}",
        r"@article{pbe0_one_erratum,",
        r"  title = {Generalized Gradient Approximation Made Simple "
        + "[Phys. Rev. Lett. 77, 3865 (1996)]},",
        r"  author = {Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias},",
        r"  journal = {Phys. Rev. Lett.},",
        r"  volume = {78},",
        r"  issue = {7},",
        r"  pages = {1396--1396},",
        r"  numpages = {0},",
        r"  year = {1997},",
        r"  month = {Feb},",
        r"  publisher = {American Physical Society},",
        r"  doi = {10.1103/PhysRevLett.78.1396},",
        r"  url = {https://link.aps.org/doi/10.1103/PhysRevLett.78.1396}",
        r"}",
        r"@article{PBE0_two,",
        r"  author = {Adamo,Carlo and Barone,Vincenzo},",
        r"  title = {Toward reliable density functional methods without "
        + "adjustable parameters: The PBE0 model},",
        r"  journal = {The Journal of Chemical Physics},",
        r"  volume = {110},",
        r"  number = {13},",
        r"  pages = {6158-6170},",
        r"  year = {1999},",
        r"  doi = {10.1063/1.478522},",
        r"  url = {https://doi.org/10.1063/1.478522},",
        r"  eprint = {https://doi.org/10.1063/1.478522}",
        r"}",
    ],
    "b3-lyp": [
        r"@article{b3lyp_one,",
        r"  author = {Becke,Axel D. },",
        r"  title = {A new mixing of Hartree–Fock and local density‐functional theories},",
        r"  journal = {The Journal of Chemical Physics},",
        r"  volume = {98},",
        r"  number = {2},",
        r"  pages = {1372-1377},",
        r"  year = {1993},",
        r"  doi = {10.1063/1.464304},",
        r"  url = {https://doi.org/10.1063/1.464304},",
        r"  eprint = {https://doi.org/10.1063/1.464304}",
        r"}",
        r"@article{b3lyp_two,",
        r"  author = {Stephens, P. J. and Devlin, F. J. and Chabalowski, C. F."
        + " and Frisch, M. J.},",
        r"  title = {Ab Initio Calculation of Vibrational Absorption and Circular"
        + " Dichroism Spectra Using Density Functional Force Fields},",
        r"  journal = {The Journal of Physical Chemistry},",
        r"  volume = {98},",
        r"  number = {45},",
        r"  pages = {11623-11627},",
        r"  year = {1994},",
        r"  doi = {10.1021/j100096a001},",
        r"  url = {https://doi.org/10.1021/j100096a001},",
        r"  eprint = {https://doi.org/10.1021/j100096a001}",
        r"}",
    ],
    "def2_basis": [
        r"@article{Weigend2005,",
        r"  author = {Weigend, Florian and Ahlrichs, Reinhart},",
        r"  doi = {10.1039/b508541a},",
        r"  journal = {Physical Chemistry Chemical Physics},",
        r"  number = {18},",
        r"  pages = {3297},",
        r"  title = {{Balanced basis sets of split valence, triple zeta valence "
        + "and quadruple zeta valence quality for H to Rn: Design and assessment of accuracy}},",
        r"  url = {http://xlink.rsc.org/?DOI=b508541a},",
        r"  volume = {7},",
        r"  year = {2005}",
        r"}",
    ],
    "def2_auxbasis": [
        r"@article{Eichkorn1997,",
        r"author = {Eichkorn, Karin and Weigend, Florian and Treutler, Oliver "
        + "and Ahlrichs, Reinhart},",
        r"doi = {10.1007/s002140050244},",
        r"journal = {Theoretical Chemistry Accounts: Theory, Computation, and "
        + "Modeling (Theoretica Chimica Acta)},",
        r"number = {1-4},",
        r"pages = {119--124},",
        r"publisher = {Springer-Verlag},",
        r"title = {{Auxiliary basis sets for main row atoms and transition metals"
        + " and their use to approximate Coulomb potentials}},",
        r"url = {http://link.springer.com/10.1007/s002140050244},",
        r"volume = {97},",
        r"year = {1997}",
        r"}",
        r"@article{Weigend2006,",
        r"  author = {Weigend, Florian},",
        r"  doi = {10.1039/b515623h},",
        r"  journal = {Physical Chemistry Chemical Physics},",
        r"  number = {9},",
        r"  pages = {1057},",
        r"  title = {{Accurate Coulomb-fitting basis sets for H to Rn}},",
        r"  url = {http://xlink.rsc.org/?DOI=b515623h},",
        r"  volume = {8},",
        r"  year = {2006}",
        r"}",
    ],
    "ri": [
        r"@article{RI,",
        r"title = {Integral approximations for LCAO-SCF calculations},",
        r"journal = {Chemical Physics Letters},",
        r"volume = {213},",
        r"number = {5},",
        r"pages = {514-518},",
        r"year = {1993},",
        r"issn = {0009-2614},",
        r"doi = {https://doi.org/10.1016/0009-2614(93)89151-7},",
        r"url = {https://www.sciencedirect.com/science/article/pii/0009261493891517},",
        r"author = {O. Vahtras and J. Almlöf and M.W. Feyereisen},",
        r"}",
    ],
    "sph": [
        r"@article{SPH,",
        r"  author = {Spicher, Sebastian and Grimme, Stefan},",
        r"  doi = {10.1021/acs.jctc.0c01306},",
        r"  journal = {Journal of Chemical Theory and Computation},",
        r"  number = {3},",
        r"  pages = {1701--1714},",
        r"  pmid = {33554604},",
        r"  title = {{Single-Point Hessian Calculations for Improved Vibrational "
        + "Frequencies and Rigid-Rotor-Harmonic-Oscillator Thermodynamics}},",
        r"  url = {https://pubs.acs.org/doi/10.1021/acs.jctc.0c01306},",
        r"  volume = {17},",
        r"  year = {2021}",
        r"}",
    ],
    "dcosmors": [
        r"@article{DCOSMORS,",
        r"  author = {Klamt, Andreas and Diedenhofen, Michael},",
        r"  doi = {10.1021/jp511158y},",
        r"  journal = {The Journal of Physical Chemistry A},",
        r"  number = {21},",
        r"  pages = {5439--5445},",
        r"  title = {{Calculation of Solvation Free Energies with DCOSMO-RS}},",
        r"  url = {http://pubs.acs.org/doi/10.1021/jp511158y},",
        r"  volume = {119},",
        r"  year = {2015}",
        r"}",
    ],
    "alpb": [
        r"@article{ALPB,",
        r"  author = {Ehlert, Sebastian and Stahn, Marcel and Spicher, Sebastian and Grimme, Stefan},",
        r"  title = {Robust and Efficient Implicit Solvation Model for Fast Semiempirical Methods},",
        r"  journal = {Journal of Chemical Theory and Computation},",
        r"  volume = {0},",
        r"  number = {0},",
        r"  pages = {null},",
        r"  year = {0},",
        r"  doi = {10.1021/acs.jctc.1c00471},",
        r"  note ={PMID: 34185531},",
        r"  URL = {https://doi.org/10.1021/acs.jctc.1c00471},",
        r"}",
    ],
}
