"""
Storing constants for the use in all CENSO modules.
Storing program paths --> still in transition
Storing censo_solvent_db solvent database across all solvation models (as fallback)
"""
import os

__version__ = "1.0.0"
DESCR = f"""
         ______________________________________________________________
        |                                                              |
        |                                                              |
        |                   CENSO - Commandline ENSO                   |
        |                           v {__version__:<{19}}              |
        |    energetic sorting of CREST Conformer Rotamer Ensembles    |
        |                    University of Bonn, MCTC                  |
        |                          June 2020                           |
        |                 based on ENSO version 2.0.1                  |
        |                     F. Bohle and S. Grimme                   |
        |                                                              |
        |______________________________________________________________|
        
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""
global ENVIRON
ENVIRON = os.environ.copy()
CODING = "ISO-8859-1"
DIGILEN = 60
PLENGTH = 100
AU2J = 4.3597482e-18  # a.u.(hartree/mol) to J
KB = 1.3806485279e-23  # J/K
AU2KCAL = 627.50947428
BOHR2ANG = 0.52917721067

# program paths:

external_paths = {}
external_paths["orcapath"] = ""
external_paths["orcaversion"] = ""
external_paths["xtbpath"] = ""
external_paths["crestpath"] = ""
external_paths["cosmorssetup"] = ""
external_paths["dbpath"] = ""
external_paths["cosmothermversion"] = ""
external_paths["mpshiftpath"] = ""
external_paths["escfpath"] = ""
external_paths["cefinepath"] = ""


# censo solvent database to chose solvents across all available solvent models
censo_solvent_db = {
    "acetone": {
        "cosmors": ["propanone_c0", "propanone_c0"],
        "dcosmors": ["propanone", "propanone"],
        "xtb": ["acetone", "acetone"],
        "cpcm": ["acetone", "acetone"],
        "smd": ["ACETONE", "ACETONE"],
        "DC": 20.7,
    },
    "chcl3": {
        "cosmors": ["chcl3_c0", "chcl3_c0"],
        "dcosmors": ["chcl3", "chcl3"],
        "xtb": ["chcl3", "chcl3"],
        "cpcm": ["chloroform", "chloroform"],
        "smd": ["CHLOROFORM", "CHLOROFORM"],
        "DC": 4.8,
    },
    "acetonitrile": {
        "cosmors": ["acetonitrile_c0", "acetonitrile_c0"],
        "dcosmors": ["acetonitrile", "acetonitrile"],
        "xtb": ["acetonitrile", "acetonitrile"],
        "cpcm": ["acetonitrile", "acetonitrile"],
        "smd": ["ACETONITRILE", "ACETONITRILE"],
        "DC": 36.6,
    },
    "ch2cl2": {
        "cosmors": ["ch2cl2_c0", "ch2cl2_c0"],
        "dcosmors": [None, "chcl3"],
        "xtb": ["ch2cl2", "ch2cl2"],
        "cpcm": ["CH2Cl2", "CH2Cl2"],
        "smd": ["DICHLOROMETHANE", "DICHLOROMETHANE"],
        "DC": 9.1,
    },
    "dmso": {
        "cosmors": ["dimethylsulfoxide_c0", "dimethylsulfoxide_c0"],
        "dcosmors": ["dimethylsulfoxide", "dimethylsulfoxide"],
        "xtb": ["dmso", "dmso"],
        "cpcm": ["DMSO", "DMSO"],
        "smd": ["DIMETHYLSULFOXIDE", "DIMETHYLSULFOXIDE"],
        "DC": 47.2,
    },
    "h2o": {
        "cosmors": ["h2o_c0", "h2o_c0"],
        "dcosmors": ["h2o", "h2o"],
        "xtb": ["h2o", "h2o"],
        "cpcm": ["Water", "Water"],
        "smd": ["WATER", "WATER"],
        "DC": 80.1,
    },
    "methanol": {
        "cosmors": ["methanol_c0", "methanol_c0"],
        "dcosmors": ["methanol", "methanol"],
        "xtb": ["methanol", "methanol"],
        "cpcm": ["Methanol", "Methanol"],
        "smd": ["METHANOL", "METHANOL"],
        "DC": 32.7,
    },
    "thf": {
        "cosmors": ["thf_c0", "thf_c0"],
        "dcosmors": ["thf", "thf"],
        "xtb": ["thf", "thf"],
        "cpcm": ["THF", "THF"],
        "smd": ["TETRAHYDROFURAN", "TETRAHYDROFURAN"],
        "DC": 7.6,
    },
    "toluene": {
        "cosmors": ["toluene_c0", "toluene_c0"],
        "dcosmors": ["toluene", "toluene"],
        "xtb": ["toluene", "toluene"],
        "cpcm": ["Toluene", "Toluene"],
        "smd": ["TOLUENE", "TOLUENE"],
        "DC": 2.4,
    },
    "octanol": {
        "cosmors": ["1-octanol_c0", "1-octanol_c0"],
        "dcosmors": ["octanol", "octanol"],
        "xtb": ["octanol", "octanol"],
        "cpcm": ["Octanol", "Octanol"],
        "smd": ["1-OCTANOL", "1-OCTANOL"],
        "DC": 9.9,
    },
    "woctanol": {
        "cosmors": [None, "woctanol"],
        "dcosmors": ["wet-otcanol", "wet-octanol"],
        "xtb": ["woctanol", "woctanol"],
        "cpcm": [None, "Octanol"],
        "smd": [None, "1-OCTANOL"],
        "DC": 8.1,
    },
    "hexadecane": {
        "cosmors": ["n-hexadecane_c0", "n-hexadecane_c0"],
        "dcosmors": ["hexadecane", "hexadecane"],
        "xtb": ["hexadecane", "hexadecane"],
        "cpcm": [None, "Hexane"],
        "smd": ["N-HEXADECANE", "N-HEXADECANE"],
        "DC": 2.1,
    },
    "dmf": {
        "cosmors": ["dimethylformamide_c0", "dimethylformamide_c0"],
        "dcosmors": [None, "dimethylsulfoxide"],
        "xtb": ["dmf", "dmf"],
        "cpcm": ["DMF", "DMF"],
        "smd": ["N,N-DIMETHYLFORMAMIDE", "N,N-DIMETHYLFORMAMIDE"],
        "DC": 38.3,
    },
    "aniline": {
        "cosmors": ["aniline_c0", "aniline_c0"],
        "dcosmors": ["aniline", "aniline"],
        "xtb": ["aniline", "aniline"],
        "cpcm": [None, "Pyridine"],
        "smd": ["ANILINE", "ANILINE"],
        "DC": 6.9,
    },
    "cyclohexane": {
        "cosmors": ["cyclohexane_c0", "cyclohexane_c0"],
        "dcosmors": ["cyclohexane", "cyclohexane"],
        "xtb": [None, "hexane"],
        "cpcm": ["Cyclohexane", "Cyclohexane"],
        "smd": ["CYCLOHEXANE", "CYCLOHEXANE"],
        "DC": 2.0,
    },
    "ccl4": {
        "cosmors": ["ccl4_c0", "ccl4_c0"],
        "dcosmors": ["ccl4", "ccl4"],
        "xtb": ["ccl4", "ccl4"],
        "cpcm": ["CCl4", "CCl4"],
        "smd": ["CARBON TETRACHLORIDE", "CARBON TETRACHLORIDE"],
        "DC": 2.2,
    },
    "diethylether": {
        "cosmors": ["diethylether_c0", "diethylether_c0"],
        "dcosmors": ["diethylether", "diethylether"],
        "xtb": ["ether", "ether"],
        "cpcm": [None, "THF"],
        "smd": ["DIETHYL ETHER", "DIETHYL ETHER"],
        "DC": 4.4,
    },
    "ethanol": {
        "cosmors": ["ethanol_c0", "ethanol_c0"],
        "dcosmors": ["ethanol", "ethanol"],
        "xtb": [None, "methanol"],
        "cpcm": [None, "Methanol"],
        "smd": ["ETHANOL", "ETHANOL"],
        "DC": 24.6,
    },
    "hexane": {
        "cosmors": ["hexane_c0", "hexane_c0"],
        "dcosmors": ["hexane", "hexane"],
        "xtb": ["hexane", "hexane"],
        "cpcm": ["Hexane", "Hexane"],
        "smd": ["N-HEXANE", "N-HEXANE"],
        "DC": 1.9,
    },
    "nitromethane": {
        "cosmors": ["nitromethane_c0", "nitromethane_c0"],
        "dcosmors": ["nitromethane", "nitromethane"],
        "xtb": ["nitromethane", "nitromethane"],
        "cpcm": [None, "methanol"],
        "smd": "",
        "DC": 38.2,
    },
    "benzaldehyde": {
        "cosmors": ["benzaldehyde_c0", "benzaldehyde_c0"],
        "dcosmors": [None, "propanone"],
        "xtb": ["benzaldehyde", "benzaldehyde"],
        "cpcm": [None, "Pyridine"],
        "smd": ["BENZALDEHYDE", "BENZALDEHYDE"],
        "DC": 18.2,
    },
    "benzene": {
        "cosmors": ["benzene_c0", "benzene_c0"],
        "dcosmors": [None, "toluene"],
        "xtb": ["benzene", "benzene"],
        "cpcm": ["Benzene", "Benzene"],
        "smd": ["BENZENE", "BENZENE"],
        "DC": 2.3,
    },
    "cs2": {
        "cosmors": ["cs2_c0", "cs2_c0"],
        "dcosmors": [None, "ccl4"],
        "xtb": ["cs2", "cs2"],
        "cpcm": [None, "CCl4"],
        "smd": ["CARBON DISULFIDE", "CARBON DISULFIDE"],
        "DC": 2.6,
    },
    "dioxane": {
        "cosmors": ["dioxane_c0", "dioxane_c0"],
        "dcosmors": [None, "diethylether"],
        "xtb": ["dioxane", "dioxane"],
        "cpcm": [None, "Cyclohexane"],
        "smd": ["1,4-DIOXANE", "1,4-DIOXANE"],
        "DC": 2.2,
    },
    "ethylacetate": {
        "cosmors": ["ethylacetate_c0", "ethylacetate_c0"],
        "dcosmors": [None, "diethylether"],
        "xtb": ["ethylacetate", "ethylacetate"],
        "cpcm": [None, "THF"],
        "smd": ["ETHYL ETHANOATE", "ETHYL ETHANOATE"],
        "DC": 5.9,
    },
    "furan": {
        "cosmors": ["furane_c0", "furane_c0"],
        "dcosmors": [None, "diethylether"],
        "xtb": ["furane", "furane"],
        "cpcm": [None, "THF"],
        "smd": [None, "THF"],
        "DC": 3.0,
    },
    "phenol": {
        "cosmors": ["phenol_c0", "phenol_c0"],
        "dcosmors": [None, "thf"],
        "xtb": ["phenol", "phenol"],
        "cpcm": [None, "THF"],
        "smd": [None, "THIOPHENOL"],
        "DC": 8.0,
    },
}
