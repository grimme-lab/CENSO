"""
Storing constants for the use in all CENSO modules.
Storing program paths --> still in transition
Storing censo_solvent_db solvent database across all solvation models (as fallback)
"""
import os
import sys

__version__ = "1.1.1"

DESCR = f"""
         ______________________________________________________________
        |                                                              |
        |                                                              |
        |                   CENSO - Commandline ENSO                   |
        |                           v {__version__:<{19}}              |
        |    energetic sorting of CREST Conformer Rotamer Ensembles    |
        |                    University of Bonn, MCTC                  |
        |                           Feb 2021                           |
        |                 based on ENSO version 2.0.1                  |
        |                     F. Bohle and S. Grimme                   |
        |                                                              |
        |______________________________________________________________|

        Please cite: 
        S. Grimme, F. Bohle, A. Hansen, P. Pracht, S. Spicher, and M. Stahn 
        J. Phys. Chem. A 2021, 125, 19, 4039-4054.
        DOI: https://doi.org/10.1021/acs.jpca.1c00971
        
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

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
DIGILEN = 60
PLENGTH = 100
AU2J = 4.3597482e-18  # a.u.(hartree/mol) to J
KB = 1.3806485279e-23  # J/K
R = 1.987203585e-03  # kcal/(mol*K)
AU2KCAL = 627.50947428
BOHR2ANG = 0.52917721067
WARNLEN = max([len(i) for i in ["WARNING:", "ERROR:", "INFORMATION:"]]) + 1


class dfa_settings:
    """Contains information on available functionals, per part, qm code, dispersion
    correction..."""

    # definitions:
    # gga_dfa = ("tpss", "pbe", "kt1", "kt2", "b97-d3", "b97-d")
    # hybrid_dfa = (
    #     "pbe0",
    #     "pw6b95",
    #     "wb97x-d3",
    #     "wb97x-d3bj",
    #     "wb97x-v",
    #     "cam-b3lyp",
    #     "b3-lyp",  # tm
    #     "b3lyp",  # orca
    #     "pbeh-3c",
    #     "m06x",
    #     "bh-lyp",
    #     "tpssh",
    # )
    # dh_dfa = ("dsd-blyp",)
    # disp_already_included_in_func = composite_dfa + (
    #     "b97-d3",  # orca
    #     "b97-d",  # tm
    #     "wb97x-d3",
    #     "wb97x-d3bj",
    #     "wb97x-v",
    # )
    # func_orca = ["pbeh-3c", "b97-3c", "tpss", "b97-d3", "pbe"]
    # func_tm = ["pbeh-3c", "b97-3c", "tpss", "r2scan-3c", "b97-d", "pbe"]
    # func3_orca = [
    #     "pw6b95",
    #     "pbe0",
    #     "b97-d3",
    #     "wb97x-d3",
    #     "wb97x-d3bj",
    #     "wb97x-v",
    #     "dsd-blyp",
    #     "b97-3c",
    #     "pbeh-3c",
    #     "tpss",
    #     "pbe",
    #     "b3lyp",
    # ]
    # func3_tm = [
    #     "pw6b95",
    #     "pbe0",
    #     "b97-d",
    #     "r2scan-3c",
    #     "b97-3c",
    #     "wb97x-v",
    #     "pbeh-3c",
    #     "tpss",
    #     "pbe",
    #     "b3-lyp",
    # ]
    # func_j_tm = ["tpss", "pbe0", "pbeh-3c", "r2scan-3c", "pbe"]
    # func_j_orca = ["tpss", "pbe0", "pbeh-3c", "pbe"]
    # func_s_tm = ["tpss", "pbe0", "pbeh-3c", "b97-3c", "kt1", "kt2", "r2scan-3c", "pbe"]
    # func_s_orca = ["tpss", "pbe0", "dsd-blyp", "pbeh-3c", "kt2", "pbe"]
    # "func0"
    # "func"
    # "func3"
    # "func_j"
    # "func_s"
    # "func_or"
    # "func_or_scf"

    # d3
    # d3(0)
    # d4
    # nl
    # novdw
    # included
    # composite

    # type: {'global_hybrid',
    #        'mGGA',
    #        'composite_mGGA',
    #        'doublehybrid',
    #        'rsh_hybrid',
    #        'composite_hybrid',
    #        'GGA',
    #        'composite_gga'
    # }
    composite_method_basis = {
        "pbeh-3c": "def2-mSVP",
        "pbeh3c": "def2-mSVP",
        "b97-3c": "def2-mTZVP",
        "b973c": "def2-mTZVP",
        "hf3c": "minix",
        "hf-3c": "minix",
        "r2scan-3c": "def2-mTZVPP",
        "r2scan3c": "def2-mTZVPP",
        "b3lyp-3c": "def2-mSVP",
        "b3lyp3c": "def2-mSVP",
    }

    relay_functionals = {
        "pbeh-3c": "pbeh-3c",
        "b97-3c": "b97-3c",
        "r2scan-3c": "r2scan-3c",
        "b3lyp-3c": "b3lyp-3c",
        "pbe": "pbe-d4",
        "tpss": "tpss-d4",
        "b97-d": "b97-d3(0)",
        "kt1": "kt1-novdw",
        "kt2": "kt2-novdw",
        "pbe0": "pbe0-d4",
        "pw6b95": "pw6b95-d4",
        "b3lyp": "b3lyp-d4",
        "b3-lyp": "b3lyp-d4",
        "wb97x-v": "wb97x-v",
        "wb97x-d3": "wb97x-d3",
        "wb97x-d3bj": "wb97x-d3bj",
        "dsd-blyp": "dsd-blyp-d3",
    }
    functionals = {
        "pbeh-3c": {
            "tm": "pbeh-3c",
            "orca": "pbeh-3c",
            "disp": "composite",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "composite_hybrid",
        },
        "b97-3c": {
            "tm": "b97-3c",
            "orca": "b97-3c",
            "disp": "composite",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "composite_gga",
        },
        "r2scan-3c": {
            "tm": "r2scan-3c",
            "orca": "r2scan-3c",
            "disp": "composite",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "composite_mGGA",
        },
        "b3lyp-3c": {
            "tm": "b3lyp-3c",
            "orca": "b3lyp",
            "disp": "composite",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "composite_hybrid",
        },
        "pbe-novdw": {
            "tm": "pbe",
            "orca": "pbe",
            "disp": "novdw",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "GGA",
        },
        "pbe-d3": {
            "tm": "pbe",
            "orca": "pbe",
            "disp": "d3bj",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "GGA",
        },
        "pbe-d3(0)": {
            "tm": "pbe",
            "orca": "pbe",
            "disp": "d3(0)",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "GGA",
        },
        "pbe-d4": {
            "tm": "pbe",
            "orca": "pbe",
            "disp": "d4",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "GGA",
        },
        "pbe-nl": {
            "tm": "pbe",
            "orca": None,
            "disp": "nl",
            "part": ["func", "func3", "func_j", "func_s", "func_or", "func_or_scf"],
            "type": "GGA",
        },
        "tpss-novdw": {
            "tm": "tpss",
            "orca": "tpss",
            "disp": "novdw",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "mGGA",
        },
        "tpss-d3": {
            "tm": "tpss",
            "orca": "tpss",
            "disp": "d3bj",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "mGGA",
        },
        "tpss-d3(0)": {
            "tm": "tpss",
            "orca": "tpss",
            "disp": "d3(0)",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "mGGA",
        },
        "tpss-d4": {
            "tm": "tpss",
            "orca": "tpss",
            "disp": "d4",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "mGGA",
        },
        "tpss-nl": {
            "tm": "tpss",
            "orca": None,
            "disp": "nl",
            "part": ["func", "func3", "func_j", "func_s", "func_or", "func_or_scf"],
            "type": "mGGA",
        },
        "b97-d3": {
            "tm": "b97-d",
            "orca": "b97-d3",
            "disp": "d3bj",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "GGA",
        },
        "b97-d3(0)": {
            "tm": "b97-d",
            "orca": None,
            "disp": "d3(0)",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "GGA",
        },
        "kt1-novdw": {
            "tm": "kt1",
            "orca": None,
            "disp": "novdw",
            "part": ["func_j", "func_s"],
            "type": "GGA",
        },
        "kt2-novdw": {
            "tm": "kt2",
            "orca": "kt2",
            "disp": "novdw",
            "part": ["func_j", "func_s"],
            "type": "GGA",
        },
        "pbe0-novdw": {
            "tm": "pbe0",
            "orca": "pbe0",
            "disp": "novdw",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pbe0-d3": {
            "tm": "pbe0",
            "orca": "pbe0",
            "disp": "d3bj",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pbe0-d3(0)": {
            "tm": "pbe0",
            "orca": "pbe0",
            "disp": "d3(0)",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pbe0-d4": {
            "tm": "pbe0",
            "orca": "pbe0",
            "disp": "d4",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pbe0-nl": {
            "tm": "pbe0",
            "orca": None,
            "disp": "nl",
            "part": ["func", "func3", "func_j", "func_s", "func_or", "func_or_scf"],
            "type": "global_hybrid",
        },
        "pw6b95-novdw": {
            "tm": "pw6b95",
            "orca": "pw6b95",
            "disp": "novdw",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pw6b95-d3": {
            "tm": "pw6b95",
            "orca": "pw6b95",
            "disp": "d3bj",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pw6b95-d3(0)": {
            "tm": "pw6b95",
            "orca": "pw6b95",
            "disp": "d3(0)",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "pw6b95-d4": {
            "tm": "pw6b95",
            "orca": "pw6b95",
            "disp": "d4",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        # "pw6b95-nl":{
        #     "tm": "pw6b95",
        #     "orca": None,
        #     "disp": "nl",
        #     "part": ["func", "func3", "func_j", "func_s", "func_or", "func_or_scf"],
        #     "type": "global_hybrid",
        # },
        "b3lyp-novdw": {
            "tm": "b3-lyp",
            "orca": "b3lyp",
            "disp": "novdw",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "b3lyp-d3": {
            "tm": "b3-lyp",
            "orca": "b3lyp",
            "disp": "d3bj",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "b3lyp-d3(0)": {
            "tm": "b3-lyp",
            "orca": "b3lyp",
            "disp": "d3(0)",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "b3lyp-d4": {
            "tm": "b3-lyp",
            "orca": "b3lyp",
            "disp": "d4",
            "part": [
                "func0",
                "func",
                "func3",
                "func_j",
                "func_s",
                "func_or",
                "func_or_scf",
            ],
            "type": "global_hybrid",
        },
        "b3lyp-nl": {
            "tm": "b3-lyp",
            "orca": "b3lyp",
            "disp": "nl",
            "part": ["func", "func3", "func_j", "func_s", "func_or", "func_or_scf"],
            "type": "global_hybrid",
        },
        "wb97x-v": {
            "tm": "wb97x-v",
            "orca": "wb97x-v",
            "disp": "included",
            "part": ["func0", "func", "func3", "func_j", "func_s"],
            "type": "rsh_hybrid",
        },
        "wb97x-d3": {
            "tm": None,
            "orca": "wb97x-d3",
            "disp": "included",
            "part": ["func0", "func", "func3", "func_j", "func_s"],
            "type": "rsh_hybrid",
        },
        "wb97x-d3bj": {
            "tm": None,
            "orca": "wb97x-d3bj",
            "disp": "included",
            "part": ["func0", "func", "func3", "func_j", "func_s"],
            "type": "rsh_hybrid",
        },
        "wb97x-d4": {
            "tm": None,
            "orca": "wb97x-d4",
            "disp": "included",
            "part": ["func0", "func", "func3", "func_j", "func_s"],
            "type": "rsh_hybrid",
        },
        "dsd-blyp-d3": {
            "tm": None,
            "orca": "ri-dsd-blyp",
            "disp": "d3bj",
            "part": ["func3", "func_s"],
            "type": "doublehybrid",
        },
    }

    def disp_already_included_in_func(self):
        """return list of all density functionals which inherently include the description
        of (long range) London dispersion"""
        return [
            func
            for func in self.functionals
            if self.functionals[func]["disp"] in ("included", "composite")
        ]

    def dh_dfa(self):
        """return a list of all double hybrid density functionals"""
        return [
            func
            for func in self.functionals
            if self.functionals[func]["type"] in ("doublehybrid")
        ]

    def hybrid_dfa(self):
        """return a list of all hybrid density functionals"""
        return [
            func
            for func in self.functionals
            if self.functionals[func]["type"]
            in ("global_hybrid", "rsh_hybrid", "composite_hybrid")
        ]

    def infos(self, request, prog=None):
        """get information on available functionals for parts and for QM codes"""
        if prog is None and request in (
            "func0",
            "func",
            "func3",
            "func_j",
            "func_s",
            "func_or",
            "func_or_scf",
        ):
            # ("func0_available", "func_available", "fun3_available",
            # "func_j_available", "func_s_available"):
            tmp = []
            for func in self.functionals:
                if request in self.functionals[func].get("part", []):
                    tmp.append(func)
            for func, relay in self.relay_functionals.items():
                if request in self.functionals.get(relay, {}).get("part", []):
                    tmp.append(func)
            return list(set(tmp))
        elif prog in ("tm", "orca") and request in (
            "func0",
            "func",
            "func3",
            "func_j",
            "func_s",
            "func_or",
            "func_or_scf",
        ):
            # list with functional available in either tm or orca
            tmp = []
            for func in self.functionals:
                if (
                    request in self.functionals[func].get("part", [])
                    and self.functionals[func].get(prog, None) is not None
                ):
                    tmp.append(func)
            for func, relay in self.relay_functionals.items():
                if (
                    request in self.functionals.get(relay, {}).get("part", [])
                    and self.functionals.get(relay, {}).get(prog, None) is not None
                ):
                    tmp.append(func)
            return list(set(tmp))
        else:
            print(f"request {request} not known!")


# list with known basissets for warning printing
knownbasissets = [
    "SVP",
    "SV(P)",
    "TZVP",
    "TZVPP",
    "QZVP",
    "QZVPP",
    "def2-SV(P)",
    "def2-mSVP",
    "def2-SVP",
    "def2-TZVP",
    "def2-TZVPP",
    "def2-mTZVP",
    "def2-mTZVPP",
    "def2-TZVPD",
    "def2-SVPD",
    "def-SVP",
    "def-SV(P)",
    "def2-QZVP",
    "DZ",
    "QZV",
    "cc-pVDZ",
    "cc-pVTZ",
    "cc-pVQZ",
    "cc-pV5Z",
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "aug-cc-pVQZ",
    "aug-cc-pV5Z",
    "def2-QZVPP",
    "minix",
    "pcJ-0",
    "pcJ-1",
    "pcJ-2",
    "pcSseg-0",
    "pcSseg-1",
    "pcSseg-2",
    "pcSseg-3",
    "x2c-SVPall-s",
    "x2c-TZVPall-s",
    "def2-TZVP(-f)",  # only ORCA
    "def2-QZVP(-gf)",  # tested for TM
    "def2-TZVPD(-f)",  # tested for TM
]

# program paths:
external_paths = {}
external_paths["orcapath"] = ""
external_paths["orcaversion"] = ""
external_paths["xtbpath"] = ""
external_paths["crestpath"] = ""
external_paths["cosmorssetup"] = ""
external_paths["dbpath"] = ""  # without "DATABASE-COSMO/XXXXXXXXXX"
external_paths["dbpath_normal"] = ""  # with "DATABASE-COSMO/BP-TZVP-COSMO"
external_paths["dbpath_fine"] = ""  # with "DATABASE-COSMO/BP-TZVPD-FINE"
external_paths["cosmothermversion"] = ""
external_paths["mpshiftpath"] = ""
external_paths["escfpath"] = ""
external_paths["cefinepath"] = ""

# information on cosmors parametrizations
cosmors_param = {
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
    "dichloroethane": {
        "cosmors": ["1,2-dichloroethane_c0", "1,2-dichloroethane_c0"],
        "dcosmors": [None, "chcl3"],
        "xtb": [None, "ch2cl2"],
        "cpcm": [None, "CH2Cl2"],
        "smd": ["1,2-DICHLOROETHANE", "1,2-DICHLOROETHANE"],
        "DC": 10.125,
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
    "octane": {
        "cosmors": ["octane_c0", "octane_c0"],
        "dcosmors": [None, "octanol"],
        "xtb": [None, "hexane"],
        "cpcm": [None, "hexane"],
        "smd": ["N-OCTANE", "N-OCTANE"],
        "DC": 1.94,
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
        "xtb": ["ethanol", "ethanol"],
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


class NmrRef:
    """nmrreference data in the format:
    [reference-molecule][func-geometry][funcS][basisS][solvent]
    calculated only for the def2-TZVP basis set (or the respective basis
    set of the composite methods)

    geoopt @ func + basis (def2-TZVP or automatic) + DCOSMO-RS or SMD
    funcS/ + basisS (def2-TZVP or automatic) + DCOSMO-RS or SMD

    h_tm_shieldings
    c_tm_shieldings
    f_tm_shieldings
    si_tm_shieldings
    p_tm_shieldings
    h_orca_shieldings
    c_orca_shieldings
    f_orca_shieldings
    si_orca_shieldings
    p_orca_shieldings
    """

    h_tm_shieldings = {
        "TMS": {
            "r2scan-3c": {
                "r2scan-3c": {
                    "def2-mTZVPP": {
                        "gas": 31.54782223333333,
                        "acetone": 31.507362541666666,
                        "chcl3": 31.51165496666667,
                        "acetonitrile": 31.503900624999996,
                        "ch2cl2": 31.523814674999997,
                        "dmso": 31.506092600000002,
                        "h2o": 31.51254989166667,
                        "methanol": 31.520660391666667,
                        "thf": 31.509407100000004,
                        "toluene": 31.489292383333336,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 32.174366075,
                        "acetone": 32.137970108333334,
                        "chcl3": 32.14006955,
                        "acetonitrile": 32.13512790833333,
                        "ch2cl2": 32.150912858333335,
                        "dmso": 32.137865033333334,
                        "h2o": 32.142147875000006,
                        "methanol": 32.14963604166667,
                        "thf": 32.13958754166667,
                        "toluene": 32.121203466666664,
                    }
                },
                "b97-3c": {
                    "def2-mTZVP": {
                        "gas": 32.174913100000005,
                        "acetone": 32.03599786666667,
                        "chcl3": 32.14600130833333,
                        "acetonitrile": 32.14376344166667,
                        "ch2cl2": 32.15591694166667,
                        "dmso": 32.146624700000004,
                        "h2o": 32.14978633333334,
                        "methanol": 32.15600726666666,
                        "thf": 32.14655236666667,
                        "toluene": 32.128767625,
                    }
                },
                "tpss": {
                    "def2-TZVP": {
                        "gas": 31.898731900000012,
                        "acetone": 31.868752283333336,
                        "chcl3": 31.869320491666667,
                        "acetonitrile": 31.86620341666666,
                        "ch2cl2": 31.87919225,
                        "dmso": 31.870728866666667,
                        "h2o": 31.87253611666667,
                        "methanol": 31.87929998333333,
                        "thf": 31.870970700000004,
                        "toluene": 31.851654483333338,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.66924141666666,
                        "acetone": 31.631933983333337,
                        "chcl3": 31.63644031666667,
                        "acetonitrile": 31.629251491666665,
                        "ch2cl2": 31.64774088333333,
                        "dmso": 31.630312091666667,
                        "h2o": 31.637786341666665,
                        "methanol": 31.644428750000003,
                        "thf": 31.633562058333336,
                        "toluene": 31.6153103,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 31.667160058333337,
                        "acetone": 31.629431858333334,
                        "chcl3": 31.634026916666674,
                        "acetonitrile": 31.627110300000002,
                        "ch2cl2": 31.645328925,
                        "dmso": 31.629352433333338,
                        "h2o": 31.635154825,
                        "methanol": 31.642592158333333,
                        "thf": 31.631994574999997,
                        "toluene": 31.612610516666663,
                    }
                },
            },
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 32.0512048,
                        "acetone": 32.03971003333333,
                        "chcl3": 32.041133316666674,
                        "acetonitrile": 32.03617056666667,
                        "ch2cl2": 32.04777176666666,
                        "dmso": 32.039681316666666,
                        "h2o": 32.036860174999994,
                        "methanol": 32.04573335,
                        "thf": 32.04154705833333,
                        "toluene": 32.02829061666666,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.820450258333327,
                        "acetone": 31.801199816666667,
                        "chcl3": 31.807363400000003,
                        "acetonitrile": 31.797744033333334,
                        "ch2cl2": 31.815502166666665,
                        "dmso": 31.797286500000002,
                        "h2o": 31.801018416666665,
                        "methanol": 31.809920125,
                        "thf": 31.802681225,
                        "toluene": 31.790892416666665,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 32.32369869999999,
                        "acetone": 32.30552229166667,
                        "chcl3": 32.30850654166667,
                        "acetonitrile": 32.3015773,
                        "ch2cl2": 32.31627083333333,
                        "dmso": 32.303862816666665,
                        "h2o": 32.30345545833333,
                        "methanol": 32.3130819,
                        "thf": 32.306951225,
                        "toluene": 32.29417180833333,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 32.099305599999994,
                        "acetone": 32.07685382499999,
                        "chcl3": 32.078372550000005,
                        "acetonitrile": 32.067920741666676,
                        "ch2cl2": 32.0876576,
                        "dmso": 32.07713496666667,
                        "h2o": 32.07222951666666,
                        "methanol": 32.085467083333334,
                        "thf": 32.07950451666667,
                        "toluene": 32.06162065,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.869211950000004,
                        "acetone": 31.83879448333333,
                        "chcl3": 31.845031441666663,
                        "acetonitrile": 31.829924375,
                        "ch2cl2": 31.855811533333338,
                        "dmso": 31.835178675000005,
                        "h2o": 31.83680665833334,
                        "methanol": 31.850090208333338,
                        "thf": 31.841073758333337,
                        "toluene": 31.824697675,
                    }
                },
                "pbeh-3c": {
                    "def2-TZVP": {
                        "gas": 32.37107341666667,
                        "acetone": 32.341934458333334,
                        "chcl3": 32.34503841666666,
                        "acetonitrile": 32.332714675,
                        "ch2cl2": 32.35537393333334,
                        "dmso": 32.34058045833333,
                        "h2o": 32.338073200000004,
                        "methanol": 32.35207416666667,
                        "thf": 32.34418670833334,
                        "toluene": 32.32693729166667,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 31.86774000000001,
                        "acetone": 31.848927016666664,
                        "chcl3": 31.851003891666664,
                        "acetonitrile": 31.843538541666664,
                        "ch2cl2": 31.860415141666664,
                        "dmso": 31.849057266666673,
                        "h2o": 31.844762508333332,
                        "methanol": 31.857667625,
                        "thf": 31.851878716666665,
                        "toluene": 31.833541825,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.636587116666664,
                        "acetone": 31.60924136666667,
                        "chcl3": 31.616506625,
                        "acetonitrile": 31.604173191666664,
                        "ch2cl2": 31.62743169166667,
                        "dmso": 31.604975658333334,
                        "h2o": 31.607992624999994,
                        "methanol": 31.620864658333335,
                        "thf": 31.611675816666665,
                        "toluene": 31.59546233333333,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 32.14311896666666,
                        "acetone": 32.11710325,
                        "chcl3": 32.12106585833333,
                        "acetonitrile": 32.11156126666667,
                        "ch2cl2": 32.1315459,
                        "dmso": 32.114840533333336,
                        "h2o": 32.11376850833333,
                        "methanol": 32.127508733333336,
                        "thf": 32.11950190833333,
                        "toluene": 32.1023676,
                    }
                },
            },
        }
    }
    h_orca_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 32.17000000000001,
                        "acetone": 32.09433333333334,
                        "chcl3": 32.10649999999999,
                        "acetonitrile": 32.09366666666667,
                        "ch2cl2": 32.099,
                        "dmso": 32.09466666666666,
                        "h2o": 32.10341666666666,
                        "methanol": 32.09250000000001,
                        "thf": 32.10183333333333,
                        "toluene": 32.122833333333325,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.819000000000003,
                        "acetone": 31.732666666666663,
                        "chcl3": 31.747000000000003,
                        "acetonitrile": 31.73166666666667,
                        "ch2cl2": 31.738416666666666,
                        "dmso": 31.732666666666663,
                        "h2o": 31.741500000000002,
                        "methanol": 31.73066666666666,
                        "thf": 31.74116666666667,
                        "toluene": 31.765999999999995,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 31.91416666666667,
                        "acetone": 31.83541666666667,
                        "chcl3": 31.84766666666667,
                        "acetonitrile": 31.834666666666667,
                        "ch2cl2": 31.839916666666667,
                        "dmso": 31.835583333333332,
                        "h2o": 31.844166666666666,
                        "methanol": 31.833166666666667,
                        "thf": 31.842583333333334,
                        "toluene": 31.86475,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 31.952,
                        "acetone": 31.867499999999996,
                        "chcl3": 31.880999999999997,
                        "acetonitrile": 31.866666666666664,
                        "ch2cl2": 31.872666666666664,
                        "dmso": 31.86758333333333,
                        "h2o": 31.876083333333337,
                        "methanol": 31.86533333333333,
                        "thf": 31.8755,
                        "toluene": 31.89966666666666,
                    }
                },
                "pbeh-3c": {
                    "def2-TZVP": {
                        "gas": 32.324999999999996,
                        "acetone": 32.23866666666667,
                        "chcl3": 32.25299999999999,
                        "acetonitrile": 32.23783333333333,
                        "ch2cl2": 32.24466666666667,
                        "dmso": 32.23866666666667,
                        "h2o": 32.24733333333333,
                        "methanol": 32.23666666666667,
                        "thf": 32.24733333333333,
                        "toluene": 32.272,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 31.817999999999998,
                        "acetone": 31.73233333333333,
                        "chcl3": 31.746333333333336,
                        "acetonitrile": 31.73133333333333,
                        "ch2cl2": 31.737666666666666,
                        "dmso": 31.73233333333333,
                        "h2o": 31.740666666666666,
                        "methanol": 31.73,
                        "thf": 31.740499999999994,
                        "toluene": 31.765666666666664,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 32.21800000000001,
                        "acetone": 32.140166666666666,
                        "chcl3": 32.152166666666666,
                        "acetonitrile": 32.140499999999996,
                        "ch2cl2": 32.145,
                        "dmso": 32.14183333333333,
                        "h2o": 32.175000000000004,
                        "methanol": 32.13766666666667,
                        "thf": 32.148,
                        "toluene": 32.168833333333325,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.868,
                        "acetone": 31.778999999999996,
                        "chcl3": 31.792583333333337,
                        "acetonitrile": 31.778666666666663,
                        "ch2cl2": 31.784333333333336,
                        "dmso": 31.78016666666667,
                        "h2o": 31.815166666666666,
                        "methanol": 31.77633333333333,
                        "thf": 31.787500000000005,
                        "toluene": 31.812,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 31.962999999999997,
                        "acetone": 31.881250000000005,
                        "chcl3": 31.89325,
                        "acetonitrile": 31.881583333333335,
                        "ch2cl2": 31.886000000000006,
                        "dmso": 31.882583333333333,
                        "h2o": 31.916833333333333,
                        "methanol": 31.878500000000003,
                        "thf": 31.889,
                        "toluene": 31.910750000000004,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 32.00091666666666,
                        "acetone": 31.913416666666663,
                        "chcl3": 31.9265,
                        "acetonitrile": 31.9135,
                        "ch2cl2": 31.918499999999995,
                        "dmso": 31.914666666666665,
                        "h2o": 31.94883333333333,
                        "methanol": 31.910666666666668,
                        "thf": 31.921500000000005,
                        "toluene": 31.94516666666667,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 32.373,
                        "acetone": 32.28366666666667,
                        "chcl3": 32.29716666666666,
                        "acetonitrile": 32.28333333333333,
                        "ch2cl2": 32.288666666666664,
                        "dmso": 32.284499999999994,
                        "h2o": 32.317166666666665,
                        "methanol": 32.28066666666667,
                        "thf": 32.29183333333334,
                        "toluene": 32.31616666666667,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 31.868,
                        "acetone": 31.778666666666663,
                        "chcl3": 31.792500000000004,
                        "acetonitrile": 31.778666666666663,
                        "ch2cl2": 31.784333333333336,
                        "dmso": 31.78033333333333,
                        "h2o": 31.794583333333332,
                        "methanol": 31.77633333333333,
                        "thf": 31.787500000000005,
                        "toluene": 31.812,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 31.97300000000001,
                        "acetone": 31.898,
                        "chcl3": 31.909500000000005,
                        "acetonitrile": 31.897833333333338,
                        "ch2cl2": 31.902666666666665,
                        "dmso": 31.898999999999997,
                        "h2o": 31.910666666666668,
                        "methanol": 31.89566666666667,
                        "thf": 31.90516666666667,
                        "toluene": 31.925,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 31.625,
                        "acetone": 31.537166666666668,
                        "chcl3": 31.550499999999996,
                        "acetonitrile": 31.536666666666665,
                        "ch2cl2": 31.542500000000004,
                        "dmso": 31.537666666666667,
                        "h2o": 31.549500000000005,
                        "methanol": 31.53458333333334,
                        "thf": 31.545499999999993,
                        "toluene": 31.569,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 31.718000000000004,
                        "acetone": 31.639666666666667,
                        "chcl3": 31.651416666666663,
                        "acetonitrile": 31.639499999999998,
                        "ch2cl2": 31.644083333333338,
                        "dmso": 31.640416666666667,
                        "h2o": 31.65216666666667,
                        "methanol": 31.636916666666664,
                        "thf": 31.64683333333333,
                        "toluene": 31.667833333333334,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 31.757,
                        "acetone": 31.672250000000002,
                        "chcl3": 31.68516666666667,
                        "acetonitrile": 31.67166666666667,
                        "ch2cl2": 31.6775,
                        "dmso": 31.67266666666666,
                        "h2o": 31.68466666666666,
                        "methanol": 31.66966666666667,
                        "thf": 31.680166666666665,
                        "toluene": 31.703,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 32.13400000000001,
                        "acetone": 32.047333333333334,
                        "chcl3": 32.06066666666667,
                        "acetonitrile": 32.04666666666666,
                        "ch2cl2": 32.05266666666666,
                        "dmso": 32.047666666666665,
                        "h2o": 32.059,
                        "methanol": 32.044666666666664,
                        "thf": 32.05566666666666,
                        "toluene": 32.079,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 31.622999999999994,
                        "acetone": 31.536666666666665,
                        "chcl3": 31.55,
                        "acetonitrile": 31.5365,
                        "ch2cl2": 31.54183333333333,
                        "dmso": 31.537666666666667,
                        "h2o": 31.548666666666666,
                        "methanol": 31.533833333333334,
                        "thf": 31.544833333333333,
                        "toluene": 31.56866666666667,
                    }
                },
            },
        }
    }
    c_tm_shieldings = {
        "TMS": {
            "r2scan-3c": {
                "r2scan-3c": {
                    "def2-mTZVPP": {
                        "gas": 194.0424432,
                        "acetone": 194.63496962500002,
                        "chcl3": 194.51409719999998,
                        "acetonitrile": 194.719151575,
                        "ch2cl2": 194.35334055,
                        "dmso": 194.712879125,
                        "h2o": 194.73630830000002,
                        "methanol": 194.47241222500003,
                        "thf": 194.5933908,
                        "toluene": 194.76792010000003,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 197.51118245,
                        "acetone": 198.07638194999998,
                        "chcl3": 197.9935314,
                        "acetonitrile": 198.194965175,
                        "ch2cl2": 197.81458097499998,
                        "dmso": 198.137269575,
                        "h2o": 198.22405855,
                        "methanol": 197.911550875,
                        "thf": 198.02789195,
                        "toluene": 198.28027294999998,
                    }
                },
                "b97-3c": {
                    "def2-mTZVP": {
                        "gas": 184.57099812500002,
                        "acetone": 183.56356677499997,
                        "chcl3": 185.08358385,
                        "acetonitrile": 185.26145292500001,
                        "ch2cl2": 184.92964515000003,
                        "dmso": 185.1924247,
                        "h2o": 185.317873325,
                        "methanol": 185.01346424999997,
                        "thf": 185.100088675,
                        "toluene": 185.3094498,
                    }
                },
                "tpss": {
                    "def2-TZVP": {
                        "gas": 185.706348925,
                        "acetone": 186.20184215,
                        "chcl3": 186.13177655,
                        "acetonitrile": 186.29212465,
                        "ch2cl2": 186.000285975,
                        "dmso": 186.2141097,
                        "h2o": 186.333736775,
                        "methanol": 186.07168465,
                        "thf": 186.13180945,
                        "toluene": 186.32699157500002,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 187.894285475,
                        "acetone": 188.51151167499998,
                        "chcl3": 188.36626362500002,
                        "acetonitrile": 188.58446245,
                        "ch2cl2": 188.21852824999996,
                        "dmso": 188.60960655,
                        "h2o": 188.5966455,
                        "methanol": 188.355007975,
                        "thf": 188.4689729,
                        "toluene": 188.594034275,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 189.78494644999998,
                        "acetone": 190.329502875,
                        "chcl3": 190.204013175,
                        "acetonitrile": 190.397052075,
                        "ch2cl2": 190.06665505,
                        "dmso": 190.4107424,
                        "h2o": 190.40970589999998,
                        "methanol": 190.188391875,
                        "thf": 190.2872299,
                        "toluene": 190.41299607500002,
                    }
                },
            },
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 186.6465687,
                        "acetone": 187.27903107499998,
                        "chcl3": 187.238498325,
                        "acetonitrile": 187.372512775,
                        "ch2cl2": 187.0771589,
                        "dmso": 187.243299225,
                        "h2o": 187.37157565,
                        "methanol": 187.10988087500002,
                        "thf": 187.19458635,
                        "toluene": 187.48276635,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 188.859355325,
                        "acetone": 189.6196798,
                        "chcl3": 189.4971041,
                        "acetonitrile": 189.698041075,
                        "ch2cl2": 189.318608125,
                        "dmso": 189.68253637499998,
                        "h2o": 189.65553119999998,
                        "methanol": 189.409198575,
                        "thf": 189.55889105,
                        "toluene": 189.776394325,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 198.41611147499998,
                        "acetone": 199.13367970000002,
                        "chcl3": 199.054179875,
                        "acetonitrile": 199.250248325,
                        "ch2cl2": 198.845265825,
                        "dmso": 199.185056825,
                        "h2o": 199.2289907,
                        "methanol": 198.917945675,
                        "thf": 199.076003325,
                        "toluene": 199.3931504,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 186.97419324999998,
                        "acetone": 187.496073025,
                        "chcl3": 187.45393565,
                        "acetonitrile": 187.554538075,
                        "ch2cl2": 187.31238564999998,
                        "dmso": 187.469466275,
                        "h2o": 187.57139320000002,
                        "methanol": 187.344972675,
                        "thf": 187.42200885,
                        "toluene": 187.671731225,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 189.169130675,
                        "acetone": 189.816064175,
                        "chcl3": 189.69082477499998,
                        "acetonitrile": 189.860330875,
                        "ch2cl2": 189.532330975,
                        "dmso": 189.88587445000002,
                        "h2o": 189.8368566,
                        "methanol": 189.62332455,
                        "thf": 189.76569125,
                        "toluene": 189.94371412499999,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 198.7168509,
                        "acetone": 199.3308802,
                        "chcl3": 199.25125382500002,
                        "acetonitrile": 199.41320919999998,
                        "ch2cl2": 199.06108425,
                        "dmso": 199.390014125,
                        "h2o": 199.41478467500002,
                        "methanol": 199.13192775,
                        "thf": 199.28161922500001,
                        "toluene": 199.562540575,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 185.410099625,
                        "acetone": 185.99193982499997,
                        "chcl3": 185.949648475,
                        "acetonitrile": 186.0799505,
                        "ch2cl2": 185.80363820000002,
                        "dmso": 185.97415155,
                        "h2o": 186.07484635,
                        "methanol": 185.839592875,
                        "thf": 185.9190184,
                        "toluene": 186.17204557500003,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 187.626469575,
                        "acetone": 188.34549135,
                        "chcl3": 188.212218325,
                        "acetonitrile": 188.413268225,
                        "ch2cl2": 188.04820440000003,
                        "dmso": 188.42875420000001,
                        "h2o": 188.3724699,
                        "methanol": 188.14698049999998,
                        "thf": 188.2963985,
                        "toluene": 188.46803717499998,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 197.27823677499998,
                        "acetone": 197.953274625,
                        "chcl3": 197.871683925,
                        "acetonitrile": 198.0615831,
                        "ch2cl2": 197.6764831,
                        "dmso": 198.014841225,
                        "h2o": 198.048432475,
                        "methanol": 197.75143105,
                        "thf": 197.905333025,
                        "toluene": 198.186480775,
                    }
                },
            },
        }
    }
    c_orca_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 188.604,
                        "acetone": 189.7395,
                        "chcl3": 189.5435,
                        "acetonitrile": 189.77,
                        "ch2cl2": 189.6625,
                        "dmso": 189.8015,
                        "h2o": 189.8495,
                        "methanol": 189.77,
                        "thf": 189.647,
                        "toluene": 189.30400000000003,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 188.867,
                        "acetone": 190.265,
                        "chcl3": 190.02224999999999,
                        "acetonitrile": 190.298,
                        "ch2cl2": 190.16649999999998,
                        "dmso": 190.33175,
                        "h2o": 190.38799999999998,
                        "methanol": 190.29875,
                        "thf": 190.1445,
                        "toluene": 189.73375,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 191.37099999999998,
                        "acetone": 192.606,
                        "chcl3": 192.385,
                        "acetonitrile": 192.63599999999997,
                        "ch2cl2": 192.51575000000003,
                        "dmso": 192.66625000000002,
                        "h2o": 192.7205,
                        "methanol": 192.63524999999998,
                        "thf": 192.4955,
                        "toluene": 192.12275,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 190.36075,
                        "acetone": 191.689,
                        "chcl3": 191.453,
                        "acetonitrile": 191.72175000000001,
                        "ch2cl2": 191.5935,
                        "dmso": 191.753,
                        "h2o": 191.8085,
                        "methanol": 191.72150000000002,
                        "thf": 191.57150000000001,
                        "toluene": 191.17225,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 198.458,
                        "acetone": 199.905,
                        "chcl3": 199.649,
                        "acetonitrile": 199.94,
                        "ch2cl2": 199.8025,
                        "dmso": 199.9715,
                        "h2o": 200.0265,
                        "methanol": 199.93900,
                        "thf": 199.7775,
                        "toluene": 199.3395,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 190.719,
                        "acetone": 191.988,
                        "chcl3": 191.7645,
                        "acetonitrile": 192.019,
                        "ch2cl2": 191.8965,
                        "dmso": 192.05150000000003,
                        "h2o": 192.1055,
                        "methanol": 192.02,
                        "thf": 191.8775,
                        "toluene": 191.4905,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 188.908,
                        "acetone": 190.0265,
                        "chcl3": 189.83749999999998,
                        "acetonitrile": 190.062,
                        "ch2cl2": 189.954,
                        "dmso": 190.103,
                        "h2o": 190.07774999999998,
                        "methanol": 190.0595,
                        "thf": 189.9445,
                        "toluene": 189.614,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 189.18025,
                        "acetone": 190.57025000000002,
                        "chcl3": 190.33075,
                        "acetonitrile": 190.60525,
                        "ch2cl2": 190.47,
                        "dmso": 190.65175,
                        "h2o": 190.59925000000004,
                        "methanol": 190.60775,
                        "thf": 190.456,
                        "toluene": 190.058,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 191.66199999999998,
                        "acetone": 192.88025,
                        "chcl3": 192.66174999999998,
                        "acetonitrile": 192.915,
                        "ch2cl2": 192.79025,
                        "dmso": 192.95425,
                        "h2o": 192.91275000000002,
                        "methanol": 192.91250000000002,
                        "thf": 192.77625,
                        "toluene": 192.4135,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 190.65525,
                        "acetone": 191.97199999999998,
                        "chcl3": 191.73825,
                        "acetonitrile": 192.00725,
                        "ch2cl2": 191.875,
                        "dmso": 192.04950000000002,
                        "h2o": 191.99675000000002,
                        "methanol": 192.007,
                        "thf": 191.86025,
                        "toluene": 191.47125,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 198.752,
                        "acetone": 200.196,
                        "chcl3": 199.9445,
                        "acetonitrile": 200.23250000000002,
                        "ch2cl2": 200.0925,
                        "dmso": 200.277,
                        "h2o": 200.15925,
                        "methanol": 200.23350000000002,
                        "thf": 200.075,
                        "toluene": 199.65050000000002,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 191.037,
                        "acetone": 192.29649999999998,
                        "chcl3": 192.0765,
                        "acetonitrile": 192.3275,
                        "ch2cl2": 192.20350000000002,
                        "dmso": 192.3755,
                        "h2o": 192.188,
                        "methanol": 192.33275,
                        "thf": 192.1925,
                        "toluene": 191.8175,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 187.22,
                        "acetone": 188.442,
                        "chcl3": 188.214,
                        "acetonitrile": 188.4745,
                        "ch2cl2": 188.351,
                        "dmso": 188.5115,
                        "h2o": 188.58350000000002,
                        "methanol": 188.473,
                        "thf": 188.33950000000002,
                        "toluene": 187.965,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 187.5725,
                        "acetone": 188.99225,
                        "chcl3": 188.73424999999997,
                        "acetonitrile": 189.0295,
                        "ch2cl2": 188.8875,
                        "dmso": 189.06875,
                        "h2o": 189.14175,
                        "methanol": 189.0275,
                        "thf": 188.8665,
                        "toluene": 188.4305,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 190.06825,
                        "acetone": 191.39,
                        "chcl3": 191.15425,
                        "acetonitrile": 191.42600000000002,
                        "ch2cl2": 191.29475000000002,
                        "dmso": 191.461,
                        "h2o": 191.53225,
                        "methanol": 191.4225,
                        "thf": 191.27499999999998,
                        "toluene": 190.87675000000002,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 189.04575,
                        "acetone": 190.45225000000002,
                        "chcl3": 190.20074999999997,
                        "acetonitrile": 190.4885,
                        "ch2cl2": 190.35025000000002,
                        "dmso": 190.52525,
                        "h2o": 190.5975,
                        "methanol": 190.4855,
                        "thf": 190.32899999999998,
                        "toluene": 189.904,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 197.184,
                        "acetone": 198.7195,
                        "chcl3": 198.449,
                        "acetonitrile": 198.75799999999998,
                        "ch2cl2": 198.611,
                        "dmso": 198.7955,
                        "h2o": 198.8655,
                        "methanol": 198.755,
                        "thf": 198.587,
                        "toluene": 198.1245,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 189.386,
                        "acetone": 190.7245,
                        "chcl3": 190.488,
                        "acetonitrile": 190.7585,
                        "ch2cl2": 190.6275,
                        "dmso": 190.7975,
                        "h2o": 190.87900000000002,
                        "methanol": 190.75799999999998,
                        "thf": 190.6095,
                        "toluene": 190.2095,
                    }
                },
            },
        }
    }
    f_tm_shieldings = {
        "CFCl3": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 163.5665883,
                        "acetone": 165.9168679,
                        "chcl3": 165.043061,
                        "acetonitrile": 166.377831,
                        "ch2cl2": 164.776383,
                        "dmso": 166.1839641,
                        "h2o": 166.880495,
                        "methanol": 165.4364879,
                        "thf": 165.7384153,
                        "toluene": 165.7812123,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 179.4820255,
                        "acetone": 181.9743764,
                        "chcl3": 181.1338758,
                        "acetonitrile": 182.4438224,
                        "ch2cl2": 180.8751895,
                        "dmso": 182.2224636,
                        "h2o": 182.9958356,
                        "methanol": 181.5031528,
                        "thf": 181.7669891,
                        "toluene": 181.7963177,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 225.045234,
                        "acetone": 226.6335916,
                        "chcl3": 226.0133192,
                        "acetonitrile": 226.9371636,
                        "ch2cl2": 225.8300352,
                        "dmso": 226.8061873,
                        "h2o": 227.4000142,
                        "methanol": 226.3012569,
                        "thf": 226.5247654,
                        "toluene": 226.555523,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 150.4514566,
                        "acetone": 151.5612999,
                        "chcl3": 150.5819485,
                        "acetonitrile": 151.9884593,
                        "ch2cl2": 150.2953968,
                        "dmso": 151.8818575,
                        "h2o": 151.6179136,
                        "methanol": 151.0439011,
                        "thf": 151.4207377,
                        "toluene": 151.4686522,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 167.7783433,
                        "acetone": 169.09491,
                        "chcl3": 168.1354478,
                        "acetonitrile": 169.5416871,
                        "ch2cl2": 167.8558489,
                        "dmso": 169.3950732,
                        "h2o": 169.2178304,
                        "methanol": 168.5860848,
                        "thf": 168.9136991,
                        "toluene": 168.9347931,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 213.6651892,
                        "acetone": 214.1284506,
                        "chcl3": 213.4293417,
                        "acetonitrile": 214.4297108,
                        "ch2cl2": 213.2298905,
                        "dmso": 214.366451,
                        "h2o": 214.1162368,
                        "methanol": 213.76845,
                        "thf": 214.0512078,
                        "toluene": 214.0924969,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 146.4091676,
                        "acetone": 148.7113398,
                        "chcl3": 147.7715256,
                        "acetonitrile": 149.1854535,
                        "ch2cl2": 147.4708159,
                        "dmso": 148.9781692,
                        "h2o": 148.8407317,
                        "methanol": 148.1815132,
                        "thf": 148.5140784,
                        "toluene": 148.6001306,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 163.4654205,
                        "acetone": 165.9356023,
                        "chcl3": 165.0269644,
                        "acetonitrile": 166.4188044,
                        "ch2cl2": 164.7336009,
                        "dmso": 166.1830401,
                        "h2o": 166.0858984,
                        "methanol": 165.4145633,
                        "thf": 165.7038144,
                        "toluene": 165.7726604,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 209.8752809,
                        "acetone": 211.4025693,
                        "chcl3": 210.7286529,
                        "acetonitrile": 211.7120494,
                        "ch2cl2": 210.5166504,
                        "dmso": 211.5990015,
                        "h2o": 211.4250312,
                        "methanol": 211.0321396,
                        "thf": 211.2798891,
                        "toluene": 211.3499046,
                    }
                },
            },
        }
    }
    f_orca_shieldings = {
        "CFCl3": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 166.028,
                        "acetone": 167.858,
                        "chcl3": 167.569,
                        "acetonitrile": 167.92,
                        "ch2cl2": 167.732,
                        "dmso": 167.992,
                        "h2o": 168.239,
                        "methanol": 167.889,
                        "thf": 167.737,
                        "toluene": 167.278,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 178.99,
                        "acetone": 181.086,
                        "chcl3": 180.741,
                        "acetonitrile": 181.154,
                        "ch2cl2": 180.939,
                        "dmso": 181.224,
                        "h2o": 181.464,
                        "methanol": 181.123,
                        "thf": 180.934,
                        "toluene": 180.377,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 225.542,
                        "acetone": 227.877,
                        "chcl3": 227.478,
                        "acetonitrile": 227.949,
                        "ch2cl2": 227.712,
                        "dmso": 228.007,
                        "h2o": 228.213,
                        "methanol": 227.919,
                        "thf": 227.691,
                        "toluene": 227.033,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 193.433,
                        "acetone": 195.381,
                        "chcl3": 195.059,
                        "acetonitrile": 195.445,
                        "ch2cl2": 195.245,
                        "dmso": 195.508,
                        "h2o": 195.733,
                        "methanol": 195.415,
                        "thf": 195.239,
                        "toluene": 194.719,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 224.834,
                        "acetone": 226.308,
                        "chcl3": 226.076,
                        "acetonitrile": 226.36,
                        "ch2cl2": 226.207,
                        "dmso": 226.424,
                        "h2o": 226.639,
                        "methanol": 226.333,
                        "thf": 226.215,
                        "toluene": 225.843,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 144.178,
                        "acetone": 146.15,
                        "chcl3": 145.821,
                        "acetonitrile": 146.219,
                        "ch2cl2": 146.007,
                        "dmso": 146.298,
                        "h2o": 146.569,
                        "methanol": 146.185,
                        "thf": 146.012,
                        "toluene": 145.488,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 153.325,
                        "acetone": 153.259,
                        "chcl3": 152.987,
                        "acetonitrile": 153.326,
                        "ch2cl2": 153.137,
                        "dmso": 153.425,
                        "h2o": 153.729,
                        "methanol": 153.292,
                        "thf": 153.16,
                        "toluene": 152.731,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 167.245,
                        "acetone": 167.447,
                        "chcl3": 167.121,
                        "acetonitrile": 167.52,
                        "ch2cl2": 167.31,
                        "dmso": 167.626,
                        "h2o": 167.92,
                        "methanol": 167.486,
                        "thf": 167.322,
                        "toluene": 166.785,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 216.287,
                        "acetone": 217.144,
                        "chcl3": 216.726,
                        "acetonitrile": 217.223,
                        "ch2cl2": 216.969,
                        "dmso": 217.304,
                        "h2o": 217.555,
                        "methanol": 217.19,
                        "thf": 216.957,
                        "toluene": 216.272,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 182.767,
                        "acetone": 182.921,
                        "chcl3": 182.602,
                        "acetonitrile": 182.99,
                        "ch2cl2": 182.783,
                        "dmso": 183.077,
                        "h2o": 183.351,
                        "methanol": 182.957,
                        "thf": 182.792,
                        "toluene": 182.279,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 213.421,
                        "acetone": 213.215,
                        "chcl3": 212.997,
                        "acetonitrile": 213.271,
                        "ch2cl2": 213.116,
                        "dmso": 213.36,
                        "h2o": 213.627,
                        "methanol": 213.241,
                        "thf": 213.14,
                        "toluene": 212.796,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 130.539,
                        "acetone": 130.291,
                        "chcl3": 130.081,
                        "acetonitrile": 130.364,
                        "ch2cl2": 130.242,
                        "dmso": 130.472,
                        "h2o": 130.803,
                        "methanol": 130.326,
                        "thf": 130.267,
                        "toluene": 129.808,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 148.387,
                        "acetone": 149.573,
                        "chcl3": 149.247,
                        "acetonitrile": 149.647,
                        "ch2cl2": 149.43,
                        "dmso": 149.748,
                        "h2o": 150.066,
                        "methanol": 149.609,
                        "thf": 149.446,
                        "toluene": 148.927,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 162.075,
                        "acetone": 163.638,
                        "chcl3": 163.239,
                        "acetonitrile": 163.71,
                        "ch2cl2": 163.472,
                        "dmso": 163.807,
                        "h2o": 164.125,
                        "methanol": 163.671,
                        "thf": 163.476,
                        "toluene": 162.835,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 211.635,
                        "acetone": 213.66,
                        "chcl3": 213.199,
                        "acetonitrile": 213.746,
                        "ch2cl2": 213.469,
                        "dmso": 213.828,
                        "h2o": 214.092,
                        "methanol": 213.71,
                        "thf": 213.451,
                        "toluene": 212.692,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 177.986,
                        "acetone": 179.452,
                        "chcl3": 179.093,
                        "acetonitrile": 179.528,
                        "ch2cl2": 179.299,
                        "dmso": 179.616,
                        "h2o": 179.902,
                        "methanol": 179.491,
                        "thf": 179.302,
                        "toluene": 178.721,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 208.73,
                        "acetone": 209.687,
                        "chcl3": 209.429,
                        "acetonitrile": 209.749,
                        "ch2cl2": 209.573,
                        "dmso": 209.825,
                        "h2o": 210.102,
                        "methanol": 209.716,
                        "thf": 209.592,
                        "toluene": 209.176,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 124.897,
                        "acetone": 126.154,
                        "chcl3": 125.806,
                        "acetonitrile": 126.235,
                        "ch2cl2": 126.001,
                        "dmso": 126.345,
                        "h2o": 126.689,
                        "methanol": 126.193,
                        "thf": 126.019,
                        "toluene": 125.465,
                    }
                },
            },
        }
    }
    p_tm_shieldings = {
        "PH3": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 560.9783608,
                        "acetone": 559.5567974,
                        "chcl3": 555.7297268,
                        "acetonitrile": 558.7420853,
                        "ch2cl2": 555.9207578,
                        "dmso": 559.0317956,
                        "h2o": 551.9868157,
                        "methanol": 557.7229598,
                        "thf": 559.4070044,
                        "toluene": 558.9538264,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 573.7889709,
                        "acetone": 572.6807308,
                        "chcl3": 568.6200619,
                        "acetonitrile": 572.0156003,
                        "ch2cl2": 568.6775273,
                        "dmso": 572.2984368,
                        "h2o": 564.8512663,
                        "methanol": 570.6948985,
                        "thf": 572.4491708,
                        "toluene": 572.2945282,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 622.6149401,
                        "acetone": 624.221383,
                        "chcl3": 622.2460822,
                        "acetonitrile": 624.0839458,
                        "ch2cl2": 622.3660073,
                        "dmso": 623.8685076,
                        "h2o": 622.54767,
                        "methanol": 623.1569748,
                        "thf": 623.7253948,
                        "toluene": 623.2733775,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 559.5296772,
                        "acetone": 557.5438599,
                        "chcl3": 553.7653249,
                        "acetonitrile": 556.735552,
                        "ch2cl2": 554.1613395,
                        "dmso": 557.010476,
                        "h2o": 550.1185847,
                        "methanol": 555.82703,
                        "thf": 557.2207586,
                        "toluene": 556.8427805,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 572.4232552,
                        "acetone": 570.7398164,
                        "chcl3": 566.7271447,
                        "acetonitrile": 570.0779914,
                        "ch2cl2": 566.9826221,
                        "dmso": 570.3456887,
                        "h2o": 563.05667,
                        "methanol": 568.8622417,
                        "thf": 570.3305746,
                        "toluene": 570.2507738,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 621.2286124,
                        "acetone": 622.356702,
                        "chcl3": 620.3365742,
                        "acetonitrile": 622.2263079,
                        "ch2cl2": 620.6570087,
                        "dmso": 621.9912341,
                        "h2o": 620.7021951,
                        "methanol": 621.3567408,
                        "thf": 621.7091401,
                        "toluene": 621.3088355,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 558.1589032,
                        "acetone": 556.5475548,
                        "chcl3": 553.3273579,
                        "acetonitrile": 555.6559443,
                        "ch2cl2": 553.600544,
                        "dmso": 556.0983125,
                        "h2o": 548.970911,
                        "methanol": 555.4535832,
                        "thf": 556.3191064,
                        "toluene": 555.9299261,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 571.012794,
                        "acetone": 569.7250563,
                        "chcl3": 566.2936179,
                        "acetonitrile": 568.9923465,
                        "ch2cl2": 566.4237381,
                        "dmso": 569.4236946,
                        "h2o": 561.898531,
                        "methanol": 568.4989088,
                        "thf": 569.4140377,
                        "toluene": 569.3191735,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 620.0674752,
                        "acetone": 621.5116584,
                        "chcl3": 619.9397925,
                        "acetonitrile": 621.2898165,
                        "ch2cl2": 620.15928,
                        "dmso": 621.2154327,
                        "h2o": 619.7280828,
                        "methanol": 621.0126668,
                        "thf": 620.9449236,
                        "toluene": 620.5363442,
                    }
                },
            },
        },
        "TMP": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 281.6302978,
                        "acetone": 265.4354914,
                        "chcl3": 257.5409613,
                        "acetonitrile": 263.2430698,
                        "ch2cl2": 257.0543221,
                        "dmso": 262.8752182,
                        "h2o": 242.4838211,
                        "methanol": 245.6431135,
                        "thf": 266.7188352,
                        "toluene": 269.0597797,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 277.8252556,
                        "acetone": 261.5502528,
                        "chcl3": 254.1109855,
                        "acetonitrile": 259.5059377,
                        "ch2cl2": 253.6358478,
                        "dmso": 258.7821425,
                        "h2o": 239.5329333,
                        "methanol": 242.1687948,
                        "thf": 262.8378646,
                        "toluene": 265.4050199,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 390.6073841,
                        "acetone": 378.6668397,
                        "chcl3": 373.2000393,
                        "acetonitrile": 377.1343123,
                        "ch2cl2": 372.9163524,
                        "dmso": 376.6203422,
                        "h2o": 362.7163813,
                        "methanol": 364.8220379,
                        "thf": 379.5051748,
                        "toluene": 381.2789752,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 276.8654211,
                        "acetone": 259.8829696,
                        "chcl3": 251.5648819,
                        "acetonitrile": 257.7225804,
                        "ch2cl2": 251.0880934,
                        "dmso": 256.90761,
                        "h2o": 234.4800595,
                        "methanol": 237.4630709,
                        "thf": 261.291204,
                        "toluene": 263.9827571,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 273.0911933,
                        "acetone": 256.1507446,
                        "chcl3": 248.2072561,
                        "acetonitrile": 254.0571117,
                        "ch2cl2": 247.7513367,
                        "dmso": 253.0100842,
                        "h2o": 231.7425518,
                        "methanol": 234.1695454,
                        "thf": 257.4644157,
                        "toluene": 260.3717755,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 386.2437698,
                        "acetone": 373.8145109,
                        "chcl3": 368.1719462,
                        "acetonitrile": 372.350904,
                        "ch2cl2": 367.8934403,
                        "dmso": 371.4995766,
                        "h2o": 355.9965281,
                        "methanol": 358.0517851,
                        "thf": 374.7716841,
                        "toluene": 376.8283779,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 278.0447826,
                        "acetone": 261.4382678,
                        "chcl3": 253.5317417,
                        "acetonitrile": 259.5831076,
                        "ch2cl2": 253.0735218,
                        "dmso": 258.8205488,
                        "h2o": 236.9938311,
                        "methanol": 240.0596152,
                        "thf": 262.646474,
                        "toluene": 265.5482099,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 274.1582231,
                        "acetone": 257.5976215,
                        "chcl3": 250.0455696,
                        "acetonitrile": 255.8739799,
                        "ch2cl2": 249.6032437,
                        "dmso": 254.7109046,
                        "h2o": 234.1066151,
                        "methanol": 236.6658834,
                        "thf": 258.6914971,
                        "toluene": 261.8410368,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 387.4697022,
                        "acetone": 375.2569197,
                        "chcl3": 369.9533245,
                        "acetonitrile": 374.0256406,
                        "ch2cl2": 369.6688695,
                        "dmso": 373.1520781,
                        "h2o": 358.1827766,
                        "methanol": 360.3168296,
                        "thf": 376.0015788,
                        "toluene": 378.3153047,
                    }
                },
            },
        },
    }
    p_orca_shieldings = {
        "PH3": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 578.49,
                        "acetone": 577.53,
                        "chcl3": 577.773,
                        "acetonitrile": 577.631,
                        "ch2cl2": 577.63,
                        "dmso": 577.688,
                        "h2o": 577.764,
                        "methanol": 577.506,
                        "thf": 577.671,
                        "toluene": 577.946,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 573.639,
                        "acetone": 573.637,
                        "chcl3": 573.71,
                        "acetonitrile": 573.764,
                        "ch2cl2": 573.67,
                        "dmso": 573.829,
                        "h2o": 573.914,
                        "methanol": 573.632,
                        "thf": 573.688,
                        "toluene": 573.665,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 569.431,
                        "acetone": 567.575,
                        "chcl3": 567.994,
                        "acetonitrile": 567.65,
                        "ch2cl2": 567.746,
                        "dmso": 567.695,
                        "h2o": 567.745,
                        "methanol": 567.531,
                        "thf": 567.809,
                        "toluene": 568.372,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 568.27,
                        "acetone": 568.185,
                        "chcl3": 568.261,
                        "acetonitrile": 568.31,
                        "ch2cl2": 568.218,
                        "dmso": 568.375,
                        "h2o": 568.459,
                        "methanol": 568.18,
                        "thf": 568.236,
                        "toluene": 568.231,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 622.505,
                        "acetone": 626.377,
                        "chcl3": 625.536,
                        "acetonitrile": 626.609,
                        "ch2cl2": 626.042,
                        "dmso": 626.709,
                        "h2o": 626.85,
                        "methanol": 626.48,
                        "thf": 625.933,
                        "toluene": 624.513,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 587.254,
                        "acetone": 587.821,
                        "chcl3": 587.78,
                        "acetonitrile": 587.962,
                        "ch2cl2": 587.81,
                        "dmso": 588.032,
                        "h2o": 588.129,
                        "methanol": 587.829,
                        "thf": 587.812,
                        "toluene": 587.606,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 574.673,
                        "acetone": 575.587,
                        "chcl3": 575.672,
                        "acetonitrile": 575.6,
                        "ch2cl2": 575.619,
                        "dmso": 575.662,
                        "h2o": 575.948,
                        "methanol": 575.57,
                        "thf": 575.668,
                        "toluene": 575.8,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 569.721,
                        "acetone": 571.667,
                        "chcl3": 571.577,
                        "acetonitrile": 571.703,
                        "ch2cl2": 571.631,
                        "dmso": 571.774,
                        "h2o": 572.075,
                        "methanol": 571.67,
                        "thf": 571.656,
                        "toluene": 571.48,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 565.936,
                        "acetone": 565.88,
                        "chcl3": 566.179,
                        "acetonitrile": 565.866,
                        "ch2cl2": 566.012,
                        "dmso": 565.915,
                        "h2o": 566.166,
                        "methanol": 565.843,
                        "thf": 566.084,
                        "toluene": 566.506,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 564.429,
                        "acetone": 566.244,
                        "chcl3": 566.161,
                        "acetonitrile": 566.279,
                        "ch2cl2": 566.206,
                        "dmso": 566.349,
                        "h2o": 566.646,
                        "methanol": 566.247,
                        "thf": 566.233,
                        "toluene": 566.086,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 618.99,
                        "acetone": 624.483,
                        "chcl3": 623.499,
                        "acetonitrile": 624.639,
                        "ch2cl2": 624.087,
                        "dmso": 624.744,
                        "h2o": 625.072,
                        "methanol": 624.593,
                        "thf": 623.983,
                        "toluene": 622.448,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 583.324,
                        "acetone": 585.797,
                        "chcl3": 585.592,
                        "acetonitrile": 585.848,
                        "ch2cl2": 585.715,
                        "dmso": 585.925,
                        "h2o": 586.235,
                        "methanol": 585.813,
                        "thf": 585.725,
                        "toluene": 585.371,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 574.839,
                        "acetone": 574.09,
                        "chcl3": 574.267,
                        "acetonitrile": 574.11,
                        "ch2cl2": 574.167,
                        "dmso": 574.166,
                        "h2o": 574.435,
                        "methanol": 574.084,
                        "thf": 574.22,
                        "toluene": 574.478,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 569.911,
                        "acetone": 570.088,
                        "chcl3": 570.127,
                        "acetonitrile": 570.133,
                        "ch2cl2": 570.135,
                        "dmso": 570.198,
                        "h2o": 570.482,
                        "methanol": 570.103,
                        "thf": 570.164,
                        "toluene": 570.119,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 566.08,
                        "acetone": 564.411,
                        "chcl3": 564.793,
                        "acetonitrile": 564.406,
                        "ch2cl2": 564.583,
                        "dmso": 564.448,
                        "h2o": 564.684,
                        "methanol": 564.385,
                        "thf": 564.658,
                        "toluene": 565.213,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 564.63,
                        "acetone": 564.706,
                        "chcl3": 564.726,
                        "acetonitrile": 564.75,
                        "ch2cl2": 564.72,
                        "dmso": 564.813,
                        "h2o": 565.093,
                        "methanol": 564.721,
                        "thf": 564.752,
                        "toluene": 564.742,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 619.182,
                        "acetone": 623.189,
                        "chcl3": 622.29,
                        "acetonitrile": 623.352,
                        "ch2cl2": 622.833,
                        "dmso": 623.451,
                        "h2o": 623.764,
                        "methanol": 623.308,
                        "thf": 622.734,
                        "toluene": 621.304,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 583.522,
                        "acetone": 584.278,
                        "chcl3": 584.168,
                        "acetonitrile": 584.337,
                        "ch2cl2": 584.241,
                        "dmso": 584.407,
                        "h2o": 584.701,
                        "methanol": 584.305,
                        "thf": 584.256,
                        "toluene": 584.034,
                    }
                },
            },
        },
        "TMP": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 291.33,
                        "acetone": 276.264,
                        "chcl3": 277.254,
                        "acetonitrile": 275.207,
                        "ch2cl2": 276.171,
                        "dmso": 276.988,
                        "h2o": 262.671,
                        "methanol": 263.366,
                        "thf": 278.685,
                        "toluene": 283.761,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 277.761,
                        "acetone": 262.673,
                        "chcl3": 263.634,
                        "acetonitrile": 261.631,
                        "ch2cl2": 262.58,
                        "dmso": 263.406,
                        "h2o": 249.27,
                        "methanol": 249.931,
                        "thf": 265.061,
                        "toluene": 270.123,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 299.195,
                        "acetone": 286.35,
                        "chcl3": 287.213,
                        "acetonitrile": 285.469,
                        "ch2cl2": 286.302,
                        "dmso": 286.997,
                        "h2o": 274.843,
                        "methanol": 275.42,
                        "thf": 288.362,
                        "toluene": 292.724,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 277.52,
                        "acetone": 262.317,
                        "chcl3": 263.295,
                        "acetonitrile": 261.26,
                        "ch2cl2": 262.227,
                        "dmso": 263.036,
                        "h2o": 248.805,
                        "methanol": 249.485,
                        "thf": 264.716,
                        "toluene": 269.816,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 390.602,
                        "acetone": 379.7,
                        "chcl3": 380.279,
                        "acetonitrile": 378.978,
                        "ch2cl2": 379.593,
                        "dmso": 380.317,
                        "h2o": 368.831,
                        "methanol": 369.216,
                        "thf": 381.391,
                        "toluene": 384.986,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 297.198,
                        "acetone": 281.884,
                        "chcl3": 282.896,
                        "acetonitrile": 280.816,
                        "ch2cl2": 281.794,
                        "dmso": 282.606,
                        "h2o": 268.382,
                        "methanol": 269.076,
                        "thf": 284.334,
                        "toluene": 289.473,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 286.404,
                        "acetone": 270.748,
                        "chcl3": 271.725,
                        "acetonitrile": 269.462,
                        "ch2cl2": 270.524,
                        "dmso": 271.355,
                        "h2o": 256.342,
                        "methanol": 257.122,
                        "thf": 273.469,
                        "toluene": 278.676,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 272.706,
                        "acetone": 257.164,
                        "chcl3": 258.119,
                        "acetonitrile": 255.895,
                        "ch2cl2": 256.94,
                        "dmso": 257.797,
                        "h2o": 242.92,
                        "methanol": 243.667,
                        "thf": 259.855,
                        "toluene": 264.973,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 294.405,
                        "acetone": 281.158,
                        "chcl3": 282.018,
                        "acetonitrile": 280.073,
                        "ch2cl2": 280.993,
                        "dmso": 281.703,
                        "h2o": 269.086,
                        "methanol": 269.737,
                        "thf": 283.464,
                        "toluene": 287.882,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 272.595,
                        "acetone": 256.861,
                        "chcl3": 257.836,
                        "acetonitrile": 255.578,
                        "ch2cl2": 256.643,
                        "dmso": 257.483,
                        "h2o": 242.627,
                        "methanol": 243.389,
                        "thf": 259.577,
                        "toluene": 264.773,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 385.991,
                        "acetone": 374.828,
                        "chcl3": 375.394,
                        "acetonitrile": 373.92,
                        "ch2cl2": 374.61,
                        "dmso": 375.349,
                        "h2o": 363.431,
                        "methanol": 363.874,
                        "thf": 376.762,
                        "toluene": 380.401,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 292.227,
                        "acetone": 276.414,
                        "chcl3": 277.413,
                        "acetonitrile": 275.12,
                        "ch2cl2": 276.191,
                        "dmso": 277.05,
                        "h2o": 262.135,
                        "methanol": 262.912,
                        "thf": 279.163,
                        "toluene": 284.4,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 286.331,
                        "acetone": 271.022,
                        "chcl3": 271.947,
                        "acetonitrile": 269.751,
                        "ch2cl2": 270.768,
                        "dmso": 271.616,
                        "h2o": 256.882,
                        "methanol": 257.6,
                        "thf": 273.659,
                        "toluene": 278.687,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 272.619,
                        "acetone": 257.298,
                        "chcl3": 258.198,
                        "acetonitrile": 256.053,
                        "ch2cl2": 257.051,
                        "dmso": 257.926,
                        "h2o": 243.408,
                        "methanol": 244.095,
                        "thf": 259.935,
                        "toluene": 264.977,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 294.334,
                        "acetone": 281.319,
                        "chcl3": 282.131,
                        "acetonitrile": 280.265,
                        "ch2cl2": 281.144,
                        "dmso": 281.852,
                        "h2o": 269.472,
                        "methanol": 270.068,
                        "thf": 283.556,
                        "toluene": 287.875,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 272.586,
                        "acetone": 257.148,
                        "chcl3": 258.069,
                        "acetonitrile": 255.901,
                        "ch2cl2": 256.919,
                        "dmso": 257.755,
                        "h2o": 243.195,
                        "methanol": 243.894,
                        "thf": 259.785,
                        "toluene": 264.863,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 385.897,
                        "acetone": 374.881,
                        "chcl3": 375.407,
                        "acetonitrile": 373.999,
                        "ch2cl2": 374.652,
                        "dmso": 375.391,
                        "h2o": 363.697,
                        "methanol": 364.097,
                        "thf": 376.757,
                        "toluene": 380.319,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 292.105,
                        "acetone": 276.574,
                        "chcl3": 277.519,
                        "acetonitrile": 275.313,
                        "ch2cl2": 276.339,
                        "dmso": 277.197,
                        "h2o": 262.553,
                        "methanol": 263.276,
                        "thf": 279.247,
                        "toluene": 284.37,
                    }
                },
            },
        },
    }
    si_tm_shieldings = {
        "TMS": {
            "r2scan-3c": {
                "r2scan-3c": {
                    "def2-mTZVPP": {
                        "gas": 358.4765853,
                        "acetone": 357.8505313,
                        "chcl3": 357.8603583,
                        "acetonitrile": 357.8319377,
                        "ch2cl2": 357.8963768,
                        "dmso": 357.874638,
                        "h2o": 357.9281381,
                        "methanol": 357.9240386,
                        "thf": 357.8718386,
                        "toluene": 357.8388333,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 425.427401,
                        "acetone": 425.1035618,
                        "chcl3": 425.0727607,
                        "acetonitrile": 425.1110602,
                        "ch2cl2": 425.0604369,
                        "dmso": 425.1396548,
                        "h2o": 425.1994123,
                        "methanol": 425.1168623,
                        "thf": 425.1081863,
                        "toluene": 425.1184486,
                    }
                },
                "b97-3c": {
                    "def2-mTZVP": {
                        "gas": 352.2983555,
                        "acetone": 345.1994374,
                        "chcl3": 351.6272162,
                        "acetonitrile": 351.6144709,
                        "ch2cl2": 351.635666,
                        "dmso": 351.6499871,
                        "h2o": 351.7346671,
                        "methanol": 351.6729095,
                        "thf": 351.6321795,
                        "toluene": 351.648296,
                    }
                },
                "tpss": {
                    "def2-TZVP": {
                        "gas": 334.0278062,
                        "acetone": 333.4843744,
                        "chcl3": 333.5002814,
                        "acetonitrile": 333.4113907,
                        "ch2cl2": 333.5231986,
                        "dmso": 333.5126664,
                        "h2o": 333.555991,
                        "methanol": 333.5402484,
                        "thf": 333.5146132,
                        "toluene": 333.4768581,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 331.8236884,
                        "acetone": 331.3040747,
                        "chcl3": 331.3215705,
                        "acetonitrile": 331.2634015,
                        "ch2cl2": 331.3385284,
                        "dmso": 331.3228128,
                        "h2o": 331.3437693,
                        "methanol": 331.3449225,
                        "thf": 331.3305993,
                        "toluene": 331.3160006,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 340.923509,
                        "acetone": 340.3208483,
                        "chcl3": 340.3430966,
                        "acetonitrile": 340.3155874,
                        "ch2cl2": 340.3599861,
                        "dmso": 340.3553012,
                        "h2o": 340.4180652,
                        "methanol": 340.3808603,
                        "thf": 340.348664,
                        "toluene": 340.3406705,
                    }
                },
            },
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 334.2579542,
                        "acetone": 334.1639413,
                        "chcl3": 334.1459912,
                        "acetonitrile": 334.1644763,
                        "ch2cl2": 334.143167,
                        "dmso": 334.2355086,
                        "h2o": 334.1700712,
                        "methanol": 334.1638302,
                        "thf": 334.1765686,
                        "toluene": 334.1672644,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 332.1432161,
                        "acetone": 332.0806043,
                        "chcl3": 332.027555,
                        "acetonitrile": 332.070525,
                        "ch2cl2": 332.0181509,
                        "dmso": 332.1389588,
                        "h2o": 332.0768365,
                        "methanol": 332.082777,
                        "thf": 332.0989747,
                        "toluene": 332.0655251,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 425.4500968,
                        "acetone": 425.4194168,
                        "chcl3": 425.3783658,
                        "acetonitrile": 425.4187809,
                        "ch2cl2": 425.3492293,
                        "dmso": 425.4302912,
                        "h2o": 425.4004059,
                        "methanol": 425.3865089,
                        "thf": 425.4157351,
                        "toluene": 425.4555181,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 334.5698984,
                        "acetone": 334.0803779,
                        "chcl3": 334.1093328,
                        "acetonitrile": 334.0665281,
                        "ch2cl2": 334.1280337,
                        "dmso": 334.1272572,
                        "h2o": 334.0495564,
                        "methanol": 334.1137413,
                        "thf": 334.1251606,
                        "toluene": 334.1235476,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 332.3546979,
                        "acetone": 331.9058869,
                        "chcl3": 331.8955148,
                        "acetonitrile": 331.8800833,
                        "ch2cl2": 331.9140658,
                        "dmso": 331.948424,
                        "h2o": 331.8617288,
                        "methanol": 331.9375391,
                        "thf": 331.9562723,
                        "toluene": 331.9253075,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 426.0062656,
                        "acetone": 425.7811084,
                        "chcl3": 425.7602588,
                        "acetonitrile": 425.745999,
                        "ch2cl2": 425.7473718,
                        "dmso": 425.779427,
                        "h2o": 425.7365851,
                        "methanol": 425.7713265,
                        "thf": 425.7964293,
                        "toluene": 425.8200844,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 333.7779314,
                        "acetone": 333.3511708,
                        "chcl3": 333.3794838,
                        "acetonitrile": 333.3298692,
                        "ch2cl2": 333.3946486,
                        "dmso": 333.3881767,
                        "h2o": 333.3406562,
                        "methanol": 333.3784136,
                        "thf": 333.3860666,
                        "toluene": 333.3885135,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 331.5820841,
                        "acetone": 331.1904714,
                        "chcl3": 331.1839521,
                        "acetonitrile": 331.1565218,
                        "ch2cl2": 331.1982524,
                        "dmso": 331.2347884,
                        "h2o": 331.1670301,
                        "methanol": 331.2231923,
                        "thf": 331.2383692,
                        "toluene": 331.2108329,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 425.0726297,
                        "acetone": 424.9009564,
                        "chcl3": 424.8706079,
                        "acetonitrile": 424.8831877,
                        "ch2cl2": 424.8554965,
                        "dmso": 424.9143792,
                        "h2o": 424.8579037,
                        "methanol": 424.8851226,
                        "thf": 424.9146175,
                        "toluene": 424.9330242,
                    }
                },
            },
        }
    }
    si_orca_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 344.281,
                        "acetone": 344.239,
                        "chcl3": 344.311,
                        "acetonitrile": 344.198,
                        "ch2cl2": 344.231,
                        "dmso": 344.292,
                        "h2o": 344.228,
                        "methanol": 344.291,
                        "thf": 344.283,
                        "toluene": 344.452,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 332.181,
                        "acetone": 332.067,
                        "chcl3": 332.162,
                        "acetonitrile": 332.033,
                        "ch2cl2": 332.082,
                        "dmso": 332.122,
                        "h2o": 332.048,
                        "methanol": 332.122,
                        "thf": 332.134,
                        "toluene": 332.298,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 357.874,
                        "acetone": 357.762,
                        "chcl3": 357.864,
                        "acetonitrile": 357.726,
                        "ch2cl2": 357.783,
                        "dmso": 357.798,
                        "h2o": 357.715,
                        "methanol": 357.809,
                        "thf": 357.826,
                        "toluene": 358.001,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 335.739,
                        "acetone": 335.641,
                        "chcl3": 335.74,
                        "acetonitrile": 335.606,
                        "ch2cl2": 335.659,
                        "dmso": 335.687,
                        "h2o": 335.608,
                        "methanol": 335.692,
                        "thf": 335.707,
                        "toluene": 335.879,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 425.385,
                        "acetone": 425.52,
                        "chcl3": 425.527,
                        "acetonitrile": 425.511,
                        "ch2cl2": 425.508,
                        "dmso": 425.578,
                        "h2o": 425.566,
                        "methanol": 425.557,
                        "thf": 425.54,
                        "toluene": 425.556,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 341.186,
                        "acetone": 341.197,
                        "chcl3": 341.284,
                        "acetonitrile": 341.166,
                        "ch2cl2": 341.208,
                        "dmso": 341.263,
                        "h2o": 341.201,
                        "methanol": 341.253,
                        "thf": 341.263,
                        "toluene": 341.446,
                    }
                },
            },
            "b97-3c": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 344.503,
                        "acetone": 344.558,
                        "chcl3": 344.676,
                        "acetonitrile": 344.487,
                        "ch2cl2": 344.537,
                        "dmso": 344.67,
                        "h2o": 344.542,
                        "methanol": 344.662,
                        "thf": 344.637,
                        "toluene": 344.919,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 332.338,
                        "acetone": 332.293,
                        "chcl3": 332.442,
                        "acetonitrile": 332.236,
                        "ch2cl2": 332.31,
                        "dmso": 332.4,
                        "h2o": 332.288,
                        "methanol": 332.392,
                        "thf": 332.403,
                        "toluene": 332.676,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 357.729,
                        "acetone": 357.628,
                        "chcl3": 357.774,
                        "acetonitrile": 357.578,
                        "ch2cl2": 357.655,
                        "dmso": 357.692,
                        "h2o": 357.632,
                        "methanol": 357.703,
                        "thf": 357.725,
                        "toluene": 357.985,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 335.744,
                        "acetone": 335.688,
                        "chcl3": 335.837,
                        "acetonitrile": 335.633,
                        "ch2cl2": 335.71,
                        "dmso": 335.774,
                        "h2o": 335.704,
                        "methanol": 335.776,
                        "thf": 335.792,
                        "toluene": 336.064,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 425.911,
                        "acetone": 426.14,
                        "chcl3": 426.185,
                        "acetonitrile": 426.113,
                        "ch2cl2": 426.124,
                        "dmso": 426.254,
                        "h2o": 426.162,
                        "methanol": 426.22,
                        "thf": 426.196,
                        "toluene": 426.294,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 341.631,
                        "acetone": 341.666,
                        "chcl3": 341.811,
                        "acetonitrile": 341.61,
                        "ch2cl2": 341.676,
                        "dmso": 341.798,
                        "h2o": 341.602,
                        "methanol": 341.777,
                        "thf": 341.781,
                        "toluene": 342.086,
                    }
                },
            },
            "tpss": {
                "tpss": {
                    "def2-TZVP": {
                        "gas": 343.24,
                        "acetone": 343.388,
                        "chcl3": 343.506,
                        "acetonitrile": 343.343,
                        "ch2cl2": 343.385,
                        "dmso": 343.48,
                        "h2o": 343.378,
                        "methanol": 343.47,
                        "thf": 343.449,
                        "toluene": 343.647,
                    }
                },
                "pbe0": {
                    "def2-TZVP": {
                        "gas": 331.055,
                        "acetone": 331.217,
                        "chcl3": 331.313,
                        "acetonitrile": 331.175,
                        "ch2cl2": 331.224,
                        "dmso": 331.303,
                        "h2o": 331.205,
                        "methanol": 331.296,
                        "thf": 331.293,
                        "toluene": 331.461,
                    }
                },
                "dsd-blyp": {
                    "def2-TZVPP": {
                        "gas": 357.099,
                        "acetone": 357.125,
                        "chcl3": 357.231,
                        "acetonitrile": 357.081,
                        "ch2cl2": 357.141,
                        "dmso": 357.179,
                        "h2o": 357.075,
                        "methanol": 357.188,
                        "thf": 357.195,
                        "toluene": 357.379,
                    }
                },
                "wb97x": {
                    "def2-TZVP": {
                        "gas": 334.802,
                        "acetone": 334.886,
                        "chcl3": 334.987,
                        "acetonitrile": 334.842,
                        "ch2cl2": 334.897,
                        "dmso": 334.957,
                        "h2o": 334.855,
                        "methanol": 334.958,
                        "thf": 334.959,
                        "toluene": 335.134,
                    }
                },
                "pbeh-3c": {
                    "def2-mSVP": {
                        "gas": 424.346,
                        "acetone": 424.653,
                        "chcl3": 424.66,
                        "acetonitrile": 424.64,
                        "ch2cl2": 424.633,
                        "dmso": 424.74,
                        "h2o": 424.718,
                        "methanol": 424.709,
                        "thf": 424.681,
                        "toluene": 424.701,
                    }
                },
                "kt2": {
                    "def2-TZVP": {
                        "gas": 340.026,
                        "acetone": 340.228,
                        "chcl3": 340.311,
                        "acetonitrile": 340.189,
                        "ch2cl2": 340.226,
                        "dmso": 340.332,
                        "h2o": 340.207,
                        "methanol": 340.311,
                        "thf": 340.302,
                        "toluene": 340.453,
                    }
                },
            },
        }
    }

    def NMRRef_to_dict(self):
        """Convert NMRRef data to a dict object"""
        dict_ret = dict(
            h_tm_shieldings=self.h_tm_shieldings,
            c_tm_shieldings=self.c_tm_shieldings,
            f_tm_shieldings=self.f_tm_shieldings,
            si_tm_shieldings=self.si_tm_shieldings,
            p_tm_shieldings=self.p_tm_shieldings,
            h_orca_shieldings=self.h_orca_shieldings,
            c_orca_shieldings=self.c_orca_shieldings,
            f_orca_shieldings=self.f_orca_shieldings,
            si_orca_shieldings=self.si_orca_shieldings,
            p_orca_shieldings=self.p_orca_shieldings,
        )
        return dict_ret

    def dict_to_NMRRef(self, dictionary):
        """Convert dict object to NMRRef data """
        NmrRef_object = NmrRef()
        NmrRef_object.h_tm_shieldings = dictionary.get(
            "h_tm_shieldings", NmrRef_object.h_tm_shieldings
        )
        NmrRef_object.c_tm_shieldings = dictionary.get(
            "c_tm_shieldings", NmrRef_object.c_tm_shieldings
        )
        NmrRef_object.f_tm_shieldings = dictionary.get(
            "f_tm_shieldings", NmrRef_object.f_tm_shieldings
        )
        NmrRef_object.si_tm_shieldings = dictionary.get(
            "si_tm_shieldings", NmrRef_object.si_tm_shieldings
        )
        NmrRef_object.p_tm_shieldings = dictionary.get(
            "p_tm_shieldings", NmrRef_object.p_tm_shieldings
        )
        NmrRef_object.h_orca_shieldings = dictionary.get(
            "h_orca_shieldings", NmrRef_object.h_orca_shieldings
        )
        NmrRef_object.c_orca_shieldings = dictionary.get(
            "c_orca_shieldings", NmrRef_object.c_orca_shieldings
        )
        NmrRef_object.f_orca_shieldings = dictionary.get(
            "f_orca_shieldings", NmrRef_object.f_orca_shieldings
        )
        NmrRef_object.si_orca_shieldings = dictionary.get(
            "si_orca_shieldings", NmrRef_object.si_orca_shieldings
        )
        NmrRef_object.p_orca_shieldings = dictionary.get(
            "p_orca_shieldings", NmrRef_object.p_orca_shieldings
        )
        return NmrRef_object


# rotational entropy from symmetry
# https://cccbdb.nist.gov/thermo.asp
rot_sym_num = {
    "c1": 1,
    "ci": 1,
    "cs": 1,
    "c2": 2,
    "c3": 3,
    "c4": 4,
    "c5": 5,
    "c6": 6,
    "c7": 7,
    "c8": 8,
    "c9": 9,
    "c10": 10,
    "c11": 11,
    "s4": 2,
    "s6": 3,
    "s8": 4,
    "d2": 4,
    "d3": 6,
    "d4": 8,
    "d5": 10,
    "d6": 12,
    "d7": 14,
    "d8": 16,
    "d9": 18,
    "d10": 20,
    "t": 12,
    "th": 12,
    "td": 12,
    "o": 24,
    "oh": 24,
    "ih": 60,
}


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

# qm_prepinfo: grid and scfconv settings for ORCA and TM
qm_prepinfo = {
    "orca": {
        "low": ["grid4 nofinalgrid", "loosescf"],
        "low+": ["grid4 nofinalgrid", "scfconv6"],
        "high": ["grid4 nofinalgrid", "scfconv7"],
        "high+": ["grid5 nofinalgrid", "scfconv7"],
    },
    "tm": {
        "low": ["-grid", "m3", "-scfconv", "6"],
        "low+": ["-grid", "m4", "-scfconv", "6"],
        "high": ["-grid", "m4", "-scfconv", "7"],
        "high+": ["-grid", "m5", "-scfconv", "7"],
    },
}
