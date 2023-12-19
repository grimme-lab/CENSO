import json
from types import MappingProxyType
from functools import reduce
from pprint import pprint


def main():
    settings_options = MappingProxyType(
        {
            int: MappingProxyType(
                {
                    "general": MappingProxyType(
                        {
                            "charge": {
                                "default": 0,
                                "range": (-10, 10),
                            },  # TODO - (re)move?
                            "unpaired": {
                                "default": 0,
                                "range": (0, 14),
                            },  # TODO - (re)move?
                            "maxprocs": {
                                "default": 1,
                                "range": (1, 128),
                            },  # number of processes
                            "omp": {
                                "default": 1,
                                "range": (1, 256),
                            },  # number of cores per process
                        }
                    ),
                    "prescreening": None,
                    "screening": None,
                    "optimization": MappingProxyType(
                        {
                            "optcycles": {
                                "default": 8,
                                "range": (1, 100),
                            },  # TODO - which value as min?
                            "radsize": {
                                "default": 10,
                                "range": (1, 100),
                            },  # TODO - same
                        }
                    ),
                    "refinement": None,
                    "nmr": None,
                    "optrot": None,
                    "uvvis": MappingProxyType(
                        {
                            "nroots": {
                                "default": 20,
                                "range": (1, 100),
                            },  # TODO - set dynamically
                        }
                    ),
                }
            ),
            float: MappingProxyType(
                {
                    "general": MappingProxyType(
                        {
                            "imagethr": {
                                "default": -100.0,
                                "range": (-300.0, 0.0),
                            },  # TODO - threshold for imaginary frequencies
                            "sthr": {
                                "default": 0.0,
                                "range": (0.0, 100.0),
                            },  # TODO - what is this?
                            "scale": {
                                "default": 1.0,
                                "range": (0.0, 1.0),
                            },  # TODO - what is this?
                            "temperature": {
                                "default": 298.15,
                                "range": (0.00001, 2000.0),
                            },  # TODO
                        }
                    ),
                    "prescreening": MappingProxyType(
                        {
                            "threshold": {
                                "default": 4.0,
                                "range": (1.0, 10.0),
                            },  # TODO - which value as min?
                        }
                    ),
                    "screening": MappingProxyType(
                        {
                            "threshold": {"default": 3.5, "range": (0.75, 7.5)},
                        }
                    ),
                    "optimization": MappingProxyType(
                        {
                            "threshold": {
                                "default": 2.5,
                                "range": (0.5, 5),
                            },  # TODO - rename?
                            "hlow": {"default": 0.01, "range": (0.01, 1.0)},  # TODO
                            "optimization_P_threshold": {
                                "default": 99.0,
                                "range": (1.0, 10.0),
                            },  # TODO
                            "spearmanthr": {"default": 0.0, "range": (0.0, 10.0)},
                        }
                    ),
                    "refinement": MappingProxyType(
                        {
                            "threshold": {
                                "default": 99.0,
                                "range": (1.0, 10.0),
                            },  # TODO
                        }
                    ),
                    "nmr": MappingProxyType(
                        {
                            "resonance_frequency": {
                                "default": 300.0,
                                "range": (150.0, 1000.0),
                            },  # TODO
                        }
                    ),
                    "optrot": None,
                    "uvvis": MappingProxyType(
                        {
                            "sigma": {"default": 0.1, "range": (0.1, 1.0)},
                        }
                    ),
                }
            ),
            str: MappingProxyType(
                {
                    "general": MappingProxyType(
                        {
                            # TODO - removed gas here, since there is already 'gas-phase' option (also doesn't make sense)
                            "solvent": {
                                "default": "h2o",
                                "options": tuple([k for k in solvents_db.keys()]),
                            },
                            # "prog_rrho": {"default": "xtb", "options": ("xtb")}, # TODO - keep this here?
                            "sm_rrho": {
                                "default": "alpb",
                                "options": ("alpb", "gbsa"),
                            },  # TODO - same
                            "cosmorsparam": {
                                "default": "automatic",
                                "options": tuple([k for k in cosmors_param.keys()]),
                            },
                        }
                    ),
                    "prescreening": MappingProxyType(
                        {
                            "func": {
                                "default": "pbe-d4",
                                "options": tuple(
                                    dfa_settings.find_func("prescreening")
                                ),
                            },
                            "basis": {
                                "default": "def2-SV(P)",
                                "options": ("automatic",)
                                + tuple(dfa_settings.composite_bs)
                                + ("def2-SV(P)", "def2-TZVP"),
                            },
                            "prog": {"default": "orca", "options": PROGS},
                            "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                            "grid": {"default": "low", "options": GRIDOPTIONS},
                        }
                    ),
                    "screening": MappingProxyType(
                        {
                            "func": {
                                "default": "r2scan-3c",
                                "options": tuple(dfa_settings.find_func("func1")),
                            },
                            "basis": {
                                "default": "automatic",
                                "options": ("automatic",) + basis_sets,
                            },
                            "prog": {"default": "orca", "options": PROGS},
                            "smgsolv": {"default": "smd_gsolv", "options": gsolv_mods},
                            "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                            "grid": {"default": "low+", "options": GRIDOPTIONS},
                        }
                    ),
                    "ensembleopt": MappingProxyType(
                        {
                            "func": {
                                "default": "r2scan-3c",
                                "options": tuple(dfa_settings.find_func("func2")),
                            },
                            "basis": {
                                "default": "automatic",
                                "options": ("automatic",) + basis_sets,
                            },
                            "prog": {"default": "orca", "options": PROGS},
                            "prog2opt": {
                                "default": "prog",
                                "options": PROGS,
                            },  # TODO - ??? prog2 ??? # FIXME
                            "sm": {"default": "smd", "options": solv_mods},  # FIXME
                            "smgsolv": {"default": "smd_gsolv", "options": gsolv_mods},
                            "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                            "optlevel2": {
                                "default": "automatic",
                                "options": (
                                    "crude",
                                    "sloppy",
                                    "loose",
                                    "lax",
                                    "normal",
                                    "tight",
                                    "vtight",
                                    "extreme",
                                    "automatic",
                                ),
                            },  # TODO - what does this mean?
                            "grid": {"default": "high", "options": GRIDOPTIONS},
                        }
                    ),
                    "refinement": MappingProxyType(
                        {
                            "prog": {"default": "orca", "options": PROGS},
                            "func": {
                                "default": "wb97x-v",
                                "options": tuple(dfa_settings.find_func("func3")),
                            },
                            "basis": {"default": "def2-TZVPP", "options": basis_sets},
                            "smgsolv": {"default": "smd_gsolv", "options": gsolv_mods},
                            "gfnv": {"default": "gfn2", "options": GFNOPTIONS},
                            "grid": {"default": "high+", "options": GRIDOPTIONS},
                        }
                    ),
                    "nmr": MappingProxyType(
                        {
                            "prog4_j": {"default": "tm", "options": PROGS},
                            "func_j": {
                                "default": "pbe0-d4",
                                "options": tuple(dfa_settings.find_func("func_j")),
                            },
                            "basis_j": {"default": "def2-TZVP", "options": basis_sets},
                            "sm4_j": {
                                "default": "default",
                                "options": solv_mods,
                            },  # FIXME
                            "prog4_s": {"default": "tm", "options": PROGS},
                            "func_s": {
                                "default": "pbe0-d4",
                                "options": tuple(dfa_settings.find_func("func_s")),
                            },
                            "basis_s": {"default": "def2-TZVP", "options": basis_sets},
                            "sm4_s": {
                                "default": "default",
                                "options": solv_mods,
                            },  # FIXME
                            "h_ref": {"default": "TMS", "options": ("TMS",)},
                            "c_ref": {"default": "TMS", "options": ("TMS",)},
                            "f_ref": {"default": "CFCl3", "options": ("CFCl3",)},
                            "si_ref": {"default": "TMS", "options": ("TMS",)},
                            "p_ref": {"default": "TMP", "options": ("TMP", "PH3")},
                        }
                    ),
                    "optrot": MappingProxyType(
                        {
                            "func": {
                                "default": "pbe-d4",
                                "options": tuple(dfa_settings.find_func("func_or")),
                            },
                            "func_or_scf": {
                                "default": "r2scan-3c",
                                "options": tuple(dfa_settings.find_func("func_or_scf")),
                            },
                            "basis": {"default": "def2-SVPD", "options": basis_sets},
                            "prog": {"default": "orca", "options": ("orca",)},
                        }
                    ),
                    "uvvis": None,
                }
            ),
            bool: MappingProxyType(
                {
                    "general": MappingProxyType(
                        {
                            "multitemp": {"default": True},
                            "evaluate_rrho": {"default": True},
                            "consider_sym": {"default": True},
                            "bhess": {"default": True},
                            "rmsdbias": {"default": False},
                            "progress": {"default": False},
                            "check": {"default": True},
                            "balance": {"default": True},
                            "vapor_pressure": {"default": False},
                            "nmrmode": {"default": False},
                            "gas-phase": {"default": False},
                        }
                    ),
                    "prescreening": MappingProxyType(
                        {
                            "run": {"default": True},
                            "gcp": {"default": True},
                        }
                    ),
                    "screening": MappingProxyType(
                        {
                            "run": {"default": True},
                            "gcp": {"default": True},
                        }
                    ),
                    "ensembleopt": MappingProxyType(
                        {
                            "run": {"default": True},
                            "gcp": {"default": True},
                            # "ancopt": {"default": True},
                            "opt_spearman": {"default": True},
                            "crestcheck": {"default": False},
                        }
                    ),
                    "refinement": MappingProxyType(
                        {
                            "run": {"default": False},
                            "gcp": {"default": True},
                        }
                    ),
                    "nmr": MappingProxyType(
                        {
                            "run": {"default": False},
                            "couplings": {"default": True},
                            "shieldings": {"default": True},
                            "h_active": {"default": True},
                            "c_active": {"default": True},
                            "f_active": {"default": False},
                            "si_active": {"default": False},
                            "p_active": {"default": False},
                        }
                    ),
                    "optrot": MappingProxyType(
                        {
                            "run": {"default": False},
                        }
                    ),
                    "uvvis": MappingProxyType(
                        {
                            "run": {"default": False},
                        }
                    ),
                }
            ),
            list: MappingProxyType(
                {
                    "general": MappingProxyType(
                        {
                            "trange": {"default": [273.15, 378.15, 5]},
                        }
                    ),
                    "prescreening": None,
                    "screening": None,
                    "ensembleopt": None,
                    "refinement": None,
                    "nmr": None,
                    "optrot": MappingProxyType(
                        {
                            "freq_or": {"default": [598.0]},
                        }
                    ),
                    "uvvis": None,
                }
            ),
        }
    )

    def replace_mappingproxytype(data):
        if isinstance(data, MappingProxyType):
            return dict((k, replace_mappingproxytype(v)) for k, v in data.items())
        elif isinstance(data, dict):
            return dict((k, replace_mappingproxytype(v)) for k, v in data.items())
        elif isinstance(data, list):
            return [replace_mappingproxytype(item) for item in data]
        else:
            return data

    settings_options = replace_mappingproxytype(settings_options)

    cleaned_options = {}

    for type_t, parts in settings_options.items():
        cleaned_parts = {}
        for part, settings in parts.items():
            if settings is not None:
                cleaned_parts[part] = settings

        cleaned_options[type_t] = cleaned_parts

    merged_dict = {}

    for type_key, type_value in cleaned_options.items():
        for section_key, section_value in type_value.items():
            merged_dict.setdefault(section_key, {}).update(section_value)

    import json

    print(json.dumps(merged_dict, indent=4))


if __name__ == "__main__":
    main()
