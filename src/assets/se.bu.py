from collections import OrderedDict

def main():
    # TODO - add default settings for uv/vis
    defaults_refine_ensemble_general = [
        # general settings
    ]
    defaults_refine_ensemble_part0 = [
        # part0
    ]
    defaults_refine_ensemble_part1 = [
        # part1
    ]
    # TODO - review crestcheck implementation
    defaults_refine_ensemble_part2 = [
        # part2
    ]
    defaults_refine_ensemble_part3 = [
        # part3
    ]
    defaults_nmrprop_part4 = [
        # part4
    ]
    defaults_optical_rotation_part5 = [
        # part5
    ]

    # TODO - defaults for uvvis (func: w-B97X-D3)

    internal_defaults = OrderedDict(
        defaults_refine_ensemble_general
        + defaults_refine_ensemble_part0
        + defaults_refine_ensemble_part1
        + defaults_refine_ensemble_part2
        + defaults_refine_ensemble_part3
        + defaults_nmrprop_part4
        + defaults_optical_rotation_part5
    )

    value_options = {
        "nconf": ["all", "number e.g. 10 up to all conformers"],
        "charge": ["number e.g. 0"],
        "unpaired": ["number e.g. 0"],
        "solvent": ["gas"],
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
        "part0_gfnv": "kok",
        "part1_gfnv": "kok",
        "part2_gfnv": "kok",
        "part3_gfnv": "kok",
        "temperature": ["temperature in K e.g. 298.15"],
        "multitemp": ["on", "off"],
        "trange": ["temperature range [start, end, step]"],
        "func0": "kok",
        "basis0": ["automatic"],
        "func": "kok",
        "basis": ["automatic"],
        "func3": "kok",
        "basis3": "kok",
        "part0_threshold": ["number e.g. 4.0"],
        "part1_threshold": ["number e.g. 5.0"],
        "opt_limit": ["number e.g. 4.0"],
        "part2_P_threshold": [
            "Boltzmann sum threshold in %. e.g. 95 (between 1 and 100)"
        ],
        "part3_threshold": [
            "Boltzmann sum threshold in %. e.g. 95 (between 1 and 100)"
        ],
        "sm2": "kok",
        "smgsolv3": "kok",
        "sm4_j": "kok",
        "sm4_s": "kok",
        "check": ["on", "off"],
        "crestcheck": ["on", "off"],
        "maxthreads": ["number of threads e.g. 2"],
        "omp": ["number cores per thread e.g. 4"],
        "smgsolv1": "kok",
        "smgsolv2": "kok",
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
        "prog4_j": ["tm", "orca", "prog"],
        "prog4_s": ["tm", "orca", "prog"],
        "func_j": "kok",
        "basis_j": "kok",
        "func_s": "kok",
        "basis_s": "kok",
        "h_ref": "kok",
        "c_ref": "kok",
        "f_ref": "kok",
        "si_ref": "kok",
        "p_ref": "kok",
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
        "cosmorsparam": ["automatic"],
        "balance": ["on", "off"],
        "progress": ["on", "off"],
    }

    for key in value_options.keys():
        if key not in internal_defaults.keys():
            print(key)


if __name__ == "__main__":
    main()
