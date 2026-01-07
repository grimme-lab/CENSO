#!/usr/bin/env python3

#
# Copyright (C) 2024 Leopold M. Seidler
#
# UVVISPLOT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UVVISPLOT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UVVISPLOT. If not, see <https://www.gnu.org/licenses/>.

"""
Created on Mar 17, 2024
last updated on 17-March-2024
@author: lmseidler
"""

import matplotlib.pyplot as plt
import os
import argparse
import json
import numpy as np
import pandas as pd

PLANCK = 6.62607015e-34
C = 2.998e8
COULOMB = 1.602e-19
# Physical constants for Boltzmann weighting
AU2J = 4.3597482e-18  # Hartree to Joules
KB = 1.3806485279e-23  # Boltzmann constant in J/K


descr = """
     __________________________________________________
    |                                                  |
    |                    UVVISPLOT                     |
    |        Plotting of ensemble UV/Vis spectra       |
    |             University of Bonn, MCTC             |
    |                   January 2026                   |
    |                     v 1.1.0                      |
    |                  L. M. Seidler                   |
    |__________________________________________________|
    """


def get_args():
    """
    Parse command-line arguments.

    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
    )  # argparse.RawDescriptionHelpFormatter) #,
    parser.add_argument(
        "-mode",
        dest="mode",
        action="store",
        required=False,
        default="wavenumber",
        type=str,
        choices=["wavenumber", "energy", "wavelength"],
        help="Set the unit of the x-axis. Can be wavenumber [cm-1], energy [eV] or wavelength [nm].",
    )
    parser.add_argument(
        "-start",
        dest="start",
        action="store",
        required=False,
        type=float,
        help="Start plotting from '<start>'. Default values: 300nm/1.8eV/14000cm-1.",
    )
    parser.add_argument(
        "-end",
        dest="end",
        action="store",
        required=False,
        type=float,
        help="End plotting at '<end>'. '<end>' must be larger than '<start>'. Default values: 700nm/4.2eV/33000cm-1.",
    )
    parser.add_argument(
        "-title",
        "--title",
        dest="title",
        action="store",
        required=False,
        default="UVVis-PLOT",
        type=str,
        help="Set title of entire plot. If no title is required use '<--title ''>'.",
    )
    parser.add_argument(
        "-lw",
        "--linewidth",
        dest="lw",
        action="store",
        required=False,
        default=1.6131e3,
        type=float,
        help="Set linewidth in cm-1.",
    )
    parser.add_argument(
        "-i",
        "--inp",
        dest="inp",
        action="store",
        required=True,
        help="Provide input file.",
    )
    parser.add_argument(
        "-fontsize",
        "--fontsize",
        dest="fontsize",
        action="store",
        required=False,
        default=15,
        type=float,
        help="Set fontsize for entire plot.",
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="out",
        action="store",
        required=False,
        default="nmrplot",
        help="Provide name of the output file (including ending).",
    )
    parser.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        action="store",
        required=False,
        default=298.15,
        type=float,
        help="Temperature in Kelvin for Boltzmann weighting. Default: 298.15 K.",
    )
    args = parser.parse_args()
    return args


def read_data(inp):
    """
    Read data from input file.

    :param inp: Input file path.
    :return: Loaded data.
    """
    cwd = os.getcwd()
    with open(os.path.join(cwd, inp)) as f:
        data = json.load(f)

    return data


def process_data(json_data, temperature):
    """
    Process JSON data from CENSO uvvis output into flat format for plotting.

    Converts the hierarchical DataDump format (with conformer data and excitations)
    into a flat list of [wavelength, weighted_osc_str, conf_name] entries.
    Applies Boltzmann weighting based on free energies.

    :param json_data: JSON data loaded from file (either DataDump format or legacy flat list).
    :param temperature: Temperature in Kelvin for Boltzmann weighting.
    :return: List of [wavelength, weighted_intensity, conf_name] entries.
    """
    # Check if this is the new DataDump format
    if not isinstance(json_data, dict) or "data" not in json_data:
        # Assume legacy format (already a flat list)
        print("Legacy data format detected. Using data as-is.")
        return json_data

    print(
        f"Processing CENSO DataDump format with Boltzmann weighting at {temperature} K."
    )

    conformer_data = json_data["data"]

    # Step 1: Calculate Boltzmann populations
    # Calculate total free energy (gtot) for each conformer
    gtot_values = {}
    for conf_name, conf_info in conformer_data.items():
        gtot = conf_info["energy"] + conf_info["gsolv"] + conf_info["grrho"]
        gtot_values[conf_name] = gtot

    # Find minimum gtot
    min_gtot = min(gtot_values.values())

    # Calculate Boltzmann factors
    from math import exp

    boltzmann_factors = {}
    for conf_name, gtot in gtot_values.items():
        delta_g = gtot - min_gtot
        boltzmann_factors[conf_name] = exp(-(delta_g * AU2J) / (KB * temperature))

    # Normalize to get populations
    partition_sum = sum(boltzmann_factors.values())
    populations = {
        conf_name: factor / partition_sum
        for conf_name, factor in boltzmann_factors.items()
    }

    # Print population summary
    print("\nBoltzmann populations:")
    for conf_name in sorted(populations.keys()):
        print(f"  {conf_name}: {populations[conf_name]:.4f}")
    print()

    # Step 2: Flatten excitations with Boltzmann weighting
    flat_data = []
    for conf_name, conf_info in conformer_data.items():
        population = populations[conf_name]
        for excitation in conf_info["excitations"]:
            wavelength = excitation["wavelength"]
            osc_str = excitation["osc_str"]
            weighted_osc_str = osc_str * population
            flat_data.append([wavelength, weighted_osc_str, conf_name])

    print(
        f"Processed {len(flat_data)} excitations from {len(conformer_data)} conformers."
    )

    return flat_data


def plot(data, args):
    """
    Plot the UV-Vis spectrum.

    :param data: The data to plot.
    :param args: Parsed arguments.
    :return: The figure.
    """
    mode = args.mode

    # Select start value
    if args.start is not None:
        start = args.start
    else:
        defaults = {"wavelength": 300, "wavenumber": 14000, "energy": 1.8}
        start = defaults[mode]

    # Select end value
    if args.end is not None:
        end = args.end
    else:
        defaults = {"wavelength": 700, "wavenumber": 33000, "energy": 4.2}
        end = defaults[mode]

    assert end > start
    xrange = np.linspace(start, end, 10000)

    # Dump single contributions to csv file
    confs = {d[2] for d in data}
    exc_number = {conf: 0 for conf in confs}
    contributions_dict = {}

    for exc in data:
        yrange = gaussian_signal(xrange, exc[0], exc[1], args.lw, mode=mode)
        contributions_dict[f"{exc[2]}_S{exc_number[exc[2]]}"] = yrange

        exc_number[exc[2]] += 1

    cwd = os.getcwd()
    contributions_df = pd.DataFrame.from_dict(contributions_dict)
    contributions_df.to_csv(os.path.join(cwd, "contributions.csv"))
    print("All contributions written to contributions.csv.")

    # Plot the whole spectrum
    fig, ax = plt.subplots()
    yrange = contributions_df.sum(axis=1)
    ax.plot(xrange, yrange)
    ax.set_title(args.title)
    labels = {
        "wavelength": r"$\mathrm{nm}$",
        "wavenumber": r"$\mathrm{cm-1}$",
        "energy": r"$\mathrm{eV}$",
    }
    ax.set_xlabel(f"{args.mode} [{labels[args.mode]}]")
    ax.set_ylabel(r"$\epsilon$ [a. u.]")

    return fig


def gaussian_signal(xrange, center_wl, eps_max, lw, mode="wavelength"):
    """
    Generate gaussian signal.

    :param xrange: Range of x values.
    :param center_wl: Center wavelength.
    :param eps_max: Maximum epsilon.
    :param lw: Linewidth.
    :param mode: Mode for units.
    :return: The signal array.
    """
    # <=> 1/λ = E / (h c)
    # 1 nm = 1e-7 cm
    # 1 cm-1 = 1e7 nm-1
    if mode == "wavelength":
        return eps_max * np.exp(-(((1 / xrange - 1 / center_wl) / (lw * 1e7)) ** 2))
    elif mode == "wavenumber":
        return eps_max * np.exp(-(((xrange - 1 / center_wl * 1e7) / lw) ** 2))
    elif mode == "energy":
        return eps_max * np.exp(
            -(
                ((xrange * COULOMB / (PLANCK * C) - 1 / center_wl * 1e9) / (lw * 1e2))
                ** 2
            )
        )


def save_plot(fig, out):
    """
    Save the plot to file.

    :param fig: The figure.
    :param out: Output file path.
    """
    cwd = os.getcwd()
    fig.savefig(os.path.join(cwd, out), format="pdf")


def main():
    """
    Main execution function.
    """

    # Parse cml args
    args = get_args()

    # Read data
    raw_data = read_data(args.inp)

    # Process data (handles both new DataDump format and legacy flat format)
    processed_data = process_data(raw_data, args.temperature)

    # Plot data
    figure = plot(processed_data, args)

    # Save plot
    save_plot(figure, args.out)


if __name__ == "__main__":
    main()
