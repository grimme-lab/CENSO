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


descr = """
     __________________________________________________
    |                                                  |
    |                    UVVISPLOT                     |
    |        Plotting of ensemble UV/Vis spectra       |
    |             University of Bonn, MCTC             |
    |                   March 2024                     |
    |                     v 1.0.0                      |
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
        help="Set title of entire plot. If no title is required use " "'<--title ''>'.",
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
    data = read_data(args.inp)

    # Plot data
    figure = plot(data, args)

    # Save plot
    save_plot(figure, args.out)


if __name__ == "__main__":
    main()
