#!/usr/bin/env python3

#
# Copyright (C) 2019 Fabian Bohle
#
# NMRPLOT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# NMRPLOT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with NMRPLOT. If not, see <https://www.gnu.org/licenses/>.

"""
Created on Jan 11, 2019
last updated on 16-October-2025
@author: bohle, lmseidler
"""

##try and except handling of imports
try:
    import argparse
except ImportError:
    raise ImportError("    Error while importing argparse module!")
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tck
except ImportError:
    raise ImportError("    Error while importing from matplotlib!")
try:
    from matplotlib import gridspec
except ImportError:
    raise ImportError("    Error while importing matplotlib.gridspec!")
try:
    import numpy as np
except ImportError:
    raise ImportError("    Error while importing numpy!")
try:
    from typing import Any
except ImportError:
    raise ImportError("    Error while importing typing!")
try:
    # from sys import version_info  # unused
    from sys import exit
    from sys import argv as sysargv
except ImportError:
    raise ImportError("    Error while importing from sys!")


descr = """
     __________________________________________________
    |                                                  |
    |                    NMRPLOT                       |
    |          Plotting of NMR spectral data           |
    |             University of Bonn, MCTC             |
    |                 January 2019                     |
    |                     v 1.1.1                       |
    |                   F. Bohle                       |
    |__________________________________________________|
    """

useit = """\
    End     Endremove    Startremove                 Start
    +               +    +                               +
    +---------------+----+-------------------------------+
    lower field                               higher field
                        delta /ppm
    """


def checkval(value):
    """
    Check if value is between 0.0 and 1.0.

    :param value: The value to check.
    :return: The float value if valid.
    :raises argparse.ArgumentTypeError: If value is out of range.
    """
    x = float(value)
    if x < 0 or x > 1.0:
        raise argparse.ArgumentTypeError("{!r} not in range [0.0, 1.0]".format(x))
    return x


def cml(descr):
    """
    Get args object from commandline interface.

    Needs argparse module.

    :param descr: Description for the parser.
    :return: Tuple of parsed args and defaults.
    """
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
    )  # argparse.RawDescriptionHelpFormatter) #,
    parser.add_argument(
        "-start",
        "--startppm",
        dest="start",
        action="store",
        required=False,
        default=0,
        type=float,
        help="Start plotting from '<start>' ppm.",
    )
    parser.add_argument(
        "-end",
        "--endppm",
        dest="end",
        action="store",
        required=False,
        default=11,
        type=float,
        help="End plotting at '<end>' ppm. Value of end has to be larger than "
        "value of start.",
    )
    parser.add_argument(
        "-startremove",
        "--startremove",
        dest="startremove",
        action="store",
        required=False,
        default=None,
        type=float,
        help="Start cutting from spectrum at '<startremove>' ppm.",
    )
    parser.add_argument(
        "-endremove",
        "--endremove",
        dest="endremove",
        action="store",
        required=False,
        default=None,
        type=float,
        help="End cutting from spectrum at '<endremove>' ppm. "
        "Value of endremove has to be larger than value of startremove.",
    )
    parser.add_argument(
        "-title",
        "--title",
        dest="title",
        action="store",
        required=False,
        default="NMR-PLOT",
        type=str,
        help="Set title of entire plot. If no title is required use " "'<--title ''>'.",
    )
    parser.add_argument(
        "-lw",
        "--linewidth",
        dest="linewidth",
        action="store",
        required=False,
        default=0.8,
        type=float,
        help="Set linewidth.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=True,
        nargs="+",
        help="Provide input_file(s) [max 3 files] -i input1(theory1) input2(theory2)"
        " input3(experiment/predicition); inputfiles format is two columns: "
        "column1 ppm , column2 intensity",
    )
    parser.add_argument(
        "-l",
        "--label",
        dest="label",
        action="store",
        required=False,
        default=[],
        nargs="+",
        help="Provide labels for all files provided -l label1 label2 label3"
        "... , if no labels are provided, filename is used as label",
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
        "-keybox",
        "--keybox",
        dest="keybox",
        action="store_true",
        required=False,
        default=False,
        help="Set Frame around key.",
    )
    parser.add_argument(
        "-ontop",
        "--ontop",
        dest="ontop",
        action="store_true",
        required=False,
        default=False,
        help="Plot all spectra ontop of each other.",
    )
    parser.add_argument(
        "-stacked",
        "--stacked",
        dest="stacked",
        action="store_true",
        required=False,
        default=False,
        help="Plot all spectra stacked over each other.",
    )
    parser.add_argument(
        "-orientation",
        "--orientation",
        dest="orientation",
        action="store",
        required=False,
        nargs="+",
        type=int,
        default=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        help="Up (1) or down (-1).",
    )
    parser.add_argument(
        "-c",
        "--colors",
        dest="colors",
        action="store",
        nargs="+",
        choices=["gray", "blue", "cyan", "red", "green", "magenta", "yellow", "black"],
        required=False,
        metavar="",
        default=["blue", "black", "red", "magenta", "green"],
        help="Select colors. Possible are: {}".format(
            ["gray", "blue", "cyan", "red", "green", "magenta", "yellow", "black"]
        ),
    )
    args_defaults = {"colors": parser.get_default("colors")}
    parser.add_argument(
        "-cut",
        "--cut",
        dest="cut",
        action="store",
        nargs="+",
        required=False,
        type=checkval,  ### own function
        default=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        help="Cut intensity. Accepts values from 0.0 (flat line) to 1.0 (full intensity).",
    )
    args_defaults["cut"] = parser.get_default("cut")
    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="nmrplot",
        help="Provide name of the output file without fileending.",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        required=False,
        default=False,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-s",
        "--shift",
        dest="shift",
        action="store",
        nargs="+",
        type=float,
        required=False,
        default=[],
        help="Shift ppm of each inputfile separately using: --shift float float "
        "float, e.g. --shift 10.0 0.0 -5.0, each file needs its own value",
    )
    args = parser.parse_args()
    return args, args_defaults


def readinput(filename, ppm, intensit, number):
    """
    Read input from given filename into ppm and intensit lists.

    :param filename: Path to the input file.
    :param ppm: List to append ppm values.
    :param intensit: List to append intensity values.
    :param number: Number (unused).
    :return: Updated ppm and intensit lists.
    """
    with open(filename) as inp:
        data = inp.readlines()
    x = []
    y = []
    for line in data:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))
    ppm.append(x)
    intensit.append(y)
    return ppm, intensit


def axdefaultsettings(ax, args):
    """
    Set default settings for the axis.

    :param ax: The matplotlib axis.
    :param args: Parsed arguments.
    """
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(
        axis="x",
        which="both",
        bottom=False,
        top=False,
        labelbottom=False,
        labelsize=args.fontsize,
    )
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    ax.get_yaxis().set_visible(False)
    ###(matplotlib 2.x space between spine and axis Set padding of Y data limits prior to autoscaling.
    ax.set_ymargin(0.002)
    ax.autoscale(enable=True, axis="y", tight=True)
    return


def equal_ticks(axes1, axes2, args):
    """
    Equalize ticks between two axes.

    :param axes1: First axis.
    :param axes2: Second axis.
    :param args: Parsed arguments.
    :return: Modified axes1 and axes2.
    """
    xticksax1 = abs(axes1.get_xticks()[0] - axes1.get_xticks()[1])
    xticksax2 = abs(axes2.get_xticks()[0] - axes2.get_xticks()[1])
    if round(xticksax1, 2) >= round(xticksax2, 2):
        nex = args.startremove + xticksax1
        while True:
            if nex > args.endremove:
                break
            nex += xticksax1
        t1 = np.arange(nex, args.end, step=xticksax1)
        if (abs(t1[0]) - abs(args.endremove)) < round(xticksax1 / 2, 2):
            t1 = t1[1:]
        t2 = np.arange(args.start, args.startremove, step=xticksax1)
        axes2.set_xticks(t2)
        axes1.set_xticks(t1)
    else:
        nex = args.startremove + xticksax2
        while True:
            if nex > args.endremove:
                break
            nex += xticksax2
        t1 = np.arange(nex, args.end, step=xticksax2)
        if (abs(t1[0]) - abs(args.endremove)) < round(xticksax2 / 2, 2):
            t1 = t1[1:]
        t2 = np.arange(args.start, args.startremove, step=xticksax2)
        axes2.set_xticks(t2)
        axes1.set_xticks(t1)

    return axes1, axes2


def debug():
    """
    Print debug information and exit.
    """
    print("    {}".format(plt.__file__))
    print("    {}".format(gridspec.__file__))
    print("    {}".format(np.__file__))
    exit(0)


def main():
    """
    Main execution function.

    :return: None
    """
    args, args_defaults = cml(descr)
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.debug:
        debug()

    ### errors from start end startremove endremove
    errors = []
    if not args.start < args.end:
        errors.append("    Error! Start has to be smaller than end!" " Going to exit!")
    if args.startremove is not None or args.endremove is not None:
        if args.startremove is None or args.endremove is None:
            errors.append(
                "    Error! Both startremove and endremove" " have to be set!"
            )
        else:
            if args.startremove is not None or args.endremove is not None:
                if not (args.startremove < args.endremove):
                    errors.append(
                        "    Error! Startremove has to be smaller than endremove!"
                        " Going to exit!"
                    )
                if args.start >= args.startremove:
                    errors.append(
                        "    Error! Startremove has to be larger than start! "
                        "Going to exit!"
                    )
                if args.end <= args.endremove:
                    errors.append(
                        "    Error! Endremove has to be smaller than end! Going "
                        "to exit!"
                    )
    if len(errors) >= 1:
        for item in errors:
            print(item)
        print(useit)
        exit(1)

    if len(args.label) < len(args.file):  ### labels from filenames
        print(
            "    Labels not provided or not provided for every file, "
            "using filenames as labels."
        )
        args.label = []
        for file in args.file:
            if file.split(".")[0] != "":
                args.label.append(file.split(".")[0])
            else:
                args.label.append(file.split(".")[1])

    if len(args.colors) < len(args.file):  ### colors
        print(
            "    Provided less colors than files, choosing "
            "defaults: {}.".format(args_defaults["colors"])
        )
        args.colors = args_defaults["colors"]
    if len(args.file) == 1:
        print("    Plotting 1 data file.")
    else:
        print("    Plotting {} data files.".format(len(args.file)))
    ### Get data from data files
    ppm: list[list[float]] = []
    intensit: list[Any] = []
    i = 0
    for file in args.file:
        try:
            ppm, intensit = readinput(file, ppm, intensit, i)
        except OSError:
            print(
                "    File: {} does not exist! Or Error while reading file! "
                "Terminating now!".format(file)
            )
            exit(1)
        i += 1
    ### ppm ranges of datafiles:
    for i in range(len(args.file)):
        print(
            "    ppm range goes from {:6.2f} to {:6.2f} in file: {}.".format(
                min(ppm[i]), max(ppm[i]), args.file[i]
            )
        )
    ### cut in case of too large intensities
    if args.cut:
        if len(args.cut) < len(args.file):
            print(
                "    Provided less cut values than files, therefore chosing defaults."
            )
            args.cut = args_defaults["cut"]
        for i in range(len(args.file)):
            newmax = np.amax(np.array(intensit[i])) * args.cut[i]
            # newmax= max(intensit[i]) * args.cut[i]
            intensit[i] = np.clip(np.array(intensit[i]), -10, newmax)  # .tolist()
            # intensit[i] = [newmax if x >= newmax else x for x in intensit[i]]
            # intensit[i] = [newmax if x <= newmax else x for x in intensit[i]]
    ### end cut

    ### shift manually:
    if len(args.shift) < len(args.file):
        print(
            "    Provided less or no shift values than files, therfore "
            "nothing is shifted!"
        )
    else:
        for i in range(len(args.file)):
            ppm[i] = [item + args.shift[i] for item in ppm[i]]
    figure = plt.figure(figsize=(11.69, 8.27))  ### A4 in inches
    figure.suptitle(args.title, fontsize=args.fontsize, y=0.93)
    # removed
    if args.startremove and args.endremove is not None:
        print(
            "    Plotting spectrum in range of {} to {} ppm,".format(
                args.start, args.end
            )
        )
        print(
            "    with data section between {} and {} ppm removed.".format(
                args.startremove, args.endremove
            )
        )
        # ratio:
        a = float(args.end) - float(args.endremove)
        b = float(args.startremove) - float(args.start)
        # ontop and removed
        if args.ontop:
            if len(args.file) == 1:
                print(
                    "    Plotting ontop with only one file is useless! Going to exit!"
                )
                exit(1)
            print("    Plotting ontop!")
            gs = gridspec.GridSpec(1, 2, width_ratios=[a / b, 1])
            ax1 = plt.subplot(gs[0])  # low field (high ppm)
            ax2 = plt.subplot(gs[1])  # high field (low ppm)
            maxint = []
            for i in range(len(args.file)):
                maxint.append(max([yval for yval in intensit[i]]))
                maxint[i] = maxint[i] + 0.000001 if maxint[i] == 0 else maxint[i]
            for i in range(len(args.file)):
                ax1.plot(
                    ppm[i],
                    [yval / maxint[i] for yval in intensit[i]],
                    label=args.label[i],
                    color=args.colors[i],
                    linewidth=args.linewidth,
                )
                axdefaultsettings(ax1, args)
                ax1.set_xlim(
                    left=args.end, right=args.endremove
                )  ### sets x limits must be here

                ax2.plot(
                    ppm[i],
                    [yval / maxint[i] for yval in intensit[i]],
                    label=args.label[i],
                    color=args.colors[i],
                    linewidth=args.linewidth,
                )
                axdefaultsettings(ax2, args)
                ax2.set_xlim(
                    left=args.startremove, right=args.start
                )  ### sets x limits must be here

            ax1.tick_params(axis="x", bottom=True, labelbottom=True)
            ax1.xaxis.set_tick_params(which="minor", bottom=True)
            ax1.xaxis.set_minor_locator(tck.AutoMinorLocator())
            ax1.spines["bottom"].set_visible(True)
            ax1.spines["bottom"].set_position(
                ("outward", 1)
            )  # set spine (in picture the x axis down by x points)
            ax1.get_yaxis().set_visible(False)
            ax2.tick_params(axis="x", bottom=True, labelbottom=True)
            ax2.xaxis.set_tick_params(which="minor", bottom=True)
            ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
            ax2.spines["bottom"].set_visible(True)
            ax2.spines["bottom"].set_position(
                ("outward", 1)
            )  # set spine (in picture the x axis down by x points)
            ax2.legend(loc="upper right", frameon=args.keybox, fontsize=args.fontsize)
            ### make ticks equal
            ax1, ax2 = equal_ticks(ax1, ax2, args)
            # From https://matplotlib.org/examples/pylab_examples/broken_axis.html
            d = 0.01  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            d2 = 0.01 * (a / b)
            kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False)
            ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # type: ignore[arg-type]
            kwargs = dict(transform=ax2.transAxes, color="k", clip_on=False)
            ax2.plot((-d2, +d2), (-d, +d), **kwargs)  # type: ignore[arg-type]
            # https://stackoverflow.com/questions/42045767/how-can-i-change-the-x-axis-in-matplotlib-so-there-is-no-white-space
        # ***removed and not ontop
        if not args.ontop:
            gs = gridspec.GridSpec(len(args.file), 2, width_ratios=[a / b, 1])
            axislist = []
            for i in range(0, len(args.file) * 2, 2):
                if i == 0:
                    axislist.append(plt.subplot(gs[0]))  # low field (high ppm)
                    axislist.append(plt.subplot(gs[1]))  # high field (low ppm)
                else:
                    axislist.append(plt.subplot(gs[i], sharex=axislist[0]))
                    axislist.append(plt.subplot(gs[i + 1], sharex=axislist[1]))
                if args.orientation[int(i / 2)] != 1:
                    intensit[int(i / 2)] = np.array(intensit[int(i / 2)]) * -1
            for i in range(0, len(args.file) * 2, 2):
                axislist[i].plot(
                    ppm[int(i / 2)],
                    intensit[int(i / 2)],
                    label=args.label[int(i / 2)],
                    color=args.colors[int(i / 2)],
                    linewidth=args.linewidth,
                )
                axdefaultsettings(axislist[i], args)
                axislist[i].set_xlim(
                    left=args.end, right=args.endremove
                )  ### sets x limits must be here
                axislist[i].get_yaxis().set_visible(False)
                axislist[i + 1].plot(
                    ppm[int(i / 2)],
                    intensit[int(i / 2)],
                    label=args.label[int(i / 2)],
                    color=args.colors[int(i / 2)],
                    linewidth=args.linewidth,
                )
                axdefaultsettings(axislist[i + 1], args)
                axislist[i + 1].set_xlim(
                    left=args.startremove, right=args.start
                )  ### sets x limits must be here
                axislist[i + 1].get_yaxis().set_visible(False)
                axislist[i + 1].legend(
                    loc="upper right",
                    bbox_to_anchor=(1.45, 0.6),
                    fancybox=False,
                    shadow=False,
                    frameon=args.keybox,
                    fontsize=args.fontsize,
                )
                if int((i + 2) / 2) == len(args.file):
                    axislist[i].patch.set_alpha(0.0)
                    axislist[i].tick_params(axis="x", bottom=True, labelbottom=True)
                    axislist[i].xaxis.set_tick_params(which="minor", bottom=True)
                    axislist[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
                    axislist[i].spines["bottom"].set_visible(True)
                    axislist[i].spines["bottom"].set_position(
                        ("outward", 1)
                    )  # set spine (in picture the x axis down by x points)
                    axislist[i + 1].patch.set_alpha(0.0)
                    axislist[i + 1].tick_params(axis="x", bottom=True, labelbottom=True)
                    axislist[i + 1].xaxis.set_tick_params(which="minor", bottom=True)
                    axislist[i + 1].xaxis.set_minor_locator(tck.AutoMinorLocator())
                    axislist[i + 1].spines["bottom"].set_visible(True)
                    axislist[i + 1].spines["bottom"].set_position(
                        ("outward", 1)
                    )  # set spine (in picture the x axis down by x points)
            ### make ticks equal
            axislist[0], axislist[1] = equal_ticks(axislist[0], axislist[1], args)

            numfiles = len(args.file)
            # From https://matplotlib.org/examples/pylab_examples/broken_axis.html
            d = (
                0.005 * numfiles
            )  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            d2 = 0.005 * (a / b) * numfiles

            for i in range(0, len(args.file) * 2, 2):
                if args.orientation[int(i / 2)] != 1:
                    kwargs = dict(
                        transform=axislist[i].transAxes, color="k", clip_on=False
                    )
                    axislist[i].plot((1 - d, 1 + d), (1 - d * 3, 1 + d * 3), **kwargs)  # type: ignore[arg-type]
                    axislist[i].patch.set_alpha(0.0)
                    kwargs = dict(
                        transform=axislist[i + 1].transAxes, color="k", clip_on=False
                    )
                    axislist[i + 1].plot((-d2, +d2), (1 - d * 3, 1 + d * 3), **kwargs)  # type: ignore[arg-type]
                    axislist[i + 1].patch.set_alpha(0.0)
                else:
                    kwargs = dict(
                        transform=axislist[i].transAxes, color="k", clip_on=False
                    )
                    axislist[i].plot((1 - d, 1 + d), (-d * 3, +d * 3), **kwargs)  # type: ignore[arg-type]
                    axislist[i].patch.set_alpha(0.0)
                    kwargs = dict(
                        transform=axislist[i + 1].transAxes, color="k", clip_on=False
                    )
                    axislist[i + 1].plot((-d2, +d2), (-d * 3, +d * 3), **kwargs)  # type: ignore[arg-type]
                    axislist[i + 1].patch.set_alpha(0.0)
            if args.orientation[int(len(args.file)) - 1] == -1:
                kwargs = dict(transform=axislist[i].transAxes, color="k", clip_on=False)
                axislist[i].plot((1 - d, 1 + d), (-d * 3, +d * 3), **kwargs)  # type: ignore[arg-type]
                axislist[i].patch.set_alpha(0.0)
                kwargs = dict(
                    transform=axislist[i + 1].transAxes, color="k", clip_on=False
                )
                axislist[i + 1].plot((-d2, +d2), (-d * 3, +d * 3), **kwargs)  # type: ignore[arg-type]
                axislist[i + 1].patch.set_alpha(0.0)

    elif args.startremove is None or args.endremove is None:
        # spectrum not removed!
        print(
            "    Plotting spectrum in range of {} to {} ppm.".format(
                args.start, args.end
            )
        )
        if args.ontop:
            # *** spectra ontop of each other and not removed!***
            if len(args.file) == 1:
                print(
                    "    Plotting ontop with only one file is useless! Going to exit!"
                )
                exit(1)
            print("    Plotting ontop!")
            gs = gridspec.GridSpec(1, 1)
            ax1 = plt.subplot(gs[0])
            maximum = []
            for number in range(0, len(args.file), 1):
                maximum.append(max([yval for yval in intensit[number]]))
                maximum[number] = (
                    maximum[number] + 0.000001
                    if maximum[number] == 0
                    else maximum[number]
                )
                ax1.plot(
                    ppm[number],
                    [yval / maximum[number] for yval in intensit[number]],
                    label=args.label[number],
                    color=args.colors[number],
                    linewidth=args.linewidth,
                )
            axdefaultsettings(ax1, args)
            ax1.autoscale(enable=False, axis="y")
            ax1.set_xlim(
                left=args.end, right=args.start
            )  ### sets x limits must be here
            ax1.tick_params(axis="x", bottom=True, labelbottom=True)
            ax1.xaxis.set_tick_params(which="minor", bottom=True)
            ax1.xaxis.set_minor_locator(tck.AutoMinorLocator())
            ax1.spines["bottom"].set_visible(True)
            # set spine (in picture the x axis down by x points)
            ax1.spines["bottom"].set_position(("outward", 1))
            ax1.get_yaxis().set_visible(False)
            ax1.legend(loc="upper right", frameon=args.keybox, fontsize=args.fontsize)

        if args.stacked:
            # *** teste versetztes plotten ontop (innerhalb eines subplots, sonst wie normales plotten )
            if len(args.file) == 1:
                print(
                    "    Plotting stacked with only one file is useless! Going to exit!"
                )
                exit(1)
            print("    Stacked plotting!")
            gs = gridspec.GridSpec(1, 1)
            ax1 = plt.subplot(gs[0])

            # colors which are opaque
            opaque_colors = [
                (0.027450980392156862, 0.3215686274509804, 0.6039215686274509, 1.0),
                (0.9176470588235294, 0.7254901960784313, 0.047058823529411764, 1.0),
                (0.5647058823529412, 0.5647058823529412, 0.5215686274509804, 1.0),
                (0.0, 0.8, 0.0, 1.0),
                (1.0, 0.0, 1.0, 0.6),
            ]

            maximum = []
            for number in range(0, len(args.file), 1):
                maximum.append(max([yval for yval in intensit[number]]))
                maximum[number] = (
                    maximum[number] + 0.000001
                    if maximum[number] == 0
                    else maximum[number]
                )
                ax1.plot(
                    ppm[number],
                    [
                        (
                            (yval / maximum[number]) / len(args.file)
                            + (1 / len(args.file) * number)
                        )
                        for yval in intensit[number]
                    ],
                    label=args.label[number],
                    color=opaque_colors[number],
                    linewidth=args.linewidth,
                )
            axdefaultsettings(ax1, args)
            ax1.autoscale(enable=False, axis="y")
            ax1.set_xlim(
                left=args.end, right=args.start
            )  ### sets x limits must be here
            ax1.tick_params(axis="x", bottom=True, labelbottom=True)
            ax1.xaxis.set_tick_params(which="minor", bottom=True)
            ax1.xaxis.set_minor_locator(tck.AutoMinorLocator())
            # set spine (in picture the x axis down by x points)
            ax1.spines["bottom"].set_position(("outward", 1))
            ax1.spines["bottom"].set_visible(True)
            ax1.get_yaxis().set_visible(False)
            ax1.legend(loc="upper right", frameon=args.keybox, fontsize=args.fontsize)
        elif not args.ontop:
            # *** spectra not removed and not ontop (normal case)***
            gs = gridspec.GridSpec(len(args.file), 1)
            axislist = []
            for i in range(len(args.file)):
                if i == 0:
                    axislist.append(plt.subplot(gs[0]))
                else:
                    axislist.append(plt.subplot(gs[i], sharex=axislist[0]))
                if args.orientation[i] != 1:
                    intensit[i] = np.array(intensit[i]) * -1
                axislist[i].plot(
                    ppm[i],
                    intensit[i],
                    label=args.label[i],
                    color=args.colors[i],
                    linewidth=args.linewidth,
                )
                axdefaultsettings(axislist[i], args)
                axislist[i].set_xlim(
                    left=args.end, right=args.start
                )  ### sets x limits must be here
                axislist[i].legend(
                    loc="upper right",
                    bbox_to_anchor=(1.1, 0.6),
                    fancybox=False,
                    shadow=False,
                    frameon=args.keybox,
                    fontsize=args.fontsize,
                )

                if i == len(args.file) - 1:
                    axislist[i].tick_params(axis="x", bottom=True, labelbottom=True)
                    axislist[i].xaxis.set_tick_params(which="minor", bottom=True)
                    axislist[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
                    axislist[i].spines["bottom"].set_visible(True)
                    axislist[i].get_yaxis().set_visible(False)
                    axislist[i].spines["bottom"].set_position(
                        ("outward", 1)
                    )  # set spine (in picture the x axis down by x points)

    figure.subplots_adjust(wspace=0.05, hspace=0.05)
    figure.text(0.5, 0.04, r"$\delta$ / ppm", ha="center", fontsize=args.fontsize)
    plt.savefig(args.out + ".pdf", dpi=300)
    plt.savefig(args.out + ".svg")
    print("    Plot is saved to {}.pdf !".format(args.out))
    if input("    Do you want to show matplotlib results?    ") in ("y", "yes"):
        plt.show()
    print("    All Done!")


if __name__ == "__main__":
    main()
