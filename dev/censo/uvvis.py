# user-customizable variables:
#
# sigma (width of the gaussians)
# scaling factor for the line spectrum (when a line spectrum should be plotted)
# x-range of the plot (x_min, x_max)
# y-max value
#
# steps:
# get nroots from input (hardcoded to be 20, can and should be made user-customizable later)
# get the table with excitations from orca output file where excitation wavelengths, wavenumbers and osc_str are given
import matplotlib.pyplot as plt
import numpy as np

# TODO - add options to set nroots, sigma and x-axis for plot
# TODO - set nroot default value dynamically depending on atom count
# TODO - option for emission spectra
# TODO - add user customizable pyplot settings

# plot spectra for every uvvis calculation
def plot(calculate):

    # set x_range for plot
    x_range = np.linspace(x_min, x_max, 1000) # G

    for conf in calculate:
        plt.plot(x_range, epsilon(x_range, conf.excitations, sigma))
        plt.savefig(...)
        plt.close()


def epsilon(x_range, excitations, sigma):
    """calculate absorption values"""
    # TODO - how to know what should be on x_axis (wl, energy, ...?)
    # each conf: nroot excitations, each with wl and osc_str

    res = np.zeros(len(x_range))

    for exc in excitations:
        res += 1.3062974 * 10**8 * exc["osc_str"] / sigma * np.exp(-((1 / x_range - 1 / exc["wavelength"]) / sigma)**2)