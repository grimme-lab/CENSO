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

# plot spectra for every uvvis calculation
plt.ion()
plt.plot(...)
plt.ioff()
plt.savefig(...)
plt.close(...)