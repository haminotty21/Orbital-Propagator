from matplotlib.pyplot import xlabel, ylabel
from numpy import array, sqrt, dot, cross, floor, pi, arccos, log2, sin, sinh, arcsin, arccosh, log, linspace, cos, \
    arctan, tan, arcsinh, arctan2, transpose
from numpy.array_api import atan2
from numpy.linalg import norm
from matplotlib import pyplot as plot


def Pork_Chop_Plot(X, Y, Z, levels):
    fig, ax = plot.subplots(constrained_layout=True)
    CS = ax.contour(X, Y, Z, levels)
    line_colors = ['black' for l in CS.levels]
    xlabel('Departure Dates')
    ylabel('Arrival Dates')
    ax.clabel(CS, fontsize=10, colors=line_colors)
    fig.colorbar(CS)
    fig.show()
    ax.set_xlabel('Departure Dates (Julian Days)')
    ax.set_xlabel('Arrival Dates (Julian Days)')

    return fig, ax










