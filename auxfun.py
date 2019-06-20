# Copyright (C) Aleksei Ponomarenko-Timofeev
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__author__ = 'Aleksei Ponomarenko-Timofeev'

import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage.filters import uniform_filter

colors = ['m', 'k', 'r', 'g', 'b']
markers = ['*', '.', '+', '^', 'o', 'd', '1', '2', '3', '4', '8']
styles = [':', '-', '-.', '--']


def enable_latex(pt: float = 12.0):
    plt.rc('text', usetex=True)
    font = {'family': 'serif', 'size': pt, 'serif': ['Latin Modern Roman']}
    plt.rc('font', **font)

def basint2(x: list, y: list, xc: int):
    xlim = [np.nanmin(x), np.nanmax(x)]

    xo, step_x = np.linspace(start=float(xlim[0]), stop=float(xlim[1]), num=xc, retstep=True)
    yo = np.tile(np.nanmin(y), [xo.__len__()])

    for i in range(y.__len__()):
        j = int(np.round((x[i] - xlim[0]) / step_x))
        if j == xc:
            j = j - 1

        if yo[j] < y[i]:
            if not np.isnan(y[i]):
                yo[j] = y[i]
            else:
                yo[j] = np.nanmin(y)

    return xo, yo

def basint3(x: list, y: list, z: list, xc: int, yc: int, xmin: float = np.nan, xmax: float = np.nan,
            ymin: float = np.nan, ymax: float = np.nan, zmin: float = np.nan):
    if np.isnan(xmin) or np.isnan(xmax):
        xlim = [np.nanmin(x), np.nanmax(x)]
    else:
        xlim = [xmin, xmax]

    if np.isnan(ymin) or np.isnan(ymax):
        ylim = [np.nanmin(y), np.nanmax(y)]
    else:
        ylim = [ymin, ymax]

    xo, step_x = np.linspace(start=float(xlim[0]), stop=float(xlim[1]), num=xc, retstep=True)
    yo, step_y = np.linspace(start=float(ylim[0]), stop=float(ylim[1]), num=yc, retstep=True)

    zo = np.tile(np.nanmin(z) if np.isnan(zmin) else zmin, [xo.__len__(), yo.__len__()])

    for i in range(z.__len__()):
        if step_x != 0:
            j = int(np.round((x[i] - xlim[0]) / step_x))
        else:
            j = xlim[0]

        if step_y != 0:
            k = int(np.round((y[i] - ylim[0]) / step_y))
        else:
            k = ylim[0]

        if j == xo.__len__():
            j = j - 1
        if k == yo.__len__():
            k = k - 1

        if zo[j, k] < z[i]:
            if not np.isnan(z[i]):
                zo[j, k] = z[i]
            else:
                zo[j, k] = np.nanmin(z)

    return xo, yo, zo

def square_up(array: np.ndarray, smooth: bool = False):
    if not smooth:
        return np.repeat(np.repeat(array, 2, axis=0).T, 2, axis=0).T
    else:
        return uniform_filter(np.repeat(np.repeat(array, 2, axis=0).T, 2, axis=0).T, size=(2, 2))

def l2db(val: float):
    return 10.0 * np.log10(val)

def db2l(val: float):
    return np.power(10.0, val / 10.0)

def excl_interactions(interactions: list = list()):
    c = 0
    for i in interactions:
        if i.typ not in [4]:
            c += 1

    return c
