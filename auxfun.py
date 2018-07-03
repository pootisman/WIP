# Copyright (C) Aleksei Ponomarenko-Timofeev
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

__author__ = 'Aleksei Ponomarenko-Timofeev'

import numpy as np
import numba

@numba.jit
def basint2(X, Y, xc):
    Xlim = [np.min(X), np.max(X)]

    stepX = 360.0 / xc

    Xo = np.arange(start=-180, step=stepX, stop=180)
    Yo = np.tile(np.min(Y), [xc])

    for i in range(Y.__len__()):
        j = int(np.round((X[i] - Xlim[0]) / stepX))
        if j == xc:
            j = j - 1

        if Yo[j] < Y[i]:
            Yo[j] = Y[i]

    return (Xo, Yo)

@numba.jit
def basint3(X, Y, Z, xc, yc):
    Xlim = [np.min(X), np.max(X)]
    Ylim = [np.min(Y), np.max(Y)]

    stepX = (Xlim[1] - Xlim[0]) / xc
    stepY = (Ylim[1] - Ylim[0]) / yc

    Xo = np.arange(start=Xlim[0], step=stepX, stop=Xlim[1])
    Yo = np.arange(start=Ylim[0], step=stepY, stop=Ylim[1])
    Zo = np.tile(np.min(Z), [xc, yc])

    for i in range(Z.__len__()):
        j = int(np.round((X[i] - Xlim[0]) / stepX))
        k = int(np.round((Y[i] - Ylim[0]) / stepY))
        if j == xc:
            j = j - 1
        if k == yc:
            k = k - 1

        if Zo[j, k] < Z[i]:
            Zo[j, k] = Z[i]

    return (Xo, Yo, Zo)

@numba.jit
def l2db(val: float):
    return 10.0 * np.log10(val)

@numba.jit
def db2l(val: float):
    return np.power(10.0, val / 10.0)
