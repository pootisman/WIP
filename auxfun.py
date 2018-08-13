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


def basint2(X: list, Y: list, xc: int):
    Xlim = [np.nanmin(X), np.nanmax(X)]

    Xo, stepX = np.linspace(start=Xlim[0], stop=Xlim[1], num=xc, retstep=True)
    Yo = np.tile(np.nanmin(Y), [Xo.__len__()])

    for i in range(Y.__len__()):
        j = int(np.round((X[i] - Xlim[0]) / stepX))
        if j == xc:
            j = j - 1

        if Yo[j] < Y[i]:
            if not np.isnan(Y[i]):
                Yo[j] = Y[i]
            else:
                Yo[j] = np.nanmin(Y)

    return Xo, Yo

def basint3(X: list, Y: list, Z: list, xc: float, yc: float, zmin: float = np.nan):
    Xlim = [np.nanmin(X), np.nanmax(X)]
    Ylim = [np.nanmin(Y), np.nanmax(Y)]

    Xo, stepX = np.linspace(start=Xlim[0], stop=Xlim[1], num=xc, retstep=True)
    Yo, stepY = np.linspace(start=Ylim[0], stop=Ylim[1], num=yc, retstep=True)

    Zo = np.tile(np.nanmin(Z) if np.isnan(zmin) else zmin, [Xo.__len__(), Yo.__len__()])

    for i in range(Z.__len__()):
        if stepX != 0:
            j = int(np.round((X[i] - Xlim[0]) / stepX))
        else:
            j = Xlim[0]

        if stepY != 0:
            k = int(np.round((Y[i] - Ylim[0]) / stepY))
        else:
            k = Ylim[0]

        if j == Xo.__len__():
            j = j - 1
        if k == Yo.__len__():
            k = k - 1

        if Zo[j, k] < Z[i]:
            if not np.isnan(Z[i]):
                Zo[j, k] = Z[i]
            else:
                Zo[j, k] = np.nanmin(Z)

    return Xo, Yo, Zo

def l2db(val: float):
    return 10.0 * np.log10(val)

def db2l(val: float):
    return np.power(10.0, val / 10.0)