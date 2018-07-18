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

import numpy as np


class varhist():
    def __init__(self, binc: int = 40, rstart: float = 0.0, rstop: float = 1.0, frac: float = 0.6, minbins: float = 1.0,
                 addfun: callable = None):
        self.bins = dict()
        self.floor = rstart
        self.ceiling = rstop
        self.tothits = 0

        if addfun is None:
            self.addfun = self.linadd
        else:
            self.addfun = addfun

        intsize = rstop - rstart

        for i in range(binc - 1):
            if intsize - intsize * frac >= minbins:
                self.bins[(rstart + intsize * frac, rstop)] = 0
                rstop = rstart + intsize * frac
                intsize = intsize * frac
            else:
                binsleft = intsize/minbins
                # We want small bins
                binsleft = int(np.ceil(binsleft))
                binsize = intsize/binsleft
                for j in range(binsleft):
                    self.bins[(rstart + binsize * j, rstart + binsize * (j + 1))] = 0
                break

    def linadd(self, binidx, val):
        self.bins[binidx] += val

    def append(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.addfun(i, 1)
                    self.tothits+=1
                    return True
        return False


class probhist(varhist):
    def __init__(self, binc, rstart, rstop, frac: float = 0.6, minbins: float = 1.0):
        varhist.__init__(self, binc, rstart, rstop, frac, minbins)
        # First element in a tuple is number of staisfactory items
        # Second is total number of hits
        for i in self.bins.keys():
            self.bins[i] = (0, 0)

    def linadd(self, binidx, tuple_val):
        self.bins[binidx] = (self.bins[binidx][0] + tuple_val[0], self.bins[binidx][1] + tuple_val[1])

    def append_succ(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.tothits += 1
                    self.linadd(i, (1, 1))

    def append_fail(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.tothits += 1
                    self.linadd(i, (0, 1))


class anghist():
    def __init__(self, binc: int = 36, rstart: float = 0.0, rstop: float = 360.0, addfun: callable = None):
        self.bins = dict()
        self.floor = rstart
        self.ceiling = rstop
        self.tothits = 0

        if addfun is None:
            self.addfun = self.linadd
        else:
            self.addfun = addfun

        for i in range(binc):
            self.bins[(rstart + i * (rstop - rstart) / binc, rstart + (i + 1) * (rstop - rstart) / binc)] = 0.0

    def linadd(self, binidx, val):
        self.bins[binidx] += val

    def append(self, ang, linpow):
        if ang < self.ceiling and ang > self.floor:
            for i in self.bins.keys():
                if ang > i[0] and ang < i[1]:
                    self.addfun(i, linpow)
                    self.tothits += 1
                    return True
        return False