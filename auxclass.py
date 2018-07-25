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
        if self.floor <= val <= self.ceiling:
            for i in self.bins.keys():
                if i[0] < val < i[1]:
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
        if self.floor <= val <= self.ceiling:
            for i in self.bins.keys():
                if i[0] < val < i[1]:
                    self.tothits += 1
                    self.linadd(i, (1, 1))

    def append_fail(self, val):
        if self.floor <= val <= self.ceiling:
            for i in self.bins.keys():
                if i[0] < val < i[1]:
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
        if self.floor <= ang <= self.ceiling:
            for i in self.bins.keys():
                if i[0] < ang < i[1]:
                    self.addfun(i, linpow)
                    self.tothits += 1
                    return True
        return False


class anghist2():
    def __init__(self, azbinc: int = 36, elbinc: int = 18, azstart: float = -180, azstop: float = 180.0, elstart=0.0,
                 elstop=180.0, addfun: callable = None):
        self.bins = dict()
        self.azfloor = azstart
        self.azceiling = azstop
        self.elfloor = elstart
        self.elceiling = elstop
        self.tothits = 0
        self.azbinc = azbinc
        self.elbinc = elbinc

        if addfun is None:
            self.addfun = self.linadd
        else:
            self.addfun = addfun

        for i in range(elbinc):
            for j in range(azbinc):
                self.bins[(azstart + j * (azstop - azstart) / azbinc, azstart + (j + 1) * (azstop - azstart) / azbinc,
                           elstart + i * (elstop - elstart) / elbinc, elstart + (i + 1) * (elstop - elstart) / elbinc)]\
                    = 0.0

    def linadd(self, binidx, val):
        self.bins[binidx] += val

    def append(self, azang, elang, linpow):
        if self.azfloor <= azang <= self.azceiling and self.elfloor <= elang <= self.elceiling:
            for i in self.bins.keys():
                if i[0] < azang < i[1] and i[2] < elang < i[3]:
                    self.addfun(i, linpow)
                    self.tothits += 1
                    return True
        return False