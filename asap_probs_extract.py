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

import pairdata
from auxfun import *
import matplotlib.pyplot as mpl
import matplotlib.patches as mpp
import matplotlib.collections as mpc
import numba

__author__ = 'Aleksei Ponomarenko-Timofeev'

class varhist():
    def __init__(self, binc, rstart, rstop, frac: float = 0.6):
        self.bins = dict()
        self.samps = 0
        self.floor = rstart
        self.ceiling = rstop
        self.tothits = 0

        intsize = rstop - rstart

        for i in range(binc - 1):
            self.bins[(rstart + intsize * frac, rstop)] = 0
            rstop = rstart + intsize * frac
            intsize = intsize * frac


    def append(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.bins[i]+=1
                    self.tothits+=1

class distanced_hist_extractor():
    def __init__(self, src: pairdata.data_stor, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7):
        self.hist = varhist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac)
        self.source = src
        self.type = None

    #@numba.jit
    def build(self, txgrp: int = -1, rxgrp: int = -1, thresh: float = -95, typ: str = 'LOS'):
        self.type = typ
        for i in self.source.txs.items():
            if i[1].setid == txgrp or txgrp == -1:
                for j in i[1].chan_to_pairs.items():
                    if j[1].dest.setid == rxgrp or rxgrp == -1:
                        for k in j[1].paths.items():
                            # No interactions indicate LOS link
                            if k[1].interactions.__len__() == 0 and l2db(k[1].pow) > thresh and typ == 'LOS':
                                self.hist.append(j[1].dist)
                            elif k[1].interactions.__len__() > 0 and l2db(k[1].pow) > thresh and typ == 'NLOS':
                                self.hist.append(j[1].dist)
                            elif k[1].interactions.__len__() == 0 and l2db(k[1].pow) < thresh and typ == 'noLOS':
                                self.hist.append(j[1].dist)
                            elif l2db(k[1].pow) < thresh and typ == 'nolink':
                                self.hist.append(j[1].dist)
                            elif l2db(k[1].pow) > thresh and typ == 'link':
                                self.hist.append(j[1].dist)

    def CDF_mutate(self):
        cumulator = 0.0
        for i in sorted(DE.hist.bins.items(), key=lambda t: t[0][0]):
            self.hist.bins[i[0]] = self.hist.bins[i[0]] / self.hist.tothits
            cumulator += self.hist.bins[i[0]]
            self.hist.bins[i[0]] = cumulator

    def plot_hist(self):
        fig = mpl.figure()
        ax = fig.add_subplot(111)
        ax.grid()
        ax.set_xlim([self.hist.floor, self.hist.ceiling])
        ax.set_ylim([0, 1])
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('{} probability'.format(self.type))
        ax.set_title('Total hits {}'.format(self.hist.tothits))
        bars = []
        tcks = []

        for i in self.hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1]))
            tcks.append(i[0][0])
        tcks.append(self.hist.ceiling)

        mpl.xticks(tcks)

        bc = mpc.PatchCollection(bars)

        ax.add_collection(bc)

        mpl.show()


class distanced_delta_hist_extractor():
    def __init__(self, src: pairdata.data_stor, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7):
        self.hist = varhist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac)
        self.source = src
        self.type = None

    def hist_keying(self, path: pairdata.path, typ: str = 'LOS', thresh: float = -95.0, crx: pairdata.Node = None,
                    nextf: callable = None, **kwargs):
        if path.interactions.__len__() == 0 and l2db(path.pow) > thresh and typ == 'LOS':
            self.hist.append(np.linalg.norm(path.chan.dest.coords - crx.coords))
        elif path.interactions.__len__() > 0 and l2db(path.pow) > thresh and typ == 'NLOS':
            self.hist.append(np.linalg.norm(path.chan.dest.coords - crx.coords))
        elif path.interactions.__len__() == 0 and l2db(path.pow) < thresh and typ == 'noLOS':
            self.hist.append(np.linalg.norm(path.chan.dest.coords - crx.coords))
        elif path.interactions.__len__() > 0 and l2db(path.pow) < thresh and typ == 'noNLOS':
            self.hist.append(np.linalg.norm(path.chan.dest.coords - crx.coords))

        if nextf is not None:
            nextf(kwargs)

    #@numba.jit
    def build(self, txgrp: int = -1, rxgrp: int = -1, thresh: float = -95, start_typ: str = 'LOS', sw_typ: str = 'LOS'):
        excluderx = []
        self.type = '{}->{}'.format(start_typ, sw_typ)

        for i in self.source.txs.items():
            if i[1].setid == txgrp or txgrp == -1:
                for j in i[1].chan_to_pairs.items():
                    excluderx.append(j[1].dest)
                    if j[1].dest.setid == rxgrp or rxgrp == -1:
                        for k in j[1].paths.items():
                            # No interactions indicate LOS link
                            if k[1].interactions.__len__() == 0 and l2db(k[1].pow) > thresh and start_typ == 'LOS':
                                for l in self.source.rxs.items():
                                    if (l[1].setid == rxgrp or rxgrp == -1) and l[1] not in excluderx:
                                        self.hist_keying(crx=j[1].dest, path=k[1], typ=sw_typ, thresh=thresh)
                                        print('Analyzing transfer from {} to {}'.format(j[1].dest.node_id, l[0]))
                            elif k[1].interactions.__len__() > 0 and l2db(k[1].pow) > thresh and start_typ == 'NLOS':
                                self.hist.append(j[1].dist)
                            elif k[1].interactions.__len__() == 0 and l2db(k[1].pow) < thresh and start_typ == 'noLOS':
                                self.hist.append(j[1].dist)
                            elif l2db(k[1].pow) < thresh and start_typ == 'nolink':
                                self.hist.append(j[1].dist)
                            elif l2db(k[1].pow) > thresh and start_typ == 'link':
                                self.hist.append(j[1].dist)

    def CDF_mutate(self):
        cumulator = 0.0
        for i in sorted(DE.hist.bins.items(), key=lambda t: t[0][0]):
            self.hist.bins[i[0]] = self.hist.bins[i[0]] / self.hist.tothits
            cumulator += self.hist.bins[i[0]]
            self.hist.bins[i[0]] = cumulator


    def plot_hist(self):
        fig = mpl.figure()
        ax = fig.add_subplot(111)
        ax.grid()
        ax.set_xlim([self.hist.floor, self.hist.ceiling])
        ax.set_ylim([0, 1])
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('{} probability'.format(self.type))
        ax.set_title('Total hits {}'.format(self.hist.tothits))
        bars = []
        tcks = []

        for i in self.hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1]))
            tcks.append(i[0][0])
        tcks.append(self.hist.ceiling)

        mpl.xticks(tcks)

        bc = mpc.PatchCollection(bars)

        ax.add_collection(bc)

        mpl.show()




if __name__ == "__main__":
    DS = pairdata.data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/Mobile_TXRX/Class@60GHz.Mobile_TXRX.sqlite')
    DS.load_path()
    DE = distanced_hist_extractor(DS, range=(0, 18), histbins=50, frac=0.85)
    DE.build(txgrp=-1, rxgrp=4, thresh=-125, typ='LOS')
    DE.CDF_mutate()
    DE.plot_hist()
    exit()