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
import threading as trd
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
                    return

class probhist(varhist):
    def __init__(self, binc, rstart, rstop, frac: float = 0.6):
        varhist.__init__(self, binc, rstart, rstop, frac)
        self.bins.clear()
        intsize = rstop - rstart
        # First element in a tuple is number of staisfactory items
        # Second is total number of hits
        for i in range(binc - 1):
            self.bins[(rstart + intsize * frac, rstop)] = (0, 0)
            rstop = rstart + intsize * frac
            intsize = intsize * frac

    def append(self, val, cond: callable = None, kwargs: dict = None):
        assert cond != None and kwargs != None
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.tothits += 1
                    if cond(kwargs['val']):
                        self.bins[i] = (self.bins[i][0] + 1, self.bins[i][1] + 1)
                        return True
                    else:
                        self.bins[i] = (self.bins[i][0], self.bins[i][1] + 1)
                        return False

class distanced_hist_extractor():
    def __init__(self, src: pairdata.data_stor, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7, thrs: float = -115):
        self.hist = probhist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac)
        self.source = src
        self.type = None
        self.mut = False
        self.thresh = thrs

    def sig_more(self, val):
        return l2db(val) > self.thresh

    def sig_less(self, val):
        return l2db(val) < self.thresh

    @numba.jit
    def build(self, txgrp: int = -1, rxgrp: int = -1, typ: str = 'LOS'):
        self.type = typ
        for i in self.source.txs.items():
            if i[1].setid == txgrp or txgrp == -1:
                for j in i[1].chan_to_pairs.items():
                    if j[1].dest.setid == rxgrp or rxgrp == -1:
                        for k in j[1].paths.items():
                            # No interactions indicate LOS link, check if it fits
                            if k[1].interactions.__len__() == 0 and typ == 'LOS':
                                self.hist.append(j[1].dist, cond=self.sig_more, kwargs={'self': self, 'val': k[1].pow})
                            elif k[1].interactions.__len__() > 0 and typ == 'NLOS':
                                self.hist.append(j[1].dist, cond=self.sig_more, kwargs={'self': self, 'val': k[1].pow})
                            elif k[1].interactions.__len__() == 0 and typ == 'noLOS':
                                self.hist.append(j[1].dist, cond=self.sig_less, kwargs={'self': self, 'val': k[1].pow})
                            elif typ == 'nolink':
                                self.hist.append(j[1].dist, cond=self.sig_less, kwargs={'self': self, 'val': k[1].pow})
                            elif typ == 'link':
                                self.hist.append(j[1].dist, cond=self.sig_more, kwargs={'self': self, 'val': k[1].pow})

        #for i in sorted(self.hist.bins.items(), key=lambda t: t[0][0]):
        #    self.hist.bins[i[0]] = self.hist.bins[i[0]]

    def CDF_mutate(self):
        cumulator = 0.0
        self.mut = True
        for i in sorted(self.hist.bins.items(), key=lambda t: t[0][0]):
            cumulator += self.hist.bins[i[0]]
            self.hist.bins[i[0]] = cumulator

    def plot_hist(self):
        fig = mpl.figure()
        ax = fig.add_subplot(111)
        ax.grid()
        ax.set_xlim([self.hist.floor, self.hist.ceiling])
        if self.mut:
            ax.set_ylim([0, 1])
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('{} probability'.format(self.type))
        ax.set_title('Total hits {}'.format(self.hist.tothits))
        bars = []
        tcks = []

        for i in self.hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1][0]/i[1][1] if i[1][1] > 0 else 0.0))
            tcks.append(i[0][0])
            ax.text(i[0][0] + (i[0][1] - i[0][0]), 0.5, s='{}'.format(i[1][1]))
        tcks.append(self.hist.ceiling)

        mpl.xticks(tcks, rotation='vertical')

        bc = mpc.PatchCollection(bars)

        ax.add_collection(bc)

        mpl.tight_layout()
        mpl.show()


class distanced_delta_hist_extractor():
    def __init__(self, src: pairdata.data_stor, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7):
        self.hist = varhist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac)
        self.source = src
        self.type = None
        self.thread_pool = []
        self.mut = False

    @numba.jit
    def hist_keying(self, path: pairdata.path, typ: str = 'LOS', thresh: float = -95.0, crx: pairdata.Node = None):
        if path.interactions.__len__() == 0 and l2db(path.pow) > thresh and typ == 'LOS':
            self.hist.append(np.linalg.norm(path.chan.src.coords - crx.coords))
            return True
        elif path.interactions.__len__() > 0 and l2db(path.pow) > thresh and typ == 'NLOS':
            self.hist.append(np.linalg.norm(path.chan.src.coords - crx.coords))
            return True
        elif path.interactions.__len__() == 0 and l2db(path.pow) < thresh and typ == 'noLOS':
            self.hist.append(np.linalg.norm(path.chan.src.coords - crx.coords))
            return True
        elif path.interactions.__len__() >= 0 and l2db(path.pow) < thresh and typ == 'noNLOS':
            self.hist.append(np.linalg.norm(path.chan.src.coords - crx.coords))
            return True
        elif path.interactions.__len__() >= 0 and l2db(path.pow) > thresh and typ == 'link':
            self.hist.append(np.linalg.norm(path.chan.src.coords - crx.coords))
            return True
        elif path.interactions.__len__() >= 0 and l2db(path.pow) < thresh and typ == 'nolink':
            self.hist.append(np.linalg.norm(path.chan.src.coords - crx.coords))
            return True

        return False

    @numba.jit
    def build(self, txgrp: int = -1, rxgrp: int = -1, thresh: float = -95, start_typ: str = 'LOS', sw_typ: str = 'LOS'):
        excluderx = []
        self.type = '{}->{}'.format(start_typ, sw_typ)

        # Iterate over all TXs
        for i in self.source.txs.items():
            # Do we filter TXs? What group TX is in?
            if i[1].setid == txgrp or txgrp == -1:
                print('TX[{}]'.format(i[0]))
                # Iterate over channels of the TX
                for j in i[1].chan_to_pairs.items():
                    # Exclude destination of each channel from stats, avoid double counting!
                    excluderx.append(j[1].dest)
                    # Do we filter destinations? What group dest has?
                    if j[1].dest.setid == rxgrp or rxgrp == -1:
                        print('CHAN[{}->{}]'.format(j[1].src.node_id, j[1].dest.node_id))
                        # Iterate over all paths in a channel
                        for k in j[1].paths.items():
                            # No interactions and power over threshold -> LOS
                            if k[1].interactions.__len__() == 0 and l2db(k[1].pow) > thresh and start_typ == 'LOS':
                                # Iterate over all RXs
                                for l in self.source.rxs.items():
                                    # Check RX group and if it was already counted
                                    if (l[1].setid == rxgrp or rxgrp == -1) and l[1] not in excluderx:
                                        # Iterate over all paths
                                        for m in l[1].chan_to_pairs[i[1]].paths.items():
                                            if self.hist_keying(crx=j[1].dest, path=m[1], typ=sw_typ, thresh=thresh):
                                                break
                                            print('Analyzing {}->{} path transfer from {} to {}'.format(start_typ, sw_typ, k[1].pathid, m[1].pathid))
                            elif k[1].interactions.__len__() > 0 and l2db(k[1].pow) > thresh and start_typ == 'NLOS':
                                for l in self.source.rxs.items():
                                    if (l[1].setid == rxgrp or rxgrp == -1) and l[1] not in excluderx:
                                        # Iterate over all paths
                                        for m in l[1].chan_to_pairs[i[1]].paths.items():
                                            if self.hist_keying(crx=j[1].dest, path=m[1], typ=sw_typ, thresh=thresh):
                                                break
                                            print('Analyzing {}->{} path transfer from {} to {}'.format(start_typ, sw_typ, k[1].pathid, m[1].pathid))
                            elif k[1].interactions.__len__() == 0 and l2db(k[1].pow) < thresh and start_typ == 'noLOS':
                                for l in self.source.rxs.items():
                                    if (l[1].setid == rxgrp or rxgrp == -1) and l[1] not in excluderx:
                                        # Iterate over all paths
                                        for m in l[1].chan_to_pairs[i[1]].paths.items():
                                            if self.hist_keying(crx=j[1].dest, path=m[1], typ=sw_typ, thresh=thresh):
                                                break
                                            print('Analyzing {}->{} path transfer from {} to {}'.format(start_typ, sw_typ, k[1].pathid, m[1].pathid))
                            elif l2db(k[1].pow) < thresh and start_typ == 'nolink':
                                for l in self.source.rxs.items():
                                    if (l[1].setid == rxgrp or rxgrp == -1) and l[1] not in excluderx:
                                        # Iterate over all paths
                                        for m in l[1].chan_to_pairs[i[1]].paths.items():
                                            if self.hist_keying(crx=j[1].dest, path=m[1], typ=sw_typ, thresh=thresh):
                                                break
                                            print('Analyzing {}->{} path transfer from {} to {}'.format(start_typ, sw_typ, k[1].pathid, m[1].pathid))
                            elif l2db(k[1].pow) > thresh and start_typ == 'link':
                                for l in self.source.rxs.items():
                                    if (l[1].setid == rxgrp or rxgrp == -1) and l[1] not in excluderx:
                                        # Iterate over all paths
                                        for m in l[1].chan_to_pairs[i[1]].paths.items():
                                            if self.hist_keying(crx=j[1].dest, path=m[1], typ=sw_typ, thresh=thresh):
                                                break
                                            print('Analyzing {}->{} path transfer from {} to {}'.format(start_typ, sw_typ, k[1].pathid, m[1].pathid))

        for i in sorted(self.hist.bins.items(), key=lambda t: t[0][0]):
            self.hist.bins[i[0]] = self.hist.bins[i[0]] / self.hist.tothits

    def CDF_mutate(self):
        self.mut = True
        cumulator = 0.0
        for i in sorted(self.hist.bins.items(), key=lambda t: t[0][0]):
            cumulator += self.hist.bins[i[0]]
            self.hist.bins[i[0]] = cumulator

    def plot_hist(self):
        fig = mpl.figure()
        ax = fig.add_subplot(111)
        ax.grid()
        ax.set_xlim([self.hist.floor, self.hist.ceiling])
        if self.mut:
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

        mpl.xticks(tcks, rotation='vertical')

        bc = mpc.PatchCollection(bars)

        ax.add_collection(bc)

        mpl.tight_layout()
        mpl.show()

    class probplot():
        def __init__(self, src: pairdata.data_stor, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7):
            self.hist = varhist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac)
            self.source = src
            self.type = None
            self.thread_pool = []
            self.mut = False




if __name__ == "__main__":
    DS = pairdata.data_stor()
    #DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/HumanCrawl/Human_crawl_X3D_Control/Human_crawl.Human_crawl_X3D_Control.sqlite')
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TESTe/Class@60GHz.TESTe.sqlite')

    DS.load_path()

    DE = distanced_hist_extractor(DS, range=(0.08, 18), histbins=50, frac=0.85, thrs=-125)
    DE.build(txgrp=-1, rxgrp=-1, typ='LOS')
    #    DE.CDF_mutate()
    DE.plot_hist()

    #DDE = distanced_delta_hist_extractor(DS, range=(0, 18), histbins=50, frac=0.85)
    #DDE.build(txgrp=-1, rxgrp=-1, thresh=-115, start_typ='nolink', sw_typ='link')
    #DDE.CDF_mutate()
    #DDE.plot_hist()
    exit()