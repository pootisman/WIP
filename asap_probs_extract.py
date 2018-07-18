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

__author__ = 'Aleksei Ponomarenko-Timofeev'

class varhist():
    def __init__(self, binc, rstart, rstop, frac: float = 0.6, minbins: float = 1.0):
        self.bins = dict()
        self.floor = rstart
        self.ceiling = rstop
        self.tothits = 0

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

    def append(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.bins[i]+=1
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

    def append_succ(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.tothits += 1
                    self.bins[i] = (self.bins[i][0] + 1, self.bins[i][1] + 1)

    def append_fail(self, val):
        if val < self.ceiling and val > self.floor:
            for i in self.bins.keys():
                if val > i[0] and val < i[1]:
                    self.tothits += 1
                    self.bins[i] = (self.bins[i][0], self.bins[i][1] + 1)


class distanced_hist_extractor():
    def __init__(self, src: pairdata.data_stor, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7, thrs: float = -115, minbins: float = 0.2):
        self.hist = probhist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac, minbins=minbins)
        self.source = src
        self.type = None
        self.mut = False
        self.thresh = thrs

    def has_path(self, src: pairdata.Node, dest: pairdata.Node, typ: str):
        if dest in src.chan_to_pairs:
            # Go over all paths
            for i in src.chan_to_pairs[dest].paths.items():
                if typ == 'LOS' and i[1].interactions.__len__() == 0 and l2db(i[1].pow) >= self.thresh:
                    return True
                elif typ == 'NLOS' and i[1].interactions.__len__() > 0 and l2db(i[1].pow) >= self.thresh:
                    return True
                elif typ == 'NLOS-1' and i[1].interactions.__len__() == 1 and l2db(i[1].pow) >= self.thresh:
                    return True
                elif typ == 'NLOS-2' and i[1].interactions.__len__() == 2 and l2db(i[1].pow) >= self.thresh:
                    return True
                elif typ == 'NLOS-3' and i[1].interactions.__len__() == 3 and l2db(i[1].pow) >= self.thresh:
                    return True
                elif typ == 'noLOS' and i[1].interactions.__len__() == 0 and l2db(i[1].pow) >= self.thresh:
                    return False
                elif typ == 'noNLOS' and i[1].interactions.__len__() > 0 and l2db(i[1].pow) >= self.thresh:
                    return False
                elif typ == 'link' and l2db(i[1].pow) >= self.thresh:
                    return True
                elif typ == 'nolink' and l2db(i[1].pow) >= self.thresh:
                    return False
        else:
            return False

        if typ in ['LOS', 'NLOS', 'link', 'NLOS-1', 'NLOS-2', 'NLOS-3']:
            return False
        else:
            return True

    def build(self, txgrp: int = -1, rxgrp: int = -1, typ: str = 'LOS'):
        self.type = typ
        for i in self.source.txs.items():
            #print(i[0])
            if i[1].setid == txgrp or txgrp == -1:
                for j in i[1].chan_to_pairs.items():
                    #print('{} {}'.format(j[0].node_id, l2db(j[1].pow)))
                    if j[0].setid == rxgrp or rxgrp == -1:
                        if self.has_path(i[1], j[0], typ):
                            self.hist.append_succ(j[1].dist)
                        else:
                            self.hist.append_fail(j[1].dist)

    def plot_hist(self):
        fig = mpl.figure()
        ax = fig.add_subplot(211)
        ax.grid()
        ax.set_xlim([self.hist.floor, self.hist.ceiling])
        ax.set_ylim([0, 1])
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('{} probability'.format(self.type))
        ax.set_title('Total hits {} @ {} dBm threshold'.format(self.hist.tothits, self.thresh))
        bars = []
        tcks = []

        for i in self.hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1][0]/i[1][1] if i[1][1] > 0 else 0.0))
            tcks.append(i[0][0])
            ax.text(i[0][0] + (i[0][1] - i[0][0])/2.0, 0.5, s='{}'.format(i[1][0]), rotation='vertical', horizontalalignment='center', verticalalignment='center')
        tcks.append(self.hist.ceiling)

        mpl.xticks(tcks, rotation='vertical')

        bc = mpc.PatchCollection(bars)

        ax.add_collection(bc)

        ax = fig.add_subplot(212)
        ax.grid()
        ax.set_xlim([self.hist.floor, self.hist.ceiling])
        ax.set_ylim([0, 1])
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('Hit probability'.format(self.type))
        ax.set_title('Total hits {} @ {} dBm threshold'.format(self.hist.tothits, self.thresh))
        bars2 = []

        for i in self.hist.bins.items():
            bars2.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1][1]/self.hist.tothits))
            ax.text(i[0][0] + (i[0][1] - i[0][0])/2.0, 0.5, s='{}'.format(i[1][1]), rotation='vertical', horizontalalignment='center', verticalalignment='center')

        mpl.xticks(tcks, rotation='vertical')

        bc = mpc.PatchCollection(bars2)

        ax.add_collection(bc)
        mpl.tight_layout()
        mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor()
    #DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/HumanCrawl/Human_crawl_X3D_Control/Human_crawl.Human_crawl_X3D_Control.sqlite')
    DS.load_rxtx('Human_crawl.TEST.sqlite')
    DS.load_paths(npaths=75)
    DS.load_interactions(store=False)

    DE = distanced_hist_extractor(DS, range=(0.04, 1.0), histbins=50, frac=0.95, thrs=-95, minbins=0.02)
    DE.build(txgrp=-1, rxgrp=-1, typ='LOS')
    DE.plot_hist()
    exit()