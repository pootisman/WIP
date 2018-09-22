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

from pairdata import DataStorage, Node
from auxclass import ProbHist
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.patches as mpp
import matplotlib.collections as mpc

__author__ = 'Aleksei Ponomarenko-Timofeev'


class DistancedHistExtractor:
    def __init__(self, src: DataStorage, histbins: int = 10, range: tuple = (1, 16), frac: float = 0.7,
                 thrs: float = -115, minbins: float = 0.2, nffilt: bool = True):
        self.hist = ProbHist(binc=histbins, rstart=range[0], rstop=range[1], frac=frac, minbins=minbins)
        self.source = src
        self.type = None
        self.trans_type = None
        self.mut = False
        self.thresh = thrs
        self.rx_proc = {Node: list}
        self.nffilt = nffilt

    def has_path(self, src: Node, dest: Node, typ: str):
        if dest in src.chans_to_pairs:
            # Go over all paths
            for i in src.chans_to_pairs[dest].paths.items():
                if not (i[1].near_field_failed and self.nffilt) or not self.nffilt:
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
                    elif typ == 'any':
                        return True
                else:
                    print('NF test failed, ignoring path {} in chan {}->{}'.format(i[1].pathid, src.node_id,
                                                                                   dest.node_id))
                    #pass
        else:
            return False

        if typ in ['LOS', 'NLOS', 'link', 'NLOS-1', 'NLOS-2', 'NLOS-3']:
            return False
        else:
            return True

    def build(self, txgrp: int = -1, rxgrp: int = -1, typ: str = 'LOS'):
        self.type = typ
        for i in self.source.txs.items():
            if i[1].setid == txgrp or txgrp == -1:
                for j in i[1].chans_to_pairs.items():
                    if j[0].setid == rxgrp or rxgrp == -1:
                        if self.has_path(i[1], j[0], typ):
                            self.hist.append_succ(j[1].dist)
                        else:
                            self.hist.append_fail(j[1].dist)

    def build_trans(self, txgrp: int = -1, rxgrp: int = -1, typ: str = 'LOS->LOS'):
        self.type = typ
        start, stop = typ.split('->')
        for i in self.source.txs.items():
            if i[1].setid == txgrp or txgrp == -1:
                for j in i[1].chans_to_pairs.keys():
                    if j.setid == rxgrp or rxgrp == -1:
                        if self.has_path(i[1], j, start):
                            self.build_delta(ctx=i[1], crx=j, trans_typ=stop)
                        else:
                            self.build_delta(ctx=i[1], crx=j, trans_typ=False)

    '''Builds delta histogram for one point, called multiple times for eac individual TX'''
    def build_delta(self, ctx: Node, crx: Node, trans_typ: str = 'LOS'):
        if ctx not in self.rx_proc.keys():
            self.rx_proc[ctx] = []

        self.rx_proc[ctx].append(crx)

        for i in ctx.chans_to_pairs.keys():
            if i not in self.rx_proc[ctx]:
                shift = np.linalg.norm(crx.coords - i.coords)
                if isinstance(trans_typ, str):
                    if self.has_path(ctx, i, trans_typ):
                        self.hist.append_succ(shift)
                    else:
                        self.hist.append_fail(shift)
                elif trans_typ:
                    self.hist.append_fail(shift)

    def plot_hist(self, log: bool = False):
        fig = mpl.figure()
        ax = fig.add_subplot(211)
        ax.grid(linestyle='--')
        ax.set_xlim(self.hist.floor, self.hist.ceiling)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('{} probability'.format(self.type))
        ax.set_title('Total hits {} @ {} dBm threshold'.format(self.hist.tothits, self.thresh))
        bars = []
        tcks = []

        for i in self.hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1][0]/i[1][1] if i[1][1] > 0 else 0.0))
            tcks.append(i[0][0])
            ax.text(i[0][0] + (i[0][1] - i[0][0])/2.0, 0.5, s='{}'.format(i[1][0]), rotation='vertical',
                    horizontalalignment='center', verticalalignment='center')
        tcks.append(self.hist.ceiling)

        mpl.xticks(tcks, rotation='vertical')

        bc = mpc.PatchCollection(bars)

        ax.add_collection(bc)

        ax = fig.add_subplot(212)
        ax.grid(linestyle='--')
        ax.set_xlim(self.hist.floor, self.hist.ceiling)
        ax.set_ylim(0, 1)
        if log:
            mpl.semilogy()
        ax.set_xlabel('Distance, [m]')
        ax.set_ylabel('Hit probability'.format(self.type))
        ax.set_title('Total hits {} @ {} dBm threshold'.format(self.hist.tothits, self.thresh))
        bars2 = []

        for i in self.hist.bins.items():
            bars2.append(mpp.Rectangle((i[0][0], 0), i[0][1] -  i[0][0], i[1][1]/self.hist.tothits))
            ax.text(i[0][0] + (i[0][1] - i[0][0])/2.0, 0.5, s='{}'.format(i[1][1]), rotation='vertical',
                    horizontalalignment='center', verticalalignment='center')

        mpl.xticks(tcks, rotation='vertical')

        bc = mpc.PatchCollection(bars2)

        ax.add_collection(bc)
        mpl.tight_layout()
        mpl.show()


if __name__ == "__main__":
    DS = DataStorage(conf='dbconf.txt', threaded=True)
    DS.load_rxtx(dbname='Human_sitting_legsback_Sitting_sqlite')
    DS.load_paths(npaths=250)
    DS.load_interactions(store=True)
    #DS.dump_paths(csvsav=True)

    from phys_path_procs import *
    #check_data_NF(DS)

    DE = DistancedHistExtractor(DS, range=(0.0, 1.0), histbins=50, frac=0.95, thrs=-95, minbins=0.01, nffilt=False)
    DA = DistancedHistExtractor(DS, range=(0.0, 1.0), histbins=50, frac=0.95, thrs=-95, minbins=0.01, nffilt=False)

    DA.build(txgrp=-1,rxgrp=-1, typ='any')
    DE.build_trans(txgrp=-1, rxgrp=-1, typ='LOS->LOS')
    #gen_data_clusters(DS, threshold=0.00001, nff=True)
    DE.plot_hist(log=False)
    DA.plot_hist(log=False)
    exit()
