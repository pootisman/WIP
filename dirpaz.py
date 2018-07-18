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
import matplotlib.pyplot as mpl
import pairdata
from auxfun import *
import numba

__author__ = 'Aleksei Ponomarenko-Timofeev'

class rec_pat():
    def __init__(self, source):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def draw(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, print: bool = False):
        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        for i in txrange:
            if self.source.txs[i].setid == txgrp or txgrp == -1:
                rr = 0
                for j in rxrange:
                    if self.source.rxs[j].setid == rxgrp or rxgrp == -1:
                        f = mpl.figure(rr)
                        rr += 1
                        th = []
                        tt = []
                        r = []
                        for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.keys():
                            th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].AoA, l2db(self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow)))
                        th = sorted(th, key=lambda x: x[0])
                        for t in th:
                            r.append(t[1])
                            tt.append(t[0])

                        (tt, r) = basint2(tt, r, 360)
                        ax = mpl.subplot(111, projection='polar')
                        ax.plot(np.deg2rad(tt), r)
                        ax.set_rmax(-100)
                        ax.set_rmin(-200)
                        mpl.title('RX #{}'.format(j))
                        if print:
                            mpl.savefig('Recpat_tx{0:03d}->rx{1:03d}.png'.format(i,j))
                            mpl.close(f)

        if print is False:
            mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/Class@60GHz.TEST_60_MKE_15.sqlite')
    DS.load_paths()
    cir =  rec_pat(DS)
    cir.draw(1, 120, rxgrp=4, print=True)