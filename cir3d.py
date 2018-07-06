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

class cirs():
    def __init__(self, source):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def draw(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, print: bool = False, cmap: str = 'viridis', xdim: int = 100, ydim: int = 250, zmin: float = -200):
        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        for i in txrange:
            self.xdata = []
            self.ydata = []
            self.zdata = []
            if self.source.txs[i].setid == txgrp or txgrp == -1:
                f = mpl.figure(i)
                nn = 0
                maxy = 0
                for j in rxrange:
                    if self.source.rxs[j].setid == rxgrp or rxgrp == -1:
                        nn += 1
                        if self.source.txs[i].chan_to(self.source.rxs[j]) is not None:
                            for k in self.source.txs[i].chan_to(self.source.rxs[j]).paths.keys():
                                self.xdata.append(j)
                                self.ydata.append(self.source.txs[i].chan_to(self.source.rxs[j]).paths[k].delay * 1e9)
                                self.zdata.append(l2db(self.source.txs[i].chan_to(self.source.rxs[j]).paths[k].pow))
                        else:
                            if self.ydata.__len__() > 0 and maxy == 0:
                                maxy = np.max(self.ydata)

                            for k in range(ydim):
                                self.xdata.append(j)
                                self.ydata.append(maxy)
                                self.zdata.append(np.NaN)

                if np.max(self.ydata) == 0:
                    self.ydata[-1] = 1e-9

                if np.isnan(np.nanmin(self.zdata)):
                    for z in range(self.zdata.__len__()):
                        self.zdata[z] = zmin

                (X, Y, Z) = basint3(self.xdata, self.ydata, self.zdata, nn, ydim)
                [X, Y] = np.meshgrid(X, Y)
                mpl.contourf(np.transpose(X), np.transpose(Y), Z, 20, cmap=cmap)
                mpl.clim([np.min(Z), np.max(Z) + np.abs(0.1 * np.max(Z))])
                mpl.colorbar().set_label('RX power, [dBm]')
                mpl.xlabel('RX Position')
                mpl.ylabel('Delay, [ns]')
                mpl.title('CIR@TX #{}'.format(i))
                if print:
                    mpl.savefig('CIR3D_tx{0:03d}.png'.format(i))
                    mpl.close(f)

        if print is False:
            mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TESTe/Class@60GHz.TESTe.sqlite')
    DS.load_path()
    cir = cirs(DS)
    cir.draw(-1, -1, txgrp=1, rxgrp=4, print=True)