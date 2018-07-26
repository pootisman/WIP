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

__author__ = 'Aleksei Ponomarenko-Timofeev'

class cirs():
    def __init__(self, source):
        self.source = source

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def draw(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, tofile: bool = False,
             cmap: str = 'viridis', xdim: int = 100, ydim: int = 250, zmin: float = -200, nff: bool = True):
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
                self.xdata = []
                self.ydata = []
                self.zdata = []
                f = mpl.figure(i)
                nn = 0
                maxy = 0
                for j in rxrange:
                    if self.source.rxs[j].setid == rxgrp or rxgrp == -1:
                        nn += 1
                        if self.source.txs[i].chan_to(self.source.rxs[j]) is not None:
                            for k in self.source.txs[i].chan_to(self.source.rxs[j]).paths.items():
                                if nff and not k[1].near_field_failed:
                                    self.xdata.append(j)
                                    self.ydata.append(k[1].delay * 1e9)
                                    self.zdata.append(l2db(k[1].pow))
                                elif not nff:
                                    self.xdata.append(j)
                                    self.ydata.append(k[1].delay * 1e9)
                                    self.zdata.append(l2db(k[1].pow))
                        else:
                            if self.ydata.__len__() > 0 and maxy == 0:
                                maxy = np.nanmax(self.ydata)

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
                mpl.tight_layout()
                if tofile:
                    mpl.savefig('CIR3D_tx{0:03d}.png'.format(i))
                    mpl.close(f)

        if tofile is False:
            mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor()
    #DS.load_rxtx('class.sqlite')
    DS.load_rxtx('Human_crawl.TEST.sqlite')
    #DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/Class@60GHz.TEST_60_MKE_15.sqlite')
    #DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TESTe/Class@60GHz.TESTe.sqlite')
    DS.load_paths(npaths=250)
    cir = cirs(DS)
    cir.draw(-1, -1, txgrp=-1, rxgrp=2, tofile=False, ydim=250)
    exit()