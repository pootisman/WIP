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

    def draw(self, txrange, rxrange, txgrp: int = -1, rxgrp: int = -1, print: bool = False, cmap: str = 'viridis'):
        for i in range(txrange):
            if self.source.txs[i].setid == txgrp or txgrp == -1:
                f = mpl.figure(i)
                nn = 0

                for j in range(rxrange):
                    if self.source.rxs[j].setid == rxgrp or rxgrp == -1:
                        nn += 1
                        for k in self.source.txs[i].chan_to_pairs[j].paths.keys():
                            self.xdata.append(j)
                            self.ydata.append(self.source.txs[i].chan_to_pairs[j].paths[k].delay * 1e9)
                            self.zdata.append(10.0 * np.log10(self.source.txs[i].chan_to_pairs[j].paths[k].pow))

                (X, Y, Z) = basint3(self.xdata, self.ydata, self.zdata, nn, 250)
                [X, Y] = np.meshgrid(X, Y)
                mpl.contourf(np.transpose(X), np.transpose(Y), Z, 40)
                mpl.clim([np.min(self.zdata), np.max(self.zdata) + 0.1 * np.max(self.zdata)])
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
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/Class@60GHz.TEST_60_MKE_15.sqlite')
    DS.load_path()
    cir = cirs(DS)
    cir.draw(1, 120, rxgrp=4, print=True)