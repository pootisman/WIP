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
import scipy.io as sio
import matplotlib.pyplot as mpl
import pairdata
from phys_path_procs import *
from auxclass import powhist

__author__ = 'Aleksei Ponomarenko-Timofeev'

class cirs():
    def __init__(self, source):
        self.source = source

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, topng: bool = False,
             cmap: str = 'viridis', xdim: int = 100, ydim: int = 250, zmin: float = -200.0, zmax: float = np.nan,
            nff: bool = True, matsav: bool = False, plot: bool = True):

        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        for i in txrange:
            if self.source.txs[i].setid in txgrp or txgrp[0] == -1:
                self.xdata = []
                self.ydata = []
                self.zdata = []

                nn = 0
                for j in rxrange:
                    if self.source.rxs[j].setid in rxgrp or rxgrp[0] == -1:
                        nn += 1
                        if self.source.txs[i].chan_to(self.source.rxs[j]) is not None:
                            for k in self.source.txs[i].chan_to(self.source.rxs[j]).paths.items():
                                if nff and not k[1].near_field_failed and zmax > l2db(k[1].pow) > zmin:
                                    self.xdata.append(j)
                                    self.ydata.append(k[1].delay * 1e9)
                                    self.zdata.append(l2db(k[1].pow))
                                elif not nff and zmax > l2db(k[1].pow) > zmin:
                                    self.xdata.append(j)
                                    self.ydata.append(k[1].delay * 1e9)
                                    self.zdata.append(l2db(k[1].pow))
                        else:
                            if self.ydata.__len__() > 0:
                                maxy = np.nanmax(self.ydata)
                                self.xdata.append(j)
                                self.ydata.append(maxy)
                                self.zdata.append(zmin)

                if np.max(self.ydata) == 0:
                    self.ydata[-1] = 1e-9

                if np.isnan(np.nanmin(self.zdata)):
                    for z in range(self.zdata.__len__()):
                        self.zdata[z] = zmin

                (X, Y, Z) = basint3(self.xdata, self.ydata, self.zdata, nn, ydim, zmin=zmin)
                [X, Y] = np.meshgrid(X, Y)

                if plot or topng:
                    f = mpl.figure(i)

                    if np.isnan(zmax):
                        cmin = np.min(Z)
                        cmax = np.max(Z) + np.abs(0.1 * np.max(Z))
                    else:
                        cmin = zmin
                        cmax = zmax

                    mpl.contourf(np.transpose(X), np.transpose(Y), Z, 20, cmap=cmap, vmin=cmin, vmax=cmax)

                    cb = mpl.colorbar()
                    cb.set_label('RX power, [dBm]')
                    cb.set_clim(vmin=cmin, vmax=cmax)
                    cb.set_ticks(np.linspace(start=cmin, stop=cmax, num=10, endpoint=True).tolist())
                    mpl.clim(vmin=cmin, vmax=cmax)
#                    cb.ax.set_ylim([cmin, cmax])

                    mpl.xlabel('RX Position')
                    mpl.ylabel('Delay, [ns]')
                    mpl.title('CIR@TX #{}'.format(i))
                    mpl.tight_layout()

                if topng:
                    mpl.savefig('CIR3D_tx{0:03d}_rxgrp{1:03d}.png'.format(i, rxgrp[0]))
                    mpl.close(f)

                if matsav:
                    sio.savemat('CIR3D_tx{0:03d}_rxgrp{1:03d}.mat'.format(i, rxgrp[0]), {'X': X, 'Y': Y, 'Z': Z})

        if topng is False and plot:
            mpl.show()

    def export_pdp(self, tx: list = [-1], rx: list = [-1], nff: bool = True, avg: bool = False, floor: float = -110.0,
                   matsav: bool = False, csvsav: bool = False, plot: bool = True, topng: bool = False,
                   ceil: float = -40.0, rxgrp: list = [-1], txgrp: list = [-1]):

        tx = [tx] if not isinstance(tx, list) else tx
        rx = [rx] if not isinstance(rx, list) else rx

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        if tx[0] == -1:
            txs = []
            for i in self.source.txs.keys():
                if txgrp[0] != -1 or self.source.txs[i].setid in txgrp:
                    tx.append(i)


        if rx[0] == -1:
            rxs = []
            for i in self.source.rxs.keys():
                if rxgrp[0] != -1 or self.source.rxs[i].setid in rxgrp:
                    rx.append(i)

        delay = []
        pow = []

        for i in txs:
            for j in rxs:
                if not avg:
                    delay = []
                    pow = []

                if self.source.txs[i].chan_to(self.source.rxs[j]):
                    for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.items():
                        if ceil > l2db(k[1].pow) > floor:
                            if nff and not k[1].near_field_failed:
                                delay.append(k[1].delay)
                                pow.append(l2db(k[1].pow))
                            elif not nff:
                                delay.append(k[1].delay)
                                pow.append(l2db(k[1].pow))


                else:
                    print('Error, no route between TX {} and RX {}!'.format(i, j))
                    pass

                if not avg:
                    if plot:
                        f = mpl.figure((i+1)*(j+1))
                        mpl.stem(delay, pow, bottom=-120)
                        mpl.xlabel('Delay, [s]')
                        mpl.ylabel('Power, [dBm]')
                        mpl.title('PDP@[TX{}<->RX{}]'.format(i, j))
                        offset = 0.1 * (np.nanmax(pow) - np.nanmin(pow))
                        mpl.ylim([np.nanmin(pow) - offset, np.nanmax(pow) + offset])
                        mpl.grid(linestyle='--')
                        mpl.tight_layout()

                    if matsav:
                        sio.savemat('PDP@[TX{0:02d}{1:03d}<->RX{2:02d}{3:03d}].mat'.format(i, j),
                                    {'delay': delay, 'pow': pow})

                    if csvsav:
                        file = open('PDP@[TX{0:02d}{1:03d}<->RX{2:02d}{3:03d}].csv'.format(i, j), mode='w')
                        file.write('Delay [sec],Power [dBm]\n')
                        for k in range(pow.__len__()):
                            file.write('{},{}\n'.format(delay[k], pow[k]))
                        file.close()

        if avg:
            if plot or topng:
                f = mpl.figure(0)
                mpl.stem(delay, pow, bottom=-120)
                mpl.xlabel('Delay, [s]')
                mpl.ylabel('Power, [dBm]')
                mpl.title('Average PDP@[TX<->RX]'.format(i, j))
                offset = 0.1 * (np.nanmax(pow) - np.nanmin(pow))
                mpl.ylim([np.nanmin(pow) - offset, np.nanmax(pow) + offset])
                mpl.grid(linestyle='--')
                mpl.tight_layout()

            if matsav:
                sio.savemat('PDP@[TX<->RX]_avg.mat', {'delay': delay, 'pow': pow})

            if csvsav:
                file = open('PDP@[TX<->RX]_avg.csv', mode='w')
                file.write('Delay [sec],Power [dBm]\n')
                for k in range(pow.__len__()):
                    file.write('{},{}\n'.format(delay[k], pow[k]))
                file.close()

            if topng:
                mpl.savefig('PDP@[TX<->RX]_avg.png')
                mpl.close(f)

        if topng is False and plot:
            mpl.show()

if __name__ == "__main__":
    DS = pairdata.data_stor('dbconf.txt')
    DS.load_rxtx('Human_sitting_legsback_Standing_still_ELP_sqlite')
    DS.load_paths(npaths=250)
    DS.load_interactions(store=True)
    check_data_NF(DS)
    cir = cirs(DS)
    cir.export(txgrp=-1, rxgrp=6, nff=True, matsav=True, plot=True, topng=True, zmin=-190.0, zmax=-40.0)
    cir.export(txgrp=-1, rxgrp=5, nff=True, matsav=True, plot=True, topng=True, zmin=-190.0, zmax=-40.0)
    cir.export(txgrp=-1, rxgrp=4, nff=True, matsav=True, plot=True, topng=True, zmin=-190.0, zmax=-40.0)
    cir.export(txgrp=-1, rxgrp=2, nff=True, matsav=True, plot=True, topng=True, zmin=-190.0, zmax=-40.0)
    #cir.export_pdp(csvsav=True, plot=True, avg=True, tx=0, rx=[10, 20, 30, 40])
    exit()
