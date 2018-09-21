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

__author__ = 'Aleksei Ponomarenko-Timofeev'

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as mpl
import pairdata
from phys_path_procs import *
from auxclass import PowHist

class cirs:
    def __init__(self, source):
        self.source = source

        self.xdata = []
        self.ydata = []
        self.zdata = []

        self.xmin = 0.0
        self.xmax = 0.0

        self.ymin = 0.0
        self.ymax = 0.0

        self.zmin = 0.0
        self.zmax = 0.0

        self.ydim = 0
        self.xdim = 0

        self.title = ''

        self.rxgrp = -1
        self.tx = 0

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               cmap: str = 'Blues', xdim: int = 100, ydim: int = 250, zmin: float = -200.0, zmax: float = np.nan,
               nff: bool = True, matsav: bool = False, plot: bool = True, show: bool =True, fidbase: int = 0,
               title: str = ''):

        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        self.title = title

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        for i in txrange:
            if self.source.txs[i].setid in txgrp or txgrp[0] == -1:
                self.xdata = []
                self.ydata = []
                self.zdata = []

                self.ydim = ydim
                self.xdim = 0
                for j in rxrange:
                    if self.source.rxs[j].setid in rxgrp or rxgrp[0] == -1:
                        self.xdim += 1
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
                                self.ymax = np.nanmax(self.ydata)
                                self.xdata.append(j)
                                self.ydata.append(self.ymax)
                                self.zdata.append(zmin)

                if np.max(self.ydata) == 0:
                    self.ydata[-1] = 1e-9

                if np.isnan(np.nanmin(self.zdata)):
                    for z in range(self.zdata.__len__()):
                        self.zdata[z] = zmin

                (x, y, z) = basint3(self.xdata, self.ydata, self.zdata, self.xdim, self.ydim, zmin=zmin)
                [x, y] = np.meshgrid(x, y)

                self.xmax = np.nanmax(self.xdata)
                self.xmin = np.nanmin(self.xdata)

                self.ymax = np.nanmax(self.ydata)
                self.ymin = np.nanmin(self.ydata)

                self.tx = i
                self.rxgrp = rxgrp[0]

                if np.isnan(zmax):
                    self.zmin = np.nanmin(z)
                    self.zmax = np.nanmax(z) + np.abs(0.1 * np.nanmax(z))
                else:
                    self.zmin = zmin
                    self.zmax = zmax

                if plot or mkpng:
                    f = mpl.figure(fidbase + i)

                    mpl.pcolor(np.transpose(x), np.transpose(y), z, cmap=cmap, vmin=self.zmin, vmax=self.zmax)

                    cb = mpl.colorbar(ticks=np.linspace(start=self.zmin, stop=self.zmax, num=11, endpoint=True).tolist())
                    cb.set_label('RX power, [dBm]')
                    cb.set_clim(vmin=self.zmin, vmax=self.zmax)
                    mpl.clim(vmin=self.zmin, vmax=self.zmax)

                    mpl.xlabel('RX Position')
                    mpl.ylabel('Delay, [ns]')
                    mpl.title('{}CIR\@TX \#{}'.format(title, i))
                    mpl.tight_layout()

                if mkpng:
                    mpl.savefig('{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.png'.format(i, rxgrp[0], title))
                    mpl.close(f)

                if matsav:
                    sio.savemat('{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.mat'.format(i, rxgrp[0], title), {'X': x, 'Y': y, 'Z': z})

                if not mkpng and (not plot or not show):
                    mpl.close(f)

        if mkpng is False and plot and show:
            mpl.show()

    def export_pdp(self, tx: list = [-1], rx: list = [-1], nff: bool = True, avg: bool = False, floor: float = -110.0,
                   matsav: bool = False, csvsav: bool = False, plot: bool = True, mkpng: bool = False,
                   ceil: float = -40.0, rxgrp: list = [-1], txgrp: list = [-1], title: str = ''):

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
                        mpl.title('{}PDP@[TX{}<->RX{}]'.format(title, i, j))
                        offset = 0.1 * (np.nanmax(pow) - np.nanmin(pow))
                        mpl.ylim([np.nanmin(pow) - offset, np.nanmax(pow) + offset])
                        mpl.grid(linestyle='--')
                        mpl.tight_layout()

                    if matsav:
                        sio.savemat('{4}PDP@[TX{0:02d}{1:03d}<->RX{2:02d}{3:03d}].mat'.format(i, j, title),
                                    {'delay': delay, 'pow': pow})

                    if csvsav:
                        file = open('{4}PDP@[TX{0:02d}{1:03d}<->RX{2:02d}{3:03d}].csv'.format(i, j, title), mode='w')
                        file.write('Delay [sec],Power [dBm]\n')
                        for k in range(pow.__len__()):
                            file.write('{},{}\n'.format(delay[k], pow[k]))
                        file.close()

        if avg:
            if plot or mkpng:
                f = mpl.figure(0)
                mpl.stem(delay, pow, bottom=-120)
                mpl.xlabel('Delay, [s]')
                mpl.ylabel('Power, [dBm]')
                mpl.title('Average {}PDP\@[TX<->RX]'.format(title))
                offset = 0.1 * (np.nanmax(pow) - np.nanmin(pow))
                mpl.ylim([np.nanmin(pow) - offset, np.nanmax(pow) + offset])
                mpl.grid(linestyle='--')
                mpl.tight_layout()

            if matsav:
                sio.savemat('{}PDP\@[TX<->RX]_avg.mat'.format(title), {'delay': delay, 'pow': pow})

            if csvsav:
                file = open('{}PDP\@[TX<->RX]_avg.csv'.format(title), mode='w')
                file.write('Delay [sec],Power [dBm]\n')
                for k in range(pow.__len__()):
                    file.write('{},{}\n'.format(delay[k], pow[k]))
                file.close()

            if mkpng:
                mpl.savefig('{}PDP\@[TX<->RX]_avg.png'.format(title))
                mpl.close(f)

        if mkpng is False and plot:
            mpl.show()

    def delay_spread_export(self, txgrp: list = [-1], rxgrp: list = [-1], mkpng: bool = False,
               cmap: str = 'Blues', floor: float = -200.0, ceil: float = np.nan,
               nff: bool = True, matsav: bool = False, plot: bool = True, show: bool =True, fidbase: int = 0,
               title: str = ''):

        txs = list()
        rxs = list()

        for i in self.source.txs.keys():
            if txgrp[0] != -1 or self.source.txs[i].setid in txgrp:
                txs.append(i)

        for i in self.source.rxs.keys():
            if rxgrp[0] != -1 or self.source.rxs[i].setid in rxgrp:
                rxs.append(i)

        for i in txs:
            ds_avg = list()
            ds_avg_rms = list()

            for j in rxs:
                delay = list()
                power = list()

                if self.source.txs[i].chan_to(self.source.rxs[j]):
                    for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.items():
                        if ceil > l2db(k[1].pow) > floor:
                            if nff and not k[1].near_field_failed:
                                delay.append(k[1].delay)
                                power.append(l2db(k[1].pow))
                            elif not nff:
                                delay.append(k[1].delay)
                                power.append(l2db(k[1].pow))
                else:
                    print('Error, no route between TX {} and RX {}!'.format(i, j))
                    pass

        #         ds_
        #
        #         if not avg:
        #             if plot:
        #                 f = mpl.figure((i+1)*(j+1))
        #                 mpl.stem(delay, power, bottom=-120)
        #                 mpl.xlabel('Delay, [s]')
        #                 mpl.ylabel('Power, [dBm]')
        #                 mpl.title('{}PDP@[TX{}<->RX{}]'.format(title, i, j))
        #                 offset = 0.1 * (np.nanmax(power) - np.nanmin(power))
        #                 mpl.ylim([np.nanmin(power) - offset, np.nanmax(power) + offset])
        #                 mpl.grid(linestyle='--')
        #                 mpl.tight_layout()
        #
        #             if matsav:
        #                 sio.savemat('{4}PDP@[TX{0:02d}{1:03d}<->RX{2:02d}{3:03d}].mat'.format(i, j, title),
        #                             {'delay': delay, 'pow': power})
        #
        #             if csvsav:
        #                 file = open('{4}PDP@[TX{0:02d}{1:03d}<->RX{2:02d}{3:03d}].csv'.format(i, j, title), mode='w')
        #                 file.write('Delay [sec],Power [dBm]\n')
        #                 for k in range(power.__len__()):
        #                     file.write('{},{}\n'.format(delay[k], power[k]))
        #                 file.close()
        #
        # if avg:
        #     if plot or mkpng:
        #         f = mpl.figure(0)
        #         mpl.stem(delay, power, bottom=-120)
        #         mpl.xlabel('Delay, [s]')
        #         mpl.ylabel('Power, [dBm]')
        #         mpl.title('Average {}PDP\@[TX<->RX]'.format(title))
        #         offset = 0.1 * (np.nanmax(power) - np.nanmin(power))
        #         mpl.ylim([np.nanmin(power) - offset, np.nanmax(power) + offset])
        #         mpl.grid(linestyle='--')
        #         mpl.tight_layout()
        #
        #     if matsav:
        #         sio.savemat('{}PDP\@[TX<->RX]_avg.mat'.format(title), {'delay': delay, 'pow': power})
        #
        #     if csvsav:
        #         file = open('{}PDP\@[TX<->RX]_avg.csv'.format(title), mode='w')
        #         file.write('Delay [sec],Power [dBm]\n')
        #         for k in range(power.__len__()):
        #             file.write('{},{}\n'.format(delay[k], power[k]))
        #         file.close()
        #
        #     if mkpng:
        #         mpl.savefig('{}PDP\@[TX<->RX]_avg.png'.format(title))
        #         mpl.close(f)
        #
        # if mkpng is False and plot:
        #     mpl.show()

if __name__ == "__main__":
    DS = pairdata.DataStorage('dbconf.txt')
    DS.load_rxtx('Human_sitting_legsback_Sitting_fleece_sqlite')
    DS.load_paths(npaths=250)
    DS.load_interactions(store=True)
    #check_data_NF(DS)
    cir = cirs(DS)
    cir.export(txgrp=-1, rxgrp=6, nff=True, matsav=False, plot=True, mkpng=False, zmin=-190.0, zmax=-40.0)
    cir.export(txgrp=-1, rxgrp=5, nff=True, matsav=False, plot=True, mkpng=False, zmin=-190.0, zmax=-40.0)
    cir.export(txgrp=-1, rxgrp=4, nff=True, matsav=False, plot=True, mkpng=False, zmin=-190.0, zmax=-40.0)
    cir.export(txgrp=-1, rxgrp=2, nff=True, matsav=False, plot=True, mkpng=False, zmin=-190.0, zmax=-40.0)
    #cir.export_pdp(csvsav=True, plot=True, avg=True, tx=0, rx=[10, 20, 30, 40])
    exit()
