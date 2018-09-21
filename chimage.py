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

import matplotlib.pyplot as mpl
import scipy.io as sio
import pairdata
from auxclass import *

__author__ = 'Aleksei Ponomarenko-Timofeev'


class chimage_RX():
    def __init__(self, source):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               cmap: str = 'viridis', nff: bool = True, matsav: bool = False, zmin: float = np.nan,
               zmax: float = np.nan):
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
                rr = 0
                for j in rxrange:
                    if self.source.rxs[j].setid in rxgrp or rxgrp[0] == -1:
                        f = mpl.figure(rr)
                        rr += 1
                        hist = PowHist2D(azbinc=180, elbinc=90)

                        for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.items():
                            if nff and not k[1].near_field_failed:
                                hist.append(k[1].AoA, k[1].EoA, k[1].pow)

                        X = []
                        Y = []
                        Z = []

                        for k in hist.bins.items():
                            X.append(k[0][0])
                            X.append(k[0][1])
                            Y.append(k[0][2])
                            Y.append(k[0][3])
                            Z.append(k[1] if k[1] > 0 else np.nan)
                            Z.append(Z[-1])

                        if np.isnan(zmin):
                            zmin = np.nanmin(Z) if not np.isnan(np.nanmin(Z)) else np.finfo(float).eps

                        if np.isnan(zmax):
                            zmax = l2db(np.nanmax(Z))

                        for k in range(Z.__len__()):
                            Z[k] = l2db(Z[k]) if not np.isnan(Z[k]) else zmin

                        X, Y, Z = basint3(X, Y, Z, xc=hist.azbinc, yc=hist.elbinc)

                        #mpl.contourf(X, Y, Z.T, 20, cmap=cmap, vmin=zmin, vmax=zmax)
                        mpl.pcolor(X, Y, Z.T, cmap=cmap, vmin=zmin, vmax=zmax)
                        mpl.grid(linestyle='--')
                        mpl.title('Channel RX Image\@[TX \#{} $\\rightarrow$ RX \#{}]'.format(i, j))
                        mpl.xlabel('Azimuth, [degrees]')
                        mpl.ylabel('Elevation, [degrees]')
                        cbr = mpl.colorbar()
                        cbr.set_label('RX Power, [dBm]')
                        #cbr.ax.set_xlim(zmin, zmax)
                        mpl.tight_layout()
                        if mkpng:
                            mpl.savefig('RXCImage_tx[{0:01d}.{1:03d}]->rx[{2:01d}.{3:03d}].png'.
                                        format(self.source.txs[i].setid, i, self.source.rxs[j].setid, j))
                            mpl.close(f)
                        if matsav:
                            sio.savemat('RXCImage_tx[{0:01d}.{1:03d}]->rx[{2:01d}.{3:03d}].mat'.
                                        format(self.source.txs[i].setid, i, self.source.rxs[j].setid, j),
                                        {'X': X, 'Y': Y, 'Z': Z})

        if mkpng is False:
            mpl.show()


class chimage_TX():
    def __init__(self, source):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               cmap: str = 'viridis', nff: bool = True, matsav: bool = False, show: bool = True, zmin: float = np.nan,
               zmax: float = np.nan):
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
                rr = 0
                for j in rxrange:
                    if self.source.rxs[j].setid in rxgrp or rxgrp[0] == -1:
                        if mkpng or show:
                            f = mpl.figure(rr)

                        rr += 1
                        hist = PowHist2D(azbinc=180, elbinc=90)

                        for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.items():
                            if nff and not k[1].near_field_failed:
                                hist.append(k[1].AoD, k[1].EoD, k[1].pow)

                        X = []
                        Y = []
                        Z = []

                        for k in hist.bins.items():
                            X.append(k[0][0])
                            X.append(k[0][1])
                            Y.append(k[0][2])
                            Y.append(k[0][3])
                            Z.append(k[1] if k[1] > 0 else np.nan)
                            Z.append(Z[-1])

                        if np.isnan(zmin):
                            zmin = np.nanmin(Z) if not np.isnan(np.nanmin(Z)) else np.finfo(float).eps

                        if np.isnan(zmax):
                            zmax = l2db(np.nanmax(Z))

                        for k in range(Z.__len__()):
                            Z[k] = l2db(Z[k]) if not np.isnan(Z[k]) else zmin

                        X, Y, Z = basint3(X, Y, Z, xc=hist.azbinc, yc=hist.elbinc)

                        if mkpng or show:
                            #mpl.contourf(X, Y, Z.T, 20, cmap=cmap, vmin=zmin, vmax=zmax)
                            mpl.pcolor(X, Y, Z.T,  cmap=cmap, vmin=zmin, vmax=zmax)
                            mpl.grid(linestyle='--')
                            mpl.title('Channel TX Image\@[TX \#{} $\\rightarrow$ RX \#{}]'.format(i, j))
                            mpl.xlabel('Azimuth, [degrees]')
                            mpl.ylabel('Elevation, [degrees]')
                            cbr = mpl.colorbar()
                            cbr.set_label('RX Power, [dBm]')
                            mpl.tight_layout()

                        if mkpng or not show:
                            mpl.savefig('TXCImage_tx[{0:01d}.{1:03d}]->rx[{2:01d}.{3:03d}].png'.format(self.source.txs[i].setid, i,
                                                                                       self.source.rxs[j].setid, j))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('RXCImage_tx[{0:01d}.{1:03d}]->rx[{2:01d}.{3:03d}].mat'.format(self.source.txs[i].setid, i,
                                                                                       self.source.rxs[j].setid, j),
                            {'X': X, 'Y': Y, 'Z': Z})

        if mkpng is False and show:
            mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor(conf='dbconf.txt')
    DS.load_rxtx(dbname='Human_sitting_legsback_Sitting_sqlite')
    DS.load_paths(npaths=250)
    DS.load_interactions(store=True)

    from phys_path_procs import *
    check_data_NF(DS)

    enable_latex()

    DE = chimage_RX(DS)
    DE.export(rxgrp=[2,6,5], mkpng=True, nff=True, zmin=-130.0, zmax=-40.0)
    DT = chimage_TX(DS)
    DT.export(rxgrp=[2,6,5], mkpng=True, nff=True, zmin=-130.0, zmax=-40.0)
    exit()