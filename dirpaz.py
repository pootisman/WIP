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
import mayavi.mlab as mlab
import scipy.io as sio
import numpy as np
from pairdata import DataStorage
from auxfun import l2db, basint3
from auxclass import PowHist, PowHist2D

__author__ = 'Aleksei Ponomarenko-Timofeev'


class RXPatAz:
    def __init__(self, source: DataStorage):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               nff: bool = True, csvsav: bool = False, matsav: bool = False):
        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        for i in txrange:
            if txgrp in (self.source.txs[i].setid, -1):
                rr = 0
                for j in rxrange:
                    if rxgrp in (self.source.rxs[j].setid, -1):
                        f = mpl.figure(rr)
                        rr += 1
                        th = list()
                        tt = list()
                        r = list()
                        for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.keys():
                            if nff and not self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].\
                                           near_field_failed:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].AoA,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))
                            elif not nff:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].AoA,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))

                        th = sorted(th, key=lambda x: x[0])

                        hist = PowHist(rstart=-180.0, rstop=180.0, binc=36)
                        for t in th:
                            hist.append(t[0], t[1])

                        for t in hist.bins.items():
                            tt.append(t[0][0])
                            tt.append(t[0][1])
                            r.append(l2db(t[1]) if t[1] > 0 else -180)
                            r.append(l2db(t[1]) if t[1] > 0 else -180)

                        tt.append(tt[0])
                        r.append(r[0])

                        ax = mpl.subplot(111, projection='polar')
                        ax.fill(np.deg2rad(tt), r)
                        ax.set_rmax(-60)
                        ax.set_rmin(-180)
                        ax.grid(linestyle='--')
                        mpl.title('aoa@[TX #{} -> RX #{}]'.format(i, j))
                        if mkpng:
                            mpl.savefig('RXaz_tx{0:03d}->rx{1:03d}.png'.format(i, j))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('RXaz_tx{0:03d}->rx{1:03d}.mat'.format(i, j), {'theta': tt, 'pow': r})

                        if csvsav:
                            file = open('RXaz_tx{0:03d}->rx{1:03d}.csv'.format(i, j), mode='w')
                            file.write('Ang. [deg], Pow [dBm]\n')
                            for k in range(tt.__len__()):
                                file.write('{},{}\n'.format(tt[k], r[k]))
                            file.close()

        if mkpng is False:
            mpl.show()


class RXPatEl:
    def __init__(self, source: DataStorage):
        self.source = source

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               nff: bool = True, csvsav: bool = False, matsav: bool = False):
        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        for i in txrange:
            if txgrp in (self.source.txs[i].setid, -1):
                rr = 0
                for j in rxrange:
                    if rxgrp in (self.source.rxs[j].setid, -1):
                        f = mpl.figure(rr)
                        rr += 1
                        th = list()
                        tt = list()
                        r = list()
                        for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.keys():
                            if nff and not self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].\
                                    near_field_failed:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].EoA,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))
                            elif not nff:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].AoA,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))

                        th = sorted(th, key=lambda x: x[0])

                        hist = PowHist(rstart=-180.0, rstop=180.0, binc=36)
                        for t in th:
                            hist.append(t[0], t[1])

                        for t in hist.bins.items():
                            tt.append(t[0][0])
                            tt.append(t[0][1])
                            r.append(l2db(t[1]) if t[1] > 0 else -180)
                            r.append(l2db(t[1]) if t[1] > 0 else -180)

                        tt.append(tt[0])
                        r.append(r[0])

                        ax = mpl.subplot(111, projection='polar')
                        ax.fill(np.deg2rad(tt), r)
                        ax.set_rmax(-60)
                        ax.set_rmin(-180)
                        ax.grid(linestyle='--')
                        mpl.title('aoa@[TX #{} -> RX #{}]'.format(i, j))
                        if mkpng:
                            mpl.savefig('RXel_tx{0:03d}->rx{1:03d}.png'.format(i, j))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('RXel_tx{0:03d}->rx{1:03d}.mat'.format(i, j), {'theta': tt, 'pow': r})

                        if csvsav:
                            file = open('RXel_tx{0:03d}->rx{1:03d}.csv'.format(i, j), mode='w')
                            file.write('Ang. [deg], Pow [dBm]\n')
                            for k in range(tt.__len__()):
                                file.write('{},{}\n'.format(tt[k], r[k]))
                            file.close()

        if mkpng is False:
            mpl.show()


class RXPatAll:
    def __init__(self, source: DataStorage):
        self.source = source

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               nff: bool = True, matsav: bool = False):
        if txrange == -1:
            txrange = self.source.txs.keys()
        else:
            txrange = range(txrange)

        if rxrange == -1:
            rxrange = self.source.rxs.keys()
        else:
            rxrange = range(rxrange)

        for i in txrange:
            if txgrp in (self.source.txs[i].setid, -1):
                rr = 0
                for j in rxrange:
                    if rxgrp in (self.source.rxs[j].setid, -1):
                        f = mlab.figure(j)
                        rr += 1
                        hist = PowHist2D()

                        for k in self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths.items():
                            if nff and not k[1].near_field_failed:
                                hist.append(k[1].AoD, k[1].EoD, k[1].pow)

                        x = list()
                        y = list()
                        z = list()

                        for k in hist.bins.items():
                            x.append(k[0][0])
                            x.append(k[0][1])
                            y.append(k[0][2])
                            y.append(k[0][3])
                            z.append(k[1] if k[1] > 0 else np.nan)
                            z.append(z[-1])

                        minz = float(np.nanmin(z))

                        for k in range(z.__len__()):
                            z[k] = l2db(z[k]) if not np.isnan(z[k]) else l2db(minz)

                        x, y, z = basint3(x, y, z, xc=hist.azbinc, yc=hist.elbinc)

                        xlen = x.size
                        ylen = y.size

                        plot_z = z - np.min(np.min(z))

                        mlab.mesh((plot_z.T * np.cos(np.tile(np.deg2rad(x), [ylen, 1])) *
                                   np.sin(np.tile(np.deg2rad(y), [xlen, 1])).T),
                                  (plot_z.T * np.sin(np.tile(np.deg2rad(x), [ylen, 1]))
                                   * np.sin(np.tile(np.deg2rad(y), [xlen, 1])).T),
                                  (plot_z.T * np.cos(np.tile(np.deg2rad(x), [ylen, 1]))))

                        if mkpng:
                            mlab.savefig('RXpat_tx{0:03d}->rx{1:03d}.png'.format(i, j))
                            mlab.close(f)

                        if matsav:
                            sio.savemat('RXPat_tx{0:03d}->rx{1:03d}.mat'.format(i, j), {'X': x, 'Y': y, 'Z': z})

        if mkpng is False:
            mlab.show()


if __name__ == "__main__":
    DS = DataStorage()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/'
                 'Class@60GHz.TEST_60_MKE_15.sqlite')
    DS.load_paths(npaths=250)
    DS.load_interactions()
    cir = RXPatAll(DS)
    cir.export(rxgrp=4, mkpng=False)
