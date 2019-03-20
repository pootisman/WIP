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

import matplotlib.pyplot as mpl
#import mayavi.mlab as mlab
import scipy.io as sio
import numpy as np
from pairdata import DataStorage
from auxfun import l2db, basint3, enable_latex
from auxclass import PowHist, PowHist2D

class RXPatAz:
    def __init__(self, source: DataStorage):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkimg: str='',
               nff: bool = True, csvsav: bool = False, matsav: bool = False, fnappend: str = '',
               rlims: tuple = (-110, -75), drawtit: bool = False, gennpz: bool=False):
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
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].aoa,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))
                            elif not nff:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].aoa,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))

                        th = sorted(th, key=lambda x: x[0])

                        hist = PowHist(rstart=-180.0, rstop=180.0, binc=72)
                        for t in th:
                            hist.append(t[0], t[1])

                        for t in hist.bins.items():
                            tt.append(t[0][0])
                            tt.append(t[0][1])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])

                        tt.append(tt[0])
                        r.append(r[0])

                        ax = mpl.subplot(111, projection='polar')
                        ax.fill(np.deg2rad(tt), r)
                        ax.set_rmax(rlims[1])
                        ax.set_rmin(rlims[0])
                        ax.grid(linestyle='--')
                        if drawtit:
                            mpl.title('aoa@[TX #{} - RX #{}]'.format(i, j))
                        
                        if mkimg != '':
                            mpl.savefig('{2}RXaz.tx{0:03d}-rx{1:03d}.{3}'.format(i, j, fnappend, mkimg))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('{2}RXaz.tx{0:03d}-rx{1:03d}.mat'.format(i, j, fnappend),
                                        {'theta': tt, 'pow': r})

                        if gennpz:
                            np.savez_compressed('{2}RXaz.tx{0:03d}-rx{1:03d}'.format(i, j, fnappend),
                                                theta=tt, power=r)

                        if csvsav:
                            file = open('{2}RXaz.tx{0:03d}-rx{1:03d}.csv'.format(i, j, fnappend), mode='w')
                            file.write('Ang. [deg], Pow [dBm]\n')
                            for k in range(tt.__len__()):
                                file.write('{},{}\n'.format(tt[k], r[k]))
                            file.close()

        if mkimg == '':
            mpl.show()


class RXPatEl:
    def __init__(self, source: DataStorage):
        self.source = source

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkimg: str = '',
               nff: bool = True, csvsav: bool = False, matsav: bool = False, fnappend: str = '',
               rlims: tuple = (-110, -75), drawtit: bool = False, gennpz: bool=False):
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
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].eoa,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))
                            elif not nff:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].eoa,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))

                        th = sorted(th, key=lambda x: x[0])

                        hist = PowHist(rstart=-180.0, rstop=180.0, binc=72)
                        for t in th:
                            hist.append(t[0], t[1])

                        for t in hist.bins.items():
                            tt.append(t[0][0])
                            tt.append(t[0][1])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])

                        tt.append(tt[0])
                        r.append(r[0])

                        ax = mpl.subplot(111, projection='polar')
                        ax.fill(np.deg2rad(tt), r)
                        ax.set_rmax(rlims[1])
                        ax.set_rmin(rlims[0])
                        ax.grid(linestyle='--')
                        if drawtit:
                            mpl.title('eoa@[TX #{} - RX #{}]'.format(i, j))

                        if mkimg != '':
                            mpl.savefig('{2}RXel.tx{0:03d}-rx{1:03d}.{3}'.format(i, j, fnappend, mkimg))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('{2}RXel.tx{0:03d}-rx{1:03d}.mat'.format(i, j, fnappend), {'theta': tt, 'pow': r})

                        if gennpz:
                            np.savez_compressed('{2}RXel_tx{0:03d}-rx{1:03d}'.format(i, j, fnappend),
                                                theta=tt, power=r)

                        if csvsav:
                            file = open('{2}RXel.tx{0:03d}-rx{1:03d}.csv'.format(i, j, fnappend), mode='w')
                            file.write('Ang. [deg], Pow [dBm]\n')
                            for k in range(tt.__len__()):
                                file.write('{},{}\n'.format(tt[k], r[k]))
                            file.close()

        if mkimg == '':
            mpl.show()


class TXPatAz:
    def __init__(self, source: DataStorage):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkimg: str = '',
               nff: bool = True, csvsav: bool = False, matsav: bool = False, fnappend: str = '',
               rlims: tuple = (-110, -75), drawtit: bool = False, gennpz: bool=False):
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
                            if nff and not self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k]. \
                                    near_field_failed:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].aod,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))
                            elif not nff:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].aod,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))

                        th = sorted(th, key=lambda x: x[0])

                        hist = PowHist(rstart=-180.0, rstop=180.0, binc=72)
                        for t in th:
                            hist.append(t[0], t[1])

                        for t in hist.bins.items():
                            tt.append(t[0][0])
                            tt.append(t[0][1])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])

                        tt.append(tt[0])
                        r.append(r[0])

                        ax = mpl.subplot(111, projection='polar')
                        ax.fill(np.deg2rad(tt), r)
                        ax.set_rmax(rlims[1])
                        ax.set_rmin(rlims[0])
                        ax.grid(linestyle='--')
                        if drawtit:
                            mpl.title('aoa@[TX #{} - RX #{}]'.format(i, j))

                        if mkimg:
                            mpl.savefig('{2}TXaz.tx{0:03d}-rx{1:03d}.{3}'.format(i, j, fnappend, mkimg))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('{2}TXaz.tx{0:03d}-rx{1:03d}.mat'.format(i, j, fnappend),
                                        {'theta': tt, 'pow': r})

                        if gennpz:
                            np.savez_compressed('{2}TXel.tx{0:03d}-rx{1:03d}'.format(i, j, fnappend),
                                                theta=tt, power=r)

                        if csvsav:
                            file = open('{2}TXaz_tx{0:03d}-rx{1:03d}.csv'.format(i, j, fnappend), mode='w')
                            file.write('Ang. [deg], Pow [dBm]\n')
                            for k in range(tt.__len__()):
                                file.write('{},{}\n'.format(tt[k], r[k]))
                            file.close()

        if mkimg == '':
            mpl.show()


class TXPatEl:
    def __init__(self, source: DataStorage):
        self.source = source

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkimg: str = '',
               nff: bool = True, csvsav: bool = False, matsav: bool = False, fnappend: str = '',
               rlims: tuple = (-110, -75), drawtit: bool = False, npzgen: bool = False):
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
                            if nff and not self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k]. \
                                    near_field_failed:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].eod,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))
                            elif not nff:
                                th.append((self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].eod,
                                           self.source.txs[i].chans_to_pairs[self.source.rxs[j]].paths[k].pow))

                        th = sorted(th, key=lambda x: x[0])

                        hist = PowHist(rstart=-180.0, rstop=180.0, binc=72)
                        for t in th:
                            hist.append(t[0], t[1])

                        for t in hist.bins.items():
                            tt.append(t[0][0])
                            tt.append(t[0][1])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])
                            r.append(l2db(t[1]) if t[1] > 0 and l2db(t[1]) > rlims[0] else rlims[0])

                        tt.append(tt[0])
                        r.append(r[0])

                        ax = mpl.subplot(111, projection='polar')
                        ax.fill(np.deg2rad(tt), r)
                        ax.set_rmax(rlims[1])
                        ax.set_rmin(rlims[0])
                        ax.grid(linestyle='--')

                        if drawtit:
                            mpl.title('eoa@[TX #{} - RX #{}]'.format(i, j))

                        if mkimg != '':
                            mpl.savefig('{2}TXel.tx{0:03d}-rx{1:03d}.{3}'.format(i, j, fnappend, mkimg))
                            mpl.close(f)

                        if matsav:
                            sio.savemat('{2}TXel.tx{0:03d}-rx{1:03d}.mat'.format(i, j, fnappend),
                                        {'theta': tt, 'pow': r})

                        if csvsav:
                            file = open('{2}TXel.tx{0:03d}-rx{1:03d}.csv'.format(i, j, fnappend), mode='w')
                            file.write('Ang. [deg], Pow [dBm]\n')
                            for k in range(tt.__len__()):
                                file.write('{},{}\n'.format(tt[k], r[k]))
                            file.close()

        if mkimg == '':
            mpl.show()


class RXPatAll:
    def __init__(self, source: DataStorage):
        self.source = source

    def export(self, txrange: int = -1, rxrange: int = -1, txgrp: int = -1, rxgrp: int = -1, mkpng: bool = False,
               nff: bool = True, matsav: bool = False, mkpdf: bool = False):
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
                                hist.append(k[1].aod, k[1].eod, k[1].pow)

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
                            mlab.savefig('RXpat.tx{0:03d}-rx{1:03d}.png'.format(i, j))
                            mlab.close(f)

                        if mkpdf:
                            mlab.savefig('RXpat.tx{0:03d}-rx{1:03d}.pdf'.format(i, j))
                            mlab.close(f)

                        if matsav:
                            sio.savemat('RXPat.tx{0:03d}-rx{1:03d}.mat'.format(i, j), {'X': x, 'Y': y, 'Z': z})

        if mkpng is False:
            mlab.show()


if __name__ == "__main__":
    enable_latex(22)
    DS = DataStorage('dbconf.txt')
    DS.load_rxtx('BUSMOD')
    DS.load_paths(npaths=250)
    DS.load_interactions()
    cir = TXPatAz(DS)
    cir.export(rxgrp=4, mkimg='', gennpz=True, rlims=(-100, -75), drawtit=False)
