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
import matplotlib.gridspec as gsp
from pairdata import DataStorage
from phys_path_procs import *


class CIR:
    def __init__(self, source: DataStorage):
        self.source = source

        self.xdata = list()
        self.ydata = list()
        self.zdata = list()

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

    def export(self, txrange: list = [-1], rxrange: list = [-1], txgrp: list = [-1], rxgrp: list = [-1], mkpng: bool = False,
               cmap: str = 'Blues', xdim: int = 100, ydim: int = 250, zmin: float = -200.0, zmax: float = np.nan,
               nff: bool = True, matsav: bool = False, plot: bool = True, show: bool =True, fidbase: int = 0,
               title: str = '', iconvec: dict = None, disptitle: bool = True, dispylabel: bool = True, dispxlabel: bool
               = True, yrange: tuple = (np.NINF, np.PINF), mkpdf: bool = False, fnappend: str = ''):

        accept = lambda o: ((nff and not o.near_field_failed) or not nff) and (zmin < l2db(o.pow) < zmax) and\
                 (yrange[0] < o.delay * 1e9 < yrange[1])

        self.title = title

        for i in self.source.txs.items(txrange=txrange, grprange=txgrp):
            self.xdata = list()
            self.ydata = list()
            self.zdata = list()

            self.ydim = ydim
            self.xdim = 0
            for j in self.source.rxs.items(rxrange=rxrange, grprange=rxgrp):
                self.xdim += 1
                if i[1].chan_to(j[1]) is not None:
                    for k in i[1].chan_to(j[1]).paths.items():
                        if accept(k[1]):
                            self.xdata.append(j[1].node_id)
                            self.ydata.append(k[1].delay * 1e9)
                            self.zdata.append(l2db(k[1].pow))
                else:
                    if self.ydata.__len__() > 0:
                        self.ymax = np.nanmax(self.ydata)
                        self.xdata.append(j[1].node_id)
                        self.ydata.append(self.ymax)
                        self.zdata.append(zmin)

                if yrange[1] != np.PINF:
                    self.xdata.append(j[1].node_id)
                    self.ydata.append(yrange[1])
                    self.zdata.append(zmin)

                if yrange[0] != np.NINF:
                    self.xdata.append(j[1].node_id)
                    self.ydata.append(yrange[0])
                    self.zdata.append(zmin)

            if np.nanmax(self.ydata) == 0:
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
                if iconvec is not None:
                    f = mpl.figure(fidbase + i[1].node_id, constrained_layout=True)
                else:
                    f = mpl.figure(fidbase + i[1].node_id)

                if iconvec is not None:
                    gs_top = gsp.GridSpec(2, 1, figure=f, height_ratios=[5, 1])
                    ax = f.add_subplot(gs_top[0, 0])
                else:
                    ax = mpl.gca()

                mpl.sca(ax)

                mpl.pcolor(np.transpose(x), np.transpose(y), z, cmap=cmap, vmin=self.zmin, vmax=self.zmax)
                cb = mpl.colorbar(ticks=np.linspace(start=self.zmin, stop=self.zmax, num=5, endpoint=True).tolist())
                cb.set_label('RX power, [dBm]')
                cb.set_clim(vmin=self.zmin, vmax=self.zmax)
                mpl.clim(vmin=self.zmin, vmax=self.zmax)
                mpl.yticks(np.linspace(start=self.ymin, stop=self.ymax, num=5, endpoint=True).tolist())

                if dispylabel:
                    mpl.ylabel('Delay, [ns]')
                if disptitle:
                    mpl.title('{}CIR\@[TX{} $\\rightarrow$ RXg{}]'.format(title, i[1].node_id, rxgrp))

                if iconvec is None:
                    if dispxlabel:
                        mpl.xlabel('RX Position')
                    mpl.tight_layout()
                else:
                    ax, sf = f.add_subplot(gs_top[1, 0])
                    mpl.sca(ax)
                    for i in iconvec:
                        sax = sf.add_axes(0.1, 0.1, 0.1, 0.1)
                        sax.set_frame_on(True)

            if mkpng:
                mpl.savefig('{3}{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.png'.format(i[1].node_id, rxgrp[0], title, fnappend))

            if mkpdf:
                mpl.savefig('{3}{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.pdf'.format(i[1].node_id, rxgrp[0], title, fnappend))

            if matsav:
                sio.savemat('{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.mat'.format(i[1].node_id, rxgrp[0], title), {'X': x, 'Y': y, 'Z': z})

            if not plot or not show:
                mpl.close(f)

        if mkpng is False and plot and show:
            mpl.show()

    def export_single(self, txrange: list = [-1], rxrange: list = [-1], nff: bool = True, avg: bool = False, floor: float = -110.0,
                   matsav: bool = False, csvsav: bool = False, plot: bool = True, mkpng: bool = False,
                   ceil: float = -40.0, rxgrp: list = [-1], txgrp: list = [-1], title: str = '', mksvg: bool = False):

        delay = []
        pow = []

        for i in self.source.txs.items(grprange=txgrp, txrange=txrange):
            for j in self.source.rxs.items(grprange=rxgrp, rxrange=rxrange):
                if not avg:
                    delay = []
                    pow = []

                if i[1].chan_to(j[1]):
                    for k in i[1].chan_to(j[1]).paths.items():
                        if ceil > l2db(k[1].pow) > floor:
                            if nff and not k[1].near_field_failed:
                                delay.append(k[1].delay)
                                pow.append(l2db(k[1].pow))
                            elif not nff:
                                delay.append(k[1].delay)
                                pow.append(l2db(k[1].pow))
                else:
                    print('Error, no route between TX {} and RX {}!'.format(i[1].node_id, j[1].node_id))
                    pass

                if not avg:
                    if plot:
                        f = mpl.figure((i[1].node_id+1)*(j[1].node_id+1))
                        mpl.plot(delay, pow, '.')
                        mpl.xlabel('Delay, [s]')
                        mpl.ylabel('Power, [dBm]')
                        mpl.title('{}CIR@[TX{} $\\rightarrow$ RX{}]'.format(title, i[1].node_id, j[1].node_id))
                        offset = 0.1 * (np.nanmax(pow) - np.nanmin(pow))
                        mpl.ylim([np.nanmin(pow) - offset, np.nanmax(pow) + offset])
                        mpl.grid(linestyle='--')
                        mpl.tight_layout()
                        if mksvg:
                            mpl.savefig('{}CIR@[TX{}-RX{}].svg'.format(title, i[1].node_id, j[1].node_id))

                        if mkpng:
                            mpl.savefig('{}CIR@[TX{}-RX{}].png'.format(title, i[1].node_id, j[1].node_id))

                    if matsav:
                        sio.savemat('{2}CIR@[TX{0:03d}<->RX{1:03d}].mat'.format(i[1].node_id, j[1].node_id, title),
                                    {'delay': delay, 'pow': pow})

                    if csvsav:
                        file = open('{2}CIR@[TX{0:03d}<->RX{1:03d}].csv'.format(i[1].node_id, j[1].node_id, title), mode='w')
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
                mpl.title('Average {}CIR\@[TX $\\rightarrow$ RX]'.format(title))
                offset = 0.1 * (np.nanmax(pow) - np.nanmin(pow))
                mpl.ylim([np.nanmin(pow) - offset, np.nanmax(pow) + offset])
                mpl.grid(linestyle='--')
                mpl.tight_layout()

            if matsav:
                sio.savemat('{}CIR\@[TX<->RX]_avg.mat'.format(title), {'delay': delay, 'pow': pow})

            if csvsav:
                file = open('{}CIR\@[TX<->RX]_avg.csv'.format(title), mode='w')
                file.write('Delay [sec],Power [dBm]\n')

                for k in range(pow.__len__()):
                    file.write('{},{}\n'.format(delay[k], pow[k]))

                file.close()

        if mksvg is False and mkpng is False and plot:
            mpl.show()

    def delay_spread_export(self, txgrp: list = [-1], rxgrp: list = [-1], mkpng: bool = False,
               cmap: str = 'Blues', floor: float = -200.0, ceil: float = np.nan,
               nff: bool = True, matsav: bool = False, plot: bool = True, show: bool =True, fidbase: int = 0,
               title: str = ''):


        for i in self.source.txs.items(grprange=txgrp):
            ds_avg = list()
            ds_avg_rms = list()
            for j in self.source.rxs.items(grprange=rxgrp):
                delay = list()
                power = list()

                if i.chan_to(j[1]):
                    for k in i.chan_to(j[1]).paths.items():
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
    enable_latex(18)
    DS = DataStorage(conf='dbconf.txt', dbname='BUSMOD')
    cir = CIR(DS)
    cir.export_single(txgrp=[-1], rxgrp=[2], rxrange=[0,4,8,12,16,20,24,28,32,36], nff=False, matsav=True, plot=True, mkpng=False, floor=-140.0, mksvg=False)
    cir.export_single(txgrp=[-1], rxgrp=[4], rxrange=[40,44,48,52,56,60,64,68,72,76], nff=False, matsav=True, plot=True, mkpng=False, floor=-140.0, mksvg=False)
    #cir.export(txgrp=-1, rxgrp=4, nff=True, matsav=False, plot=True, mkpng=False, zmin=-190.0, zmax=-40.0)
    #cir.export(txgrp=-1, rxgrp=2, nff=True, matsav=False, plot=True, mkpng=False, zmin=-190.0, zmax=-40.0)
    exit()
