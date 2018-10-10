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
import matplotlib.pyplot as mpl
from matplotlib.colors import ListedColormap
from cir import CIR
from auxfun import basint3


class CIRCollator:
    def __init__(self):
        self.xmin = np.PINF
        self.xmax = np.NINF
        self.ymin = np.PINF
        self.ymax = np.NINF
        self.zmin = np.PINF
        self.zmax = np.NINF

        self.xdim = 0
        self.ydim = 0

        self.datas = list()
        self.cmap_ids = ['Blues', 'Reds', 'Purples']
        self.titles = list()

    def __add__(self, other):
        if isinstance(other, CIR):
            self.datas.append(other)
        elif isinstance(other, list):
            for i in other:
                if isinstance(i, CIR):
                    self.datas.append(i)
        else:
            raise TypeError('List or CIRs needed, got {}'.format(other))

    def export_collated(self, mkpng: bool = False, plot: bool = True, show: bool =True, fidbase: int = 0,
                        title: str = '', cmaplist: list = None, idxs: list = None, csq: bool = False,
                        csqloc: list = [0.3, 0.65, 0.3, 0.3], xlabel: bool = True, ylabel: bool = True,
                        title_draw: bool = True, figsize: tuple = (8, 6), yrange: tuple = None, xtics: int = 5,
                        ytics: int = 5, mkpdf: bool = False, fnappend: str = '', xdraw: bool = True,
                        ydraw: bool = True, mainlabels_fs: int = 12):
        self.titles = list()

        if idxs is None:
            idxs = [j for j in range(self.datas.__len__())]

        if not isinstance(idxs, list):
            idxs = [idxs]

        # Determine min and max of all datasets
        # Replace with lambdas?
        for i in self.datas:
            if self.xmin > i.xmin:
                self.xmin = i.xmin
            if self.xmax < i.xmax:
                self.xmax = i.xmax
            if self.ymin > i.ymin:
                self.ymin = i.ymin
            if self.ymax < i.ymax:
                self.ymax = i.ymax
            if self.zmin > i.zmin:
                self.zmin = i.zmin
            if self.zmax < i.zmax:
                self.zmax = i.zmax
            if self.xdim < i.xdim:
                self.xdim = i.xdim
            if self.ydim < i.ydim:
                self.ydim = i.ydim

        f = mpl.figure(fidbase, figsize=figsize)

        if cmaplist is not None:
            self.cmap_ids = cmaplist

        iind = 0

        for j in idxs:
            i = self.datas[j]
            self.titles.append(i.title)

            (x, y, z) = basint3(x=i.xdata, y=i.ydata, xc=self.xdim, yc=self.ydim, zmin=self.zmin, xmin=self.xmin,
                                xmax=self.xmax, ymin=self.ymin, ymax=self.ymax,
                                z=[j * (j > i.zmin) + 2.0 * self.zmin * (j <= i.zmin) for j in i.zdata])
            [x, y] = np.meshgrid(x, y)

            if plot or mkpng:
                cmap_def = mpl.get_cmap(self.cmap_ids[iind])
                cmap_cus = cmap_def(np.arange(cmap_def.N))
                cmap_cus[:, -1] = np.linspace(0, 1 - 0.3 * (iind > 0), cmap_def.N)
                cmap_cus = ListedColormap(cmap_cus)

                mpl.pcolor(np.transpose(x), np.transpose(y), z, clip_on=True,
                           cmap=cmap_cus, vmin=self.zmin, vmax=self.zmax)

            iind += 1

        mpl.clim(vmin=self.zmin, vmax=self.zmax)

        if yrange:
            mpl.yticks(np.linspace(start=yrange[0], num=xtics, stop=yrange[1], endpoint=True).tolist())
        else:
            mpl.yticks(np.linspace(start=self.ymin, num=ytics, stop=self.ymax, endpoint=True).tolist())

        mpl.xticks(np.linspace(start=self.xmin, num=xtics, stop=self.xmax, endpoint=True).tolist())

        if xlabel:
            mpl.xlabel('RX Position', fontsize=mainlabels_fs)

        if ylabel:
            mpl.ylabel('Delay, [ns]', fontsize=mainlabels_fs)

        if title_draw:
            mpl.title('{}CIR\\@TX \\#{}'.format(title, 'vs '.join(self.titles)))

        if not xdraw:
            mpl.gca().set_xticklabels([])

        if not ydraw:
            mpl.gca().set_yticklabels([])

        mpl.tight_layout()

        if csq:
            axins = f.add_axes(csqloc)

            cmap = mpl.get_cmap(self.cmap_ids[0])
            cmap1 = cmap(np.arange(cmap.N))
            cmap1[:, -1] = np.linspace(start=0, stop=1, num=cmap.N, endpoint=True)

            cmap = mpl.get_cmap(self.cmap_ids[1])
            cmap2 = cmap(np.arange(cmap.N))
            cmap2[:, -1] = np.linspace(start=0, stop=0.7, num=cmap.N, endpoint=True)

            csq_valsx = np.tile(cmap1, [cmap.N, 1, 1])
            csq_valsy = np.tile(cmap2, [cmap.N, 1, 1])

            axins.imshow(csq_valsx)
            axins.imshow(np.rot90(csq_valsy, 1))

            axins.tick_params(direction='in')

            axins.set_xlabel('{} pow. [dBm]'.format(self.titles[0]))
            axins.set_ylabel('{} pow. [dBm]'.format(self.titles[1]))

            axins.set_xticks(np.linspace(start=0.0, num=5, stop=cmap.N, endpoint=True))
            axins.set_yticks(np.linspace(start=0.0, num=5, stop=cmap.N, endpoint=True))

            powerticks = np.linspace(start=self.zmin, num=5, stop=self.zmax, endpoint=True).tolist()
            powerticks = [str(i) for i in powerticks]
            axins.set_xticklabels(powerticks, rotation=-45, ha='left')
            axins.set_yticklabels(reversed(powerticks))

        if mkpng:
            mpl.savefig('Collated[{}] CIR3D {}.png'.format('|'.join(self.titles).strip(' '), fnappend))

        if mkpdf:
            mpl.savefig('Collated[{}] CIR3D {}.pdf'.format('|'.join(self.titles).strip(' '), fnappend))

        if show:
            mpl.show()
        else:
            mpl.close(f)

    def export_adjusted(self, mkpng: bool = False, plot: bool = True, show: bool =True, fidbase: int = 0,
                        title: str = '', cmaplist: list = None, idxs: list = [0]):

        if not isinstance(idxs, list):
            idxs = [idxs]


        # Determine min and max of all datasets
        for i in self.datas:
            if self.xmin > i.xmin:
                self.xmin = i.xmin
            if self.xmax < i.xmax:
                self.xmax = i.xmax
            if self.ymin > i.ymin:
                self.ymin = i.ymin
            if self.ymax < i.ymax:
                self.ymax = i.ymax
            if self.zmin > i.zmin:
                self.zmin = i.zmin
            if self.zmax < i.zmax:
                self.zmax = i.zmax
            if self.xdim < i.xdim:
                self.xdim = i.xdim
            if self.ydim < i.ydim:
                self.ydim = i.ydim

        if cmaplist is not None:
            self.cmap_ids = cmaplist

        iind = 0

        for k in idxs:
            i = self.datas[k]
            f = mpl.figure(fidbase + i)

            (x, y, z) = basint3(x=i.xdata, y=i.ydata, xc=self.xdim, yc=self.ydim, zmin=self.zmin, xmin=self.xmin,
                                xmax=self.xmax, ymin=self.ymin, ymax=self.ymax,
                                z=[j * (j > i.zmin) + 2.0 * self.zmin * (j <= i.zmin) for j in i.zdata])
            [x, y] = np.meshgrid(x, y)

            if plot or mkpng:
                cmap_def = mpl.get_cmap(self.cmap_ids[iind])
                cmap_cus = cmap_def(np.arange(cmap_def.N))
                cmap_cus[:, -1] = np.linspace(1, 0, cmap_def.N)
                cmap_cus = ListedColormap(cmap_cus)

                mpl.pcolor(np.transpose(x), np.transpose(y), z, clip_on=True, alpha=(1.0 if iind == 0 else 0.5),
                           cmap=cmap_cus, vmin=self.zmin, vmax=self.zmax)

                mpl.colorbar()
                mpl.clim(vmin=self.zmin, vmax=self.zmax)
                mpl.xlabel('RX Position')
                mpl.ylabel('Delay, [ns]')
                mpl.title('{}CIR\@[TX {} $\\rightarrow$ RXg {}]'.format(title, i.title, i.rxgrp))
                mpl.tight_layout()

            if mkpng:
                mpl.savefig('{2}CIR3D\_tx{0:03d}\_rxgrp{1:03d}.png'.format(i.tx, i.rxgrp, i.title))
                mpl.close(f)

            iind += 1

        if show:
            mpl.show()
