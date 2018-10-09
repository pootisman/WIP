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
import numpy as np
from pairdata import DataStorage
from auxfun import l2db
import scipy.io as sio

from itertools import product

COLORS = ['r', 'g', 'b', 'k']
DASHES = ['-', '--', ':', '-.']

LINESTYLES = list(product(DASHES, COLORS))


class DelaySpreadPlot:
    '''
    Delay spread plotter, allow you to plot, compare and build histograms of the delay spread for channels between
    different receiver/transmitter pairs
    '''
    def __init__(self, source: DataStorage):
        self.source = source
        self.xdata = list()
        self.ydata = list()
        self.data = dict()

        self.title = str()

        self.xmax = np.NINF
        self.xmin = np.PINF

        self.ymax = np.NINF
        self.ymin = np.PINF

        self.tx = None
        self.rxgrp = int()

    def plot_group(self, txgrp: list = [-1], rxgrp: list = [-1], txrange: list = [-1], rxrange: list = [-1],
                   title: str = '', plot: bool = True, mkpng: bool = True, fidbase: int = 0, matsav: bool = True,
                   csvsav: bool = True, show: bool = True, rx_name_map: dict = None, tx_name_map: dict = None,
                   mkpdf: bool = False):
        '''
        Plots delay spread value vs receiver number between specific group of RX&TX, only one group is plotted
        :param txgrp: TX group to use. Ignored if set to -1
        :param rxgrp: RX group to use. Ignored if set to -1
        :param txrange: TX ids that should be used. Set txgrp to -1
        :param rxrange: RX ids that should be used. Set rxgrp to -1
        :param title: Prefix for title of the plot
        :param plot: Should we produce the plot?
        :param mkpng: Should we make PNG of the plot and save it?
        :param fidbase: Figure ID starting number. Used in order to avoid conflicts between figures
        :param matsav: Save data to .mat file for use MATLAB later on?
        :param csvsav: Save data to .csv file?
        :param show: Show figure window?
        :param rx_name_map: Maps certain RX groups to string names [cosmetic]
        :param tx_name_map: Maps certain TX groups to string names [cosmetic]
        :return:
        '''

        assert self.source is not None, 'Source needed, refusing to plot!'

        if txrange[0] == -1:
            txrange = self.source.txs.keys()

        if rxrange[0] == -1:
            rxrange = self.source.rxs.keys()

        self.title = title

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        for i in txrange:
            if self.source.txs[i].setid in txgrp or txgrp[0] == -1:
                self.xdata = list()
                self.ydata = list()

                self.data = dict()

                for j in rxrange:
                    if self.source.rxs[j].setid in rxgrp or rxgrp[0] == -1:
                        try:
                            self.xdata.append(j)
                            self.ydata.append(self.source.txs[i].chan_to(self.source.rxs[j]).delay_spread * 1e9)
                        except AttributeError:
                            print('No path between {} and {}'.format(i, j))

                self.xmax = np.nanmax(self.xdata)
                self.xmin = np.nanmin(self.xdata)

                self.ymax = np.nanmax(self.ydata)
                self.ymin = np.nanmin(self.ydata)

                self.tx = i
                self.rxgrp = rxgrp[0]

                if plot or mkpng:
                    f = mpl.figure(fidbase + i)

                    self.data['{}_X_{}'.format(title.strip(' '), j if rx_name_map is None else rx_name_map[j])] \
                        = self.xdata
                    self.data['{}_Y_{}'.format(title.strip(' '), j if rx_name_map is None else rx_name_map[j])] \
                        = self.ydata

                    mpl.stem(self.xdata, self.ydata, label='RX grp {}'.
                             format(j if rx_name_map is None else rx_name_map[j]))
                    mpl.xlabel('RX Position')
                    mpl.ylabel('Delay Spread, [ns]')
                    mpl.title('{}Delay spread [TX {} $\\rightarrow$ RXg {}]'.format(title, i, rxgrp))
                    mpl.grid(linestyle='--')
                    mpl.legend()
                    mpl.tight_layout()

                if mkpng:
                    mpl.savefig('{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.png'.format(i, rxgrp[0], title))
                    mpl.close(f)

                if mkpdf:
                    mpl.savefig('{2}CIR3D_tx{0:03d}_rxgrp{1:03d}.pdf'.format(i, rxgrp[0], title))
                    mpl.close(f)

                if matsav:
                    sio.savemat('{2}DS_tx{0:03d}_rxgrp{1:03d}.mat'.format(i, rxgrp[0], title),
                                {'rx': self.xdata, 'delay_spread': self.ydata})

                if csvsav:
                    file = open('delay_spread.csv', mode='w')
                    file.write('Rec. num. ,  Delay spread [ns]\n')
                    for k in range(self.ydata.__len__()):
                        file.write('{},{}\n'.format(self.xdata[k], self.ydata[k]))
                    file.close()

        if mkpng is False and plot and show:
            mpl.show()

    def plot_groups(self, txgrp: list = [-1], rxgrp: list = [-1], txrange: list = [-1], rxrange: list = [-1],
                    title: str = '', mkpng: bool = False, fidbase: int = 0, matsav: bool = False, csvsav: bool = False,
                    show: bool = True, overlay: bool = True, ymin: float = np.NINF, ymax: float = np.NINF,
                    rx_name_map: dict = None, tx_name_map: dict = None, dispylabel: bool = True,
                    dispxlabel: bool = True, disptitle: bool = True, mkpdf: bool = False):
        '''
        Plot multiple delay spread curves
        '''
        assert self.source is not None, 'Source needed, refusing to plot!'

        if txrange[0] == -1:
            txrange = self.source.txs.keys()

        if rxrange[0] == -1:
            rxrange = self.source.rxs.keys()

        self.title = title

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        f = mpl.figure(fidbase)

        plot_index = 0

        for i in txrange:
            if self.source.txs[i].setid in txgrp or txgrp[0] == -1:
                for j in rxgrp:
                    self.xdata = list()
                    self.ydata = list()
                    for k in rxrange:
                        if self.source.rxs[k].setid == j or rxgrp[0] == -1:
                            try:
                                self.xdata.append(k)
                                self.ydata.append(self.source.txs[i].chan_to(self.source.rxs[k]).delay_spread * 1e9)
                            except AttributeError:
                                print('No path between {} and {}'.format(i, k))

                    self.xmax = np.max([np.nanmax(self.xdata), self.xmax])
                    self.xmin = np.min([np.nanmin(self.xdata), self.xmin])

                    self.ymax = np.max([np.nanmax(self.ydata), self.ymax])
                    self.ymin = np.min([np.nanmin(self.ydata), self.ymin])

                    self.tx = i
                    self.rxgrp = rxgrp[0]

                    if overlay:
                        mpl.plot(list(np.linspace(start=0, stop=self.ydata.__len__(), num=self.ydata.__len__(),
                                                  endpoint=True)), self.ydata, ''.join(LINESTYLES[plot_index]),
                                 label='RX grp {}'.format(j if rx_name_map is None else rx_name_map[j]))
                        self.data['{}_X_{}'.format(title.strip(' '), j if rx_name_map is None else rx_name_map[j])]\
                            = list(np.linspace(start=0, stop=self.ydata.__len__(),
                                               num=self.ydata.__len__(), endpoint=True))
                        self.data['{}_Y_{}'.format(title.strip(' '), j if rx_name_map is None else rx_name_map[j])]\
                            = self.ydata
                    else:
                        self.data['{}_X_{}'.format(title.strip(' '), j if rx_name_map is None else rx_name_map[j])]\
                            = self.xdata
                        self.data['{}_Y_{}'.format(title.strip(' '), j if rx_name_map is None else rx_name_map[j])]\
                            = self.ydata
                        mpl.plot(self.xdata, self.ydata, ''.join(LINESTYLES[plot_index]),
                                 label='RX grp {}'.format(j if rx_name_map is None else rx_name_map[j]))

                    plot_index += 1

        mpl.xlim([self.xmin, self.xmax if not overlay else self.ydata.__len__()])
        mpl.ylim([0.0 if ymin == np.NINF else ymin, (self.ymax + 0.1 * self.ymax) if ymax == np.NINF else ymax])
        if dispxlabel:
            mpl.xlabel('RX Position')
        if dispylabel:
            mpl.ylabel('Delay Spread, [ns]')
        if disptitle:
            mpl.title('{}delay spread [TXg{} $\\rightarrow$ RXg{}]'.format(title, txgrp if txgrp[0] != -1 else 'All',
                                                                           rxgrp))
        mpl.grid(linestyle='--')
        mpl.legend()
        mpl.tight_layout()

        if mkpng:
            mpl.savefig('DelSprd_TX_{}.png'.format(title))
            mpl.close(f)

        if mkpdf:
            mpl.savefig('DelSprd_TX_{}.pdf'.format(title))
            mpl.close(f)


        if matsav:
            sio.savemat('{2}DS_tx{0:03d}_rxgrp{1:03d}.mat'.format(i, rxgrp[0], title), self.data)

        if csvsav:
            file = open('delay_spread.csv', mode='w')
            file.write('Rec. num. ,  Delay spread [ns]\n')
            for k in range(self.ydata.__len__()):
                file.write('{},{}\n'.format(self.xdata[k], self.ydata[k]))
            file.close()

        if mkpng is False and show:
            mpl.show()

    def collate(self, other):
        assert isinstance(other, DelaySpreadPlot), 'Can\'t collate DelaySpreadPlot with something other!'
