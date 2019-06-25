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
import matplotlib.ticker as tck
import scipy.io as sio
from pairdata import DataStorage
from phys_path_procs import check_data_nf
from auxfun import l2db, markers, excl_interactions, styles, colors


class PLPlot:
    def __init__(self, source: DataStorage = None):
        assert source is not None
        self.xdata = list()
        self.ydata = list()
        self.raw_xdata = list()
        self.raw_ydata = list()
        self.pair_info = list()
        self.raw_pair_info = list()
        self.raw_dist_info = list()
        self.dist_info = list()
        self.nsamps = 0
        self.source = source
        self.a = 0.0
        self.b = 0.0
        self.typ = None
        self.thrshld = 0.0
        self.mean_res = 0.0
        self.var_res = np.PINF

    def __sigma_rm(self, sigmamul: float = 1):
        mask = np.abs(self.ydata - (self.a + self.b * self.xdata)) < np.sqrt(np.var(self.ydata)) * sigmamul
        self.ydata = self.ydata[mask]
        self.xdata = self.xdata[mask]
        self.pair_info = self.pair_info[mask]
        self.dist_info = np.extract(mask, self.dist_info)

    def regr_comp(self, rxgrp: list = [-1], txgrp: list = [-1], typ: str ='LOS', threshold: float = -130,
                  nff: bool = True, xrange: tuple = (0.0, 1.0), sigma_f: float = 1, sigma_rounds: int = 1,
                  additive: bool = False):
        self.typ = typ
        self.thrshld = threshold

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        # TODO: Figure out what invalid Electric fields, delays and phases mean...
        # Prepare data for making the regression
        for i in self.source.txs.items():
            # TX is in valid group or group is ignored?
            if i[1].setid in txgrp or txgrp[0] == -1:
                # Which RXes are reached by TX
                for j in i[1].chans_to_pairs.items():
                    # Destination in a valid group or group is ignored?
                    if j[0].setid in rxgrp or rxgrp[0] == -1:
                        # Check paths for the RX-TX, only pick valid ones
                        best_pow = 0.0
                        best_dist = 0.0
                        best_fspl = 0.0

                        if typ != 'FULL':
                            for k in j[1].paths.items():
                                # Near field condition check is now integrated into the if-blocks, it is ignored if nff
                                # is False
                                if k[1].interactions.__len__() == 0 and typ == 'LOS' and \
                                        l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) and \
                                        (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif k[1].interactions.__len__() == 1 and typ == 'NLOS-1' and \
                                        l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                        and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif k[1].interactions.__len__() == 2 and typ == 'NLOS-2' and \
                                        l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                        and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif k[1].interactions.__len__() >= 1 and typ == 'NLOS' and \
                                        l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                        and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif excl_interactions(k[1].interactions) == 0 and typ == 'LOS-pen' and \
                                        l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                        and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif excl_interactions(k[1].interactions) == 1 and typ == 'NLOS-1-pen' and \
                                         l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                         and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif excl_interactions(k[1].interactions) == 2 and typ == 'NLOS-2-pen' and \
                                         l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                         and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                                elif excl_interactions(k[1].interactions) >= 1 and typ == 'NLOS-pen' and \
                                         l2db(k[1].pow) >= self.thrshld and (best_pow < k[1].pow or additive) \
                                         and (not nff or not k[1].near_field_failed):
                                    best_dist = j[1].dist
                                    best_fspl = k[1].fspl
                                    if additive:
                                        best_pow += k[1].pow
                                    else:
                                        best_pow = k[1].pow
                        else:
                            if l2db(j[1].pow) >= self.thrshld:
                                best_pow = j[1].pow
                                if additive:
                                    best_pow += k[1].pow
                                else:
                                    best_pow = k[1].pow

                        if best_pow > 0.0:
                            self.raw_xdata.append(np.log10(best_dist))
                            self.raw_ydata.append(-l2db(best_pow))
                            self.nsamps += 1
                            self.raw_pair_info.append(j[0].setid)
                            self.raw_dist_info.append(xrange[0] < best_dist < xrange[1])



        # Convert NumPy arrays for further processing
        self.xdata = np.asarray(self.raw_xdata)
        self.ydata = np.asarray(self.raw_ydata)

        self.raw_xdata = np.asarray(self.raw_xdata)
        self.raw_ydata = np.asarray(self.raw_ydata)

        self.pair_info = np.asarray(self.raw_pair_info)
        self.dist_info = np.asarray(self.raw_dist_info)

        self.b, self.a = np.linalg.lstsq(np.vstack([self.xdata[self.dist_info],
                                                    np.ones(self.xdata[self.dist_info].__len__())]).T,
                                         self.ydata[self.dist_info], rcond=None)[0]
        if sigma_f < np.PINF:
            for i in range(sigma_rounds):
                self.__sigma_rm(sigma_f)
                self.b, self.a = np.linalg.lstsq(np.vstack([self.xdata[self.dist_info],
                                                        np.ones(self.xdata[self.dist_info].__len__())]).T,
                                                 self.ydata[self.dist_info], rcond=None)[0]

        self.residues = self.raw_ydata - self.a - self.b * self.raw_xdata
        self.mean_res = np.mean(self.residues)
        self.var_res = np.sqrt(np.var(self.residues))

        return (self.a, self.b, self.mean_res, self.var_res, self.nsamps, self.typ, self.xdata, self.ydata,
                np.asarray(self.raw_xdata), np.asarray(self.raw_ydata), self.residues)

    def export(self, plot: bool = True, csvsav: bool = False, matsav: bool = False, split_points: bool = True,
               trajmap: dict = dict(), title_draw: bool = False, mkpdf: bool = False, mkpng: bool = False,
               show: bool = True, fname_suffix: str = "", draw_resudues: bool = True):

        if plot:
            loc = tck.AutoLocator()
            loc.create_dummy_axis()
            ytics = loc.tick_values(np.min([self.ydata.min(), self.a + self.b * self.xdata.min()]),
                                    np.max([self.ydata.max(), self.a + self.b * self.xdata.max()]))
            xtics = loc.tick_values(np.power(10.0, self.xdata).min(), np.power(10.0, self.xdata.max()))
            x_mod_data = np.linspace(start=np.min(xtics), stop=np.max(xtics), endpoint=True, num=20)
            f = mpl.figure()
            if split_points is False:
                mpl.plot(np.power(10.0, self.xdata), self.ydata, 'r.', label='SBR data')
            else:
                grps = np.unique(self.pair_info)
                for v in grps:
                    pi = np.where(np.asarray(self.pair_info) == v)
                    mpl.plot(np.power(10.0, self.xdata[pi]), self.ydata[pi], markers[v], label='SBR data, group {}'.
                             format(trajmap[v]))

            mpl.plot(x_mod_data, self.a + self.b * np.log10(x_mod_data), 'b-', label='Fitted model')

            if title_draw:
                mpl.title('{} Regression $\\alpha$={:3.2f}dBm, $\\beta$={:3.2f}dBm,\\'
                          '\\ m={:3.2f}dBm, $\\sigma$={:3.2f}dBm, N={} samples'.
                          format(self.typ, self.a, self.b/10.0, self.mean_res, self.var_res, self.nsamps))

            mpl.xlabel('Distance, [meters]')
            mpl.ylabel('Path loss, [dB]')
            mpl.xticks(xtics)
            mpl.yticks(ytics)
            mpl.grid(linestyle='--')
            mpl.legend()
            mpl.xlim([np.min(xtics), np.max(xtics)])
            mpl.ylim([np.min(ytics), np.max(ytics)])
            mpl.tight_layout()
            if show:
                mpl.show()

            if mkpng:
                mpl.savefig('Regr_{}{}.png'.format(self.typ, fname_suffix))

            if mkpdf:
                mpl.savefig('Regr_{}{}.pdf'.format(self.typ, fname_suffix))

        if csvsav:
            file = open('fspl_regr.csv', mode='w')
            file.write('Dist. [m],fspl [dB]\n')
            for k in range(self.ydata.__len__()):
                file.write('{},{}\n'.format(self.xdata[k], self.ydata[k]))
            file.write('A,B\n')
            file.write('{},{}\n'.format(self.a, self.b))
            file.close()

        if matsav:
            sio.savemat('{}_{}_regr.mat'.format(self.typ, fname_suffix), {'A': self.a, 'B': self.b, 'xdata': self.xdata,
                                                                          'ydata': self.ydata, 'grps': self.pair_info,
                                                                          'disti': self.dist_info,
                                                                          'ydata_raw': self.raw_ydata,
                                                                          'xdata_raw': self.raw_xdata,
                                                                          'residues': self.residues})

        if not show:
            mpl.close(f)

    def plot_regr_only(self, params: dict = {}, title_draw: bool = False, dists: tuple = (0.1, 1.0),
                       show_data: bool = True, fspl_draw: bool = False, plot_raw: bool = False,
                       plot_resid: bool = True):
        resid = []

        mpl.figure()

        loc = tck.AutoLocator()
        loc.create_dummy_axis()
        xtics = loc.tick_values(dists[0], dists[1])
        x_mod_data = np.linspace(start=np.min(xtics), stop=np.max(xtics), endpoint=True, num=20)

        ylims = [np.PINF, np.NINF]

        if fspl_draw:
            mpl.plot(x_mod_data, 20 * np.log10(4.0 * np.pi * x_mod_data * 60e9 / 3e8), ':', label='FSPL equation', linewidth=4.0)

        n = 0
        for i in params.keys():
            data = params[i][0] + params[i][1] * np.log10(x_mod_data)

            ylims[0] = np.min(data) if np.min(data) < ylims[0] else ylims[0]
            ylims[1] = np.max(data) if np.max(data) > ylims[1] else ylims[1]

            for j in params[i][10]:
                resid.append(j)

            mpl.plot(x_mod_data, data, ''.join([styles[n], colors[n]]), linewidth=2, label='Fitted model [{}]'.format(i))
            if show_data:
                if plot_raw:
                    mpl.plot(10.0**(params[i][8]), params[i][9], ''.join([markers[n], colors[n]]), linewidth=2, label='Simulation [{}]'.format(i))
                else:
                    mpl.plot(10.0**(params[i][6]), params[i][7], ''.join([markers[n], colors[n]]), linewidth=2, label='Simulation [{}]'.format(i))


            if title_draw:
                mpl.title('{} Regression $\\alpha$={:3.2f}dBm, $\\beta$={:3.2f}dBm,\\'
                          '\\ m={:3.2f}dBm, $\\sigma$={:3.2f}dBm, N={} samples'.
                          format(self.typ, self.a, self.b / 10.0, self.mean_res, self.var_res, self.nsamps))

            mpl.xlabel('Distance, [meters]')
            mpl.ylabel('Path loss, [dB]')
            mpl.xticks(xtics)
            mpl.grid(linestyle='--')
            mpl.legend()
            mpl.xlim([np.min(xtics), np.max(xtics)])
            n+=1



        ytics = loc.tick_values(ylims[0], ylims[1])
        mpl.yticks(ytics)
        mpl.ylim([np.min(ytics), np.max(ytics)])
        mpl.tight_layout()

        if plot_resid:
            mpl.figure()
            mpl.hist(resid, bins=30)
            mpl.xlabel('Residue, [dB]')
            mpl.ylabel('Hits')
            mpl.grid(linestyle='--')

        mpl.show()


if __name__ == '__main__':
    from auxfun import enable_latex
    enable_latex(pt=14)

    def loadall(conf='dbconf.txt', name: str = 'BUSTRX'):
        DS = DataStorage(conf=conf, dbname=name)
        DS.load_rxtx()
        DS.load_paths()
        DS.load_interactions()
        return DS

    typd = 'NLOS-1'

    dd = {}
    dd2 = {}

    sr = 1
    add = False
    fspl_draw = False
    td = False
    drraw = True

    DS = loadall()
    plp = PLPlot(source=DS)
    data = plp.regr_comp(rxgrp=[4], typ=typd, threshold=-100, nff=False, xrange=(0.5, 12), sigma_rounds=sr,
                         additive=add)
    plp.export(csvsav=False, matsav=True, plot=True, trajmap={4: 'Low'}, mkpdf=True,
              show=False, fname_suffix='low')

    print(data[0], data[1], data[3])

    dd['l'] = data

    plp = PLPlot(source=DS)
    data = plp.regr_comp(rxgrp=[2], typ=typd, threshold=-100, nff=False, xrange=(0.5, 12), sigma_rounds=sr,
                         additive=add)
    plp.export(csvsav=False, matsav=True, plot=True, trajmap={2: 'High'}, mkpdf=True,
              show=False, fname_suffix='hi')

    print(data[0], data[1], data[3])

    dd['h'] = data

    exit()
