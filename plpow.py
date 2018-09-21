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
from phys_path_procs import *

class plPlot():
    def __init__(self, source: pairdata.data_stor = None):
        assert source is not None
        self.xdata = None
        self.ydata = None
        self.nsamps = 0
        self.source = source
        self.a = 0.0
        self.b = 0.0
        self.typ = None
        self.thrshld = 0.0

    def regr_comp(self, rxgrp: list = [-1], txgrp: list = [-1], typ: str ='LOS', threshold: float = -130, nff: bool = True):
        self.xdata = []
        self.ydata = []
        self.typ = typ
        self.thrshld = threshold

        txgrp = [txgrp] if not isinstance(txgrp, list) else txgrp
        rxgrp = [rxgrp] if not isinstance(rxgrp, list) else rxgrp

        # Prepare data for making the regression
        for i in self.source.txs.items():
            # TX is in valid group or group is ignored?
            if i[1].setid in txgrp or txgrp[0] == -1:
                # Which RXes are reached by TX
                for j in i[1].chans_to_pairs.items():
                    # Destination in a valid group or group is ignored?
                    if j[0].setid in rxgrp or rxgrp[0] == -1:
                        # Check paths for the RX-TX, only pick valid ones
                        for k in j[1].paths.items():
                            if nff and not k[1].near_field_failed:
                                if k[1].interactions.__len__() == 0 and typ == 'LOS' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                                elif k[1].interactions.__len__() == 1 and typ == 'NLOS-1' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                                elif k[1].interactions.__len__() == 2 and typ == 'NLOS-2' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                                elif k[1].interactions.__len__() >= 1 and typ == 'NLOS' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                            elif not nff:
                                if k[1].interactions.__len__() == 0 and typ == 'LOS' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                                elif k[1].interactions.__len__() == 1 and typ == 'NLOS-1' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                                elif k[1].interactions.__len__() == 2 and typ == 'NLOS-2' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1
                                elif k[1].interactions.__len__() >= 1 and typ == 'NLOS' and l2db(k[1].pow) >= self.thrshld:
                                    self.xdata.append(np.log10(j[1].dist))
                                    self.ydata.append(k[1].FSPL)
                                    self.nsamps += 1

        self.xdata = np.asarray(self.xdata)
        self.ydata = np.asarray(self.ydata)
        # Data ready, calculate ax + b regression
        self.a, self.b = np.linalg.lstsq(np.vstack([self.xdata, np.ones(self.xdata.__len__())]).T, self.ydata, rcond=None)[0]
        return self.a, self.b

    def export(self, plot: bool = True, csvsav: bool = False, matsav: bool = False):
        self.mean_res = np.mean(self.ydata - self.b - self.a * self.xdata)
        self.var_res = np.sqrt(np.var(self.ydata - self.b - self.a * self.xdata))

        if plot:
            mpl.figure()
            mpl.plot(np.power(10, self.xdata), self.ydata, 'r.', label='Experimental data')
            mpl.plot(np.power(10, self.xdata), self.b + self.a * self.xdata, 'b.', label='Fitted model')
            mpl.title('{} Regression b={:3.2}dBm, a={:3.2}dBm, m={:3.2}dBm , s={:3.2}dBm'.format(self.typ, self.a, self.b,
                                                                                        self.mean_res, self.var_res))
            mpl.xlabel('Distance, [meters]')
            mpl.ylabel('Pathloss, [dB]')
            mpl.grid(linestyle='--')
            mpl.xlim([0, np.nanmax(np.power(10.0, self.xdata))])
            mpl.show()

        if csvsav:
            file = open('FSPL_regr.csv', mode='w')
            file.write('Dist. [m],FSPL [dB]\n')
            for k in range(self.ydata.__len__()):
                file.write('{},{}\n'.format(self.xdata[k], self.ydata[k]))
            file.write('A,B\n')
            file.write('{},{}\n'.format(self.a, self.b))
            file.close()

        if matsav:
            sio.savemat('FSPL_regr.mat', {'A': self.a, 'B': self.b, 'xdata': self.xdata, 'ydata': self.ydata})


if __name__ == '__main__':
    DS = pairdata.data_stor(conf='dbconf.txt')
    DS.load_rxtx('class_sqlite')
    DS.load_paths()
    DS.load_interactions()
    check_data_NF(DS)
    plp = plPlot(source=DS)
    print(plp.regr_comp(typ='NLOS', threshold=-115, nff=True))
    plp.export(csvsav=True, matsav=True)
    exit()
