import numpy as np
import matplotlib.pyplot as mpl
import pairdata
from auxfun import *
import numba

class plPlot():
    def __init__(self, source: pairdata.data_stor = None):
        assert source is not None
        self.xdata = None
        self.ydata = None
        self.nsamps = 0
        self.source = source
        self.a = 0
        self.b = 0
        self.typ = None

    @numba.jit
    def regr_comp(self, rxgrp: int = -1, txgrp: int = -1, typ: str ='LOS', threshold: int = -130):
        self.xdata = []
        self.ydata = []
        self.typ = typ
        self.thrshld = threshold

        # Prepare data for making the regression
        for i in self.source.txs.items():
            # TX is in valid group or group is ignored?
            if i[1].setid == txgrp or txgrp == -1:
                # Which RXes are reached by TX
                for j in i[1].chan_to_pairs.items():
                    # Destination in a valid group or group is ignored?
                    if j[0].setid == rxgrp or rxgrp == -1:
                        # Check paths for the RX-TX, only pick valid ones
                        for k in j[1].paths.items():
                            if k[1].interactions.__len__() == 0 and typ == 'LOS' and k[1].pow > self.thrshld:
                                self.xdata.append(np.log10(j[1].dist))
                                self.ydata.append(k[1].FSPL)
                                self.nsamps += 1
                            elif k[1].interactions.__len__() == 1 and typ == 'nLOS-1' and k[1].pow > self.thrshld:
                                self.xdata.append(np.log10(j[1].dist))
                                self.ydata.append(k[1].FSPL)
                                self.nsamps += 1
                            elif k[1].interactions.__len__() == 2 and typ == 'nLOS-2' and k[1].pow > self.thrshld:
                                self.xdata.append(np.log10(j[1].dist))
                                self.ydata.append(k[1].FSPL)
                                self.nsamps += 1
                            elif k[1].interactions.__len__() >= 1 and typ == 'nLOS-*' and k[1].pow > self.thrshld:
                                self.xdata.append(np.log10(j[1].dist))
                                self.ydata.append(k[1].FSPL)
                                self.nsamps += 1

        self.xdata = np.asarray(self.xdata)
        self.ydata = np.asarray(self.ydata)
        # Data ready, calculate ax + b regression
        self.a, self.b = np.linalg.lstsq(np.vstack([self.xdata, np.ones(self.xdata.__len__())]).T, self.ydata, rcond=None)[0]
        return self.a, self.b

    def plot_reg(self):
        mpl.figure()
        mpl.plot(np.power(10, self.xdata), self.ydata, 'r.', label='Experimental data')
        mpl.plot(np.power(10, self.xdata), self.b + self.a * self.xdata, 'b.', label='Fitted model')
        self.mean_res = np.mean(self.ydata - self.b - self.a * self.xdata)
        self.var_res = np.sqrt(np.var(self.ydata - self.b - self.a * self.xdata))
        mpl.title('{} Regression b={}dBm, a={}dBm, m={}dBm , s={}dBm'.format(self.typ, self.a, self.b, self.mean_res, self.var_res, width=6))
        mpl.xlabel('Distance, [meters]')
        mpl.ylabel('Pathloss, [dB]')
        mpl.grid()
        mpl.xlim([0, np.max(np.power(10.0, self.xdata))])
        mpl.show()


if __name__ == '__main__':
    DS = pairdata.data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TESTe/Class@60GHz.TESTe.sqlite')
    DS.load_path()
    plp = plPlot(source=DS)
    print(plp.regr_comp(typ='LOS'))
    plp.plot_reg()
    exit()
