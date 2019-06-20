import matplotlib.pyplot as mpl
import matplotlib.patches as mpp
import matplotlib.collections as mpc
import matplotlib
import pairdata
import numpy as np
from auxclass import PowHist
from auxfun import enable_latex


def add_delay(self, binidx, val):
    if isinstance(self.bins[binidx], tuple):
        self.bins[binidx] = (1, val)
    else:
        self.bins[binidx] = (
        (self.bins[binidx][1] * self.bins[binidx][0] + val) / (self.bins[binidx] + 1), self.bins[binidx] + 1)


class DelaySpreadHist:
    def __init__(self, source: pairdata = None, binc: int = 10, brange: tuple = (1, 11)):
        self.source = source
        self.binc = binc
        self.binrange = brange

    def export(self, rxgrp: list = [-1], txgrp: list = [-1], xlims: tuple = (0, 10), ylims: tuple = (0, 10), figdims: tuple = (640, 480),
               mkpng: bool = False, matsav: bool = False, plot: bool = True):

        rms_hist = PowHist(binc=self.binc, rstart=self.binrange[0], rstop=self.binrange[1], addfun=add_delay)

        fig = mpl.figure()

        for i in self.source.txs.items():
            # TX is in valid group or group is ignored?
            if i[1].setid in txgrp or txgrp[0] == -1:
                # Which RXes are reached by TX
                for j in i[1].chans_to_pairs.items():
                    # Destination in a valid group or group is ignored?
                    if j[0].setid in rxgrp or rxgrp[0] == -1:
                        delay_vector = []
                        weights_vector = []
                        for k in j[1].paths.items():
                            delay_vector.append(k[1].delay)
                            weights_vector.append(k[1].pow)
                        mw = sum(weights_vector)
                        md = min(delay_vector)
                        weights_vector = [w/mw for w in weights_vector]
                        delay_vector = [d - md for d in delay_vector]

                        rms = np.sqrt((np.sum([w * (np.average(a=delay_vector, weights=weights_vector) - d)**2 for d, w
                                               in zip(delay_vector, weights_vector)])))
                        rms_hist.append(j[1].dist, rms)

        bars = []
        ticks = []
        eebars = []
        yebars = []
        xebars = []

        for i in rms_hist.bins.items():
            if not isinstance(i[1], list):
                bars.append(mpp.Rectangle((i[0][0], 0.0), i[0][1] - i[0][0], i[1] * 1e9))
            else:
                bars.append(mpp.Rectangle((i[0][0], 0.0), i[0][1] - i[0][0], np.mean(i[1]) * 1e9))

            ticks.append(i[0][0])

            xebars.append((i[0][0] + i[0][1]) / 2.0)
            yebars.append(np.mean(i[1]) * 1e9)
            eebars.append(np.var(np.asarray(i[1])*1e9))

        rms_bc = mpc.PatchCollection(bars, facecolor='r', alpha=0.5, edgecolors='k')

        ax = fig.add_subplot(111)
        ax.grid(linestyle='--')
        ax.set_xlim(self.binrange[0], self.binrange[1])
        ax.set_ylim(0, 10)
        mpl.xlabel('Distance [m]')
        mpl.ylabel('RMS delay spread [ns]')
        mpl.xticks(ticks)
        mpl.errorbar(xebars, yebars, eebars, linestyle='', linewidth=0)
        mpl.tight_layout()
        ax.add_collection(rms_bc)

        if plot:
            mpl.show()


if __name__ == '__main__':
    enable_latex(18)
    DS = pairdata.DataStorage(dbname='/home/nonroot/Downloads/Bus_geom.HHD.sqlite')
    DS.load_rxtx()
    DS.load_paths()
    DS.load_interactions()
    dspr = DelaySpreadHist(DS)
    dspr.export(matsav=True, plot=True, rxgrp=[5, 6])
    exit()
