import matplotlib.pyplot as mpl
import matplotlib.patches as mpp
import matplotlib.collections as mpc
import pairdata
import numpy as np
from auxclass import PowHist


class DelaySpreadHist:
    def __init__(self, source: pairdata = None, binc: int = 10, brange: tuple = (0, 10)):
        self.source = source
        self.binc = binc
        self.binrange = brange

    def export(self, rxgrp: list = [-1], txgrp: list = [-1], xlims: tuple = (0, 10), ylims: tuple = (0, 10), figdims: tuple = (640, 480),
               mkpng: bool = False, matsav: bool = False, plot: bool = True):

        rms_hist = PowHist(binc=self.binc, rstart=self.binrange[0], rstop=self.binrange[1])
        simple_hist = PowHist(binc=self.binc, rstart=self.binrange[0], rstop=self.binrange[1])

        fig = mpl.figure(figsize=figdims)

        for i in self.source.txs.items():
            # TX is in valid group or group is ignored?
            if i[1].setid in txgrp or txgrp[0] == -1:
                # Which RXes are reached by TX
                for j in i[1].chans_to_pairs.items():
                    # Destination in a valid group or group is ignored?
                    if j[0].setid in rxgrp or rxgrp[0] == -1:
                        delay_vector = []
                        for k in j[1].paths.items():
                            delay_vector.append(k[1].delay)
                        simple = np.max(delay_vector) - np.min(delay_vector)
                        rms = np.sqrt(1/len(delay_vector) * (np.sum([d**2 for d in delay_vector])))
                        print('Simple: {} | RMS: {} | Dist: {}'.format(simple, rms, j[1].dist))
                        rms_hist.append(j[1].dist, rms)
                        simple_hist.append(j[1].dist, simple)

        bars = []
        ticks = []
        for i in rms_hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1]*1e9))
            ticks.append(i[0][0])

        rms_bc = mpc.PatchCollection(bars, facecolor='r', alpha=0.5)

        bars = []
        ticks = []
        for i in simple_hist.bins.items():
            bars.append(mpp.Rectangle((i[0][0], 0), i[0][1] - i[0][0], i[1]*1e9))
            ticks.append(i[0][0])

        simple_bc = mpc.PatchCollection(bars, facecolor='g', alpha=0.5)

        ax = mpl.gca()
        ax.grid()
        ax.set_xlim(rms_hist.floor, rms_hist.ceiling)
        ax.add_collection(rms_bc)
        ax.add_collection(simple_bc)

        if plot:
            mpl.tight_layout()
            mpl.show()


if __name__ == '__main__':
    DS = pairdata.DataStorage(conf='dbconf.txt', dbname='Human_sitting_legsback_Sitting_sqlite')
    DS.load_rxtx()
    DS.load_paths()
    DS.load_interactions()
    dspr = DelaySpreadHist(DS)
    dspr.export(matsav=True, rxgrp=[-1], plot=True)
    exit()