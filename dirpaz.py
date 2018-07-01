import numpy as np
import matplotlib.pyplot as mpl
import pairdata


class rec_pat():
    def __init__(self, source):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def draw(self, txrange, rxrange, txgrp: int = -1, rxgrp: int = -1, print: bool = False):
        for i in range(txrange):
            if self.source.txs[i].setid == txgrp or txgrp == -1:
                rr = 0
                for j in range(1, rxrange):
                    if self.source.rxs[j].setid == rxgrp or rxgrp == -1:
                        mpl.figure(rr)
                        rr+=1
                        th = []
                        tt = []
                        r = []
                        for k in self.source.txs[i].chan_to_pairs[j].paths.keys():
                            th.append((self.source.txs[i].chan_to_pairs[j].paths[k].AoA, 10.0 * np.log10(self.source.txs[i].chan_to_pairs[j].paths[k].pow)))
                        th = sorted(th, key=lambda t: t[0])
                        for t in th:
                            r.append(t[1])
                            tt.append(t[0])
                        ax = mpl.subplot(111, projection='polar')
                        ax.plot(np.deg2rad(tt), r)
                        ax.set_rmax(-100)
                        ax.set_rmin(-200)
                        mpl.title('RX #{}'.format(j))
                        if print:
                            mpl.savefig('Recpat_tx{0:03d}->rx{1:03d}.png'.format(i,j))
        if print is False:
            mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/Class@60GHz.TEST_60_MKE_15.sqlite')
    DS.load_path()
    cir =  rec_pat(DS)
    cir.draw(1, 120, rxgrp=4, print=True)