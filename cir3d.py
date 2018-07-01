import numpy as np
import matplotlib.pyplot as mpl
import pairdata
import numba


class cirs():
    def __init__(self, source):
        self.source = source

        self.xlim = [-np.Inf, np.Inf]
        self.ylim = [-np.Inf, np.Inf]
        self.zlim = [-np.Inf, np.Inf]

        self.xdata = []
        self.ydata = []
        self.zdata = []

    def draw(self, txrange, rxrange, txgrp: int = -1, rxgrp: int = -1, print: bool = False, cmap: str = 'viridis'):
        @numba.jit
        def basint(X,Y,Z, xc, yc):
            Xlim = [np.min(X), np.max(X)]
            Ylim = [np.min(Y), np.max(Y)]

            stepX = (Xlim[1] - Xlim[0])/xc
            stepY = (Ylim[1] - Ylim[0])/yc

            Xo = np.arange(start=Xlim[0], step=stepX, stop=Xlim[1])
            Yo = np.arange(start=Ylim[0], step=stepY, stop=Ylim[1])
            Zo = np.tile(np.min(Z), [xc, yc])

            for i in range(Z.__len__()):
                j = int(np.round((X[i] - Xlim[0])/stepX))
                k = int(np.round((Y[i] - Ylim[0])/stepY))
                if j == xc:
                    j=j-1
                if k == yc:
                    k=k-1

                if Zo[j,k] < Z[i]:
                    Zo[j,k] = Z[i]

            return (Xo, Yo, Zo)

        for i in range(txrange):
            if self.source.txs[i].setid == txgrp or txgrp == -1:
                mpl.figure(i)
                nn = 0

                for j in range(rxrange):
                    if self.source.rxs[j].setid == rxgrp or rxgrp == -1:
                        nn += 1
                        for k in self.source.txs[i].chan_to_pairs[j].paths.keys():
                            self.xdata.append(j)
                            self.ydata.append(self.source.txs[i].chan_to_pairs[j].paths[k].delay * 1e9)
                            self.zdata.append(10.0 * np.log10(self.source.txs[i].chan_to_pairs[j].paths[k].pow))

                (X, Y, Z) = basint(self.xdata, self.ydata, self.zdata, nn, 250)
                [X, Y] = np.meshgrid(X, Y)
                mpl.contourf(np.transpose(X), np.transpose(Y), Z, 40)
                mpl.clim([np.min(self.zdata), np.max(self.zdata) + 0.1 * np.max(self.zdata)])
                mpl.colorbar().set_label('RX power, [dBm]')
                mpl.xlabel('RX Position')
                mpl.ylabel('Delay, [ns]')
                mpl.title('CIR@TX #{}'.format(i))
                if print:
                    mpl.savefig('CIR3D_tx{0:03d}.png'.format(i))

        if print is False:
            mpl.show()


if __name__ == "__main__":
    DS = pairdata.data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/Class@60GHz.TEST_60_MKE_15.sqlite')
    DS.load_path()
    cir = cirs(DS)
    cir.draw(1, 120, rxgrp=4, print=True)