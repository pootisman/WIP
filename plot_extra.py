import matplotlib.pyplot as mpl
from pairdata import DataStorage
from auxfun import l2db, enable_latex

PWRS = 0x01
PHS = 0x02
AOA = 0x03
EOA = 0x04
AOD = 0x05
EOD = 0x06
DELAY = 0x07

VALID = {PWRS: 'Power [dBm]', PHS: 'Phase [deg]', AOA: 'Arrival azimuth [deg]', EOA: 'Arrival elevation [deg]', AOD: 'Departure azimuth [deg]', EOD: 'Departure elevation [deg]', DELAY: 'Delay [ns]'}

class ValPlotter:
    def __init__(self, src):
        self.source = src

    def export(self, data_req: int = DELAY, rxgrp: list = [-1], txgrp: list = [-1],
               rxrange: list = [-1], txrange: list = [-1], mkpng: bool = False,
               mksvg: bool = True, draw: bool = True, nptoplot: int = 20,
               reverse: bool = False):

        if data_req not in VALID.keys():
            raise KeyError

        if not reverse:
            r = self.source.rxs
            t = self.source.txs
        else:
            t = self.source.rxs
            r = self.source.txs

        rxs = r.items(grprange=rxgrp, noderange=rxrange)

        for i in t.items(grprange=txgrp, noderange=txrange):
            xbuf = []
            ybuf = []
            pnsv = []

            for j in rxs:
                if i[1].chan_to(j[1]) is not None:
                    pn = 0
                    for k in i[1].chan_to(j[1]).paths.items():
                        xbuf.append(j[1].node_id)
                        pnsv.append(pn % 8)
                        pn += 1
                        if data_req == PWRS:
                            ybuf.append(l2db(k[1].pow))
                        elif data_req == PHS:
                            ybuf.append(k[1].phase)
                        elif data_req == AOA:
                            ybuf.append(k[1].aoa)
                        elif data_req == EOA:
                            ybuf.append(k[1].eoa)
                        elif data_req == AOD:
                            ybuf.append(k[1].aod)
                        elif data_req == EOD:
                            ybuf.append(k[1].eod)
                        elif data_req == DELAY:
                            ybuf.append(k[1].delay * 1e9)

                        if pn > nptoplot:
                            break

            if mkpng or mksvg or draw:
                fig = mpl.figure()
                mpl.grid()
                mpl.scatter(x=xbuf, y=ybuf, c=pnsv, s=3, cmap='Dark2')
                mpl.ylabel(VALID[data_req])
                mpl.xlabel('RX index')
                mpl.tight_layout()
                mpl.xlim([min(xbuf), max(xbuf)])

                if mksvg:
                    mpl.savefig("{}_{}_to_{}.svg".format(VALID[data_req], i[1].node_id, rxgrp))

                if mkpng:
                    mpl.savefig("{}_{}_to_{}.png".format(VALID[data_req], i[1].node_id, rxgrp))

                if draw:
                    mpl.show()


if __name__ == '__main__':
    DS = DataStorage(conf='dbconf.txt', dbname='Human_crawl_TEST_sqlite')
    VP = ValPlotter(DS)
    enable_latex(18)
    gg = [2,3,5,6]
    q = False

    for g in gg:
        VP.export(DELAY, rxgrp=[g], draw=q, reverse=False)
        VP.export(PWRS, rxgrp=[g], draw=q, reverse=False)
        VP.export(AOA, rxgrp=[g], draw=q, reverse=False)
        VP.export(EOA, rxgrp=[g], draw=q, reverse=False)
        VP.export(AOD, rxgrp=[g], draw=q, reverse=False)
        VP.export(EOD, rxgrp=[g], draw=q, reverse=False)
        VP.export(PHS, rxgrp=[g], draw=q, reverse=False)
