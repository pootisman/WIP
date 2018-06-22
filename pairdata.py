import sqlite3
import os.path
import numpy as np

TX_EXTR = "SELECT tx_id, x, y, z FROM tx"
RX_EXTR = "SELECT rx_id, x, y, z FROM rx"


class Node():
    def __init__(self, typ: str):
        self.paths_to_pairs = None
        self.pairs = None
        self.node_id = None
        self.coords = [0, 0, 0]
        self.rot = [0, 0, 0]

        if typ == 'TX':
            self.txpow = 0.0
        else:
            self.rxpow = 0.0


class path():
    def __init__(self):
        self.pow = -np.Inf
        self.dist = 0.0
        self.interactions = list()


class interaction():
    def __init__(self):
        self.typ = None
        self.coord = list()


class data_stor():
    def __init__(self):
        self.txs = list()
        self.rxs = list()
        self.paths = list()
        self.dbname = None
        self.dbconn = None
        self.dbcurs = None


    def load_rxtx(self, dbname: str = None):
        if self.dbconn is None:
            if os.path.isfile(dbname):
                self.dbconn = sqlite3.connect(dbname)
                self.dbcurs = self.dbconn.cursor()
                self.dbname = dbname

        # Read TX data
        for i in self.dbcurs.execute(TX_EXTR):
            n = Node('TX')
            n.coords = [i[1], i[2], i[3]]
            n.node_id = i[0]
            self.txs.append(n)

        # Read RX data
        for i in self.dbcurs.execute(RX_EXTR):
            n = Node('RX')
            n.coords = [i[1], i[2], i[3]]
            n.node_id = i[0]
            self.rxs.append(n)


    def load_path(self, dbname: str = None):
        if self.dbconn is None:
            if os.path.isfile(dbname):
                self.dbconn = sqlite3.connect(dbname)
                self.dbcurs = self.dbconn.cursor()
                self.dbname = dbname

        for i in self.txs:
            i.


if __name__ == '__main__':
    DS = data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/Mobile_TXRX/Class@60GHz.Mobile_TXRX.sqlite')
    exit()