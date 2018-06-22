import sqlite3
import os.path
import numpy as np

TX_EXTR = "SELECT tx_id, x, y, z FROM tx"
RX_EXTR = "SELECT rx_id, x, y, z FROM rx"
TX_PAIRS = "SELECT channel_id, rx_id FROM channel WHERE tx_id = {}"


class Node():
    def __init__(self, typ: str):
        self.chan_to_pairs = dict()
        self.node_id = 0
        self.coords = [0, 0, 0]
        self.rot = [0, 0, 0]

        if typ == 'TX':
            self.txpow = 0.0
        else:
            self.rxpow = 0.0

class chan():
    def __init__(self):
        self.paths = map()

class path():
    def __init__(self):
        self.pow = -np.Inf
        self.dist = 0.0
        self.interactions = list()


class interaction():
    def __init__(self):
        self.typ = 'TX'
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
            for j in self.dbcurs.execute(TX_PAIRS.format(i.node_id)):
                i.chan_to_pairs[j[1]] = j[0]



if __name__ == '__main__':
    DS = data_stor()
    DS.load_rxtx('/home/alexey/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/Mobile_TXRX/Class@60GHz.Mobile_TXRX.sqlite')
    exit()