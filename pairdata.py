# Copyright (C) Aleksei Ponomarenko-Timofeev
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import sqlite3
import os
import numpy as np

__author__ = 'Aleksei Ponomarenko-Timofeev'

TX_EXTR = "SELECT tx_id, x, y, z, tx_set_id FROM tx"
RX_EXTR = "SELECT rx_id, x, y, z, rx_set_id FROM rx"
TX_PAIRS = "SELECT * FROM (SELECT channel_utd.channel_id, channel_utd.received_power, " \
           "channel_utd.mean_time_of_arrival, channel_utd.delay_spread, channel_utd.channel_id " \
           "FROM channel_utd WHERE channel_utd.channel_id IN " \
           "(SELECT channel_id FROM channel WHERE tx_id = {})) utd " \
           "JOIN " \
           "(SELECT channel.channel_id ,channel.rx_id FROM channel WHERE tx_id = {}) chan " \
           "ON utd.channel_id = chan.channel_id"
CHAN_PTH = "SELECT path_utd_id, received_power, time_of_arrival, departure_phi, departure_theta, arrival_phi, arrival_theta, freespace_path_loss FROM path_utd WHERE path_id IN (SELECT path_id FROM" \
           " path WHERE channel_id = {})"
INTERS = "SELECT * FROM interaction_type"
INTERS_SPEC = "SELECT x,y,z,interaction_type_id FROM interaction WHERE path_id = {}"

class Node():
    def __init__(self, typ: str):
        self.chans_to_pairs = dict()
        self.node_id = 0
        self.coords = [0, 0, 0]
        self.rot = [0, 0, 0]
        self.setid = 0
        self.type = typ

        if typ == 'TX':
            self.txpow = 0.0
        else:
            self.rxpow = 0.0

    def chan_to(self, dest):
        for i in self.chans_to_pairs.keys():
            if self.chans_to_pairs[i].dest == dest:
                return self.chans_to_pairs[i]
        return None

class chan():
    def __init__(self, dest: Node = None, src: Node = None):
        self.paths = dict()
        self.dest = dest
        self.src = src
        self.pow = 0.0
        self.delay = 0.0
        self.ds = 0.0
        self.dist = 0.0
        self.chid = 0

class path():
    def __init__(self):
        self.pathid = 0
        self.pow = 0.0
        self.delay = 0.0
        self.len = 0.0
        self.interactions = list()
        self.AoA = 0.0
        self.EoA = 0.0
        self.AoD = 0.0
        self.EoD = 0.0
        self.FSPL = 0.0
        self.chan = None


class interaction():
    def __init__(self):
        self.typ = 'TX'
        self.coord = [0, 0, 0]
        self.path = None


class data_stor():
    def __init__(self):
        self.txs = dict()
        self.rxs = dict()
        self.dbname = None
        self.dbconn = None
        self.dbcurs = None
        self.possible_inters = dict()
        self.threshold_chan = -115
        self.threshold_path = -125

    def load_rxtx(self, dbname: str = None):
        if self.dbconn is None:
            if os.path.isfile(dbname):
                self.dbconn = sqlite3.connect(dbname)
                self.dbcurs = self.dbconn.cursor()
                self.dbname = dbname

        # Read TX data
        for i in self.dbcurs.execute(TX_EXTR):
            n = Node('TX')
            n.coords = np.asarray([i[1], i[2], i[3]])
            n.node_id = i[0]
            n.setid = i[4]
            self.txs[i[0]] = n

        # Read RX data
        for i in self.dbcurs.execute(RX_EXTR):
            n = Node('RX')
            n.coords = np.asarray([i[1], i[2], i[3]])
            n.node_id = i[0]
            n.setid = i[4]
            self.rxs[i[0]] = n

    def load_paths(self, npaths: int = 25):
        self.npaths = npaths

        if self.dbconn is None:
            print('Error: connect to DB and load txs/rxs first!')
            exit(1)

        for i in self.txs.keys():
            for j in self.dbcurs.execute(TX_PAIRS.format(i, i)):
                dst = self.rxs[j[-1]]
                self.txs[i].chans_to_pairs[dst] = chan(dst, self.txs[i])
                self.rxs[j[-1]].chans_to_pairs[self.txs[i]] = self.txs[i].chans_to_pairs[dst]
                self.txs[i].chans_to_pairs[dst].pow = j[1] * 1e3
                self.txs[i].chans_to_pairs[dst].delay = j[2]
                self.txs[i].chans_to_pairs[dst].ds = j[3]
                self.txs[i].chans_to_pairs[dst].dist = np.linalg.norm(self.txs[i].coords - self.rxs[j[-1]].coords)
                self.txs[i].chans_to_pairs[dst].chid = j[0]

            for j in self.txs[i].chans_to_pairs.keys():
                self.dbcurs.execute(CHAN_PTH.format(self.txs[i].chans_to_pairs[j].chid))
                d = self.dbcurs.fetchall()
                d = sorted(d, key=lambda t: t[1], reverse=True)
                for k in d[0:(self.npaths if d.__len__() >= self.npaths else d.__len__())]:
                    self.txs[i].chans_to_pairs[j].paths[k[0]] = path()
                    self.txs[i].chans_to_pairs[j].paths[k[0]].pathid = k[0]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].chan = self.txs[i].chans_to_pairs[j]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].pow = k[1] * 1e3
                    self.txs[i].chans_to_pairs[j].paths[k[0]].FSPL = k[7]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].delay = k[2]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].AoD = k[3]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].EoD = k[4]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].AoA = k[5]
                    self.txs[i].chans_to_pairs[j].paths[k[0]].EoA = k[6]

    def load_interactions(self, store: bool = True):
        if self.dbconn is None:
            print('Error: connect to DB and load txs/rxs first!')
            exit(1)

        for i in self.dbcurs.execute(INTERS):
            self.possible_inters[i[0]] = i[1]

        for i in self.txs.items():
            for j in i[1].chans_to_pairs.items():
                for k in j[1].paths.items():
                    for l in self.dbcurs.execute(INTERS_SPEC.format(k[0])):
                        if store:
                            intr = interaction()
                            intr.path = k[1]
                            intr.coord = [l[0], l[1], l[2]]
                            intr.typ = self.possible_inters[l[3]]
                        else:
                            intr = False
                        k[1].interactions.append(intr)




if __name__ == '__main__':
    DS = data_stor()
    DS.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TEST_60_MKE_15/Class@60GHz.TEST_60_MKE_15.sqlite')
    DS.load_paths()
    DS.load_interactions()
    exit()