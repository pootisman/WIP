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
import mysql.connector as msqlc
from mysql.connector.constants import ClientFlag
import os
import numpy as np
import concurrent.futures as cof
from multiprocessing import cpu_count


__author__ = 'Aleksei Ponomarenko-Timofeev'

TX_EXTR = 'SELECT tx_id, x, y, z, tx_set_id FROM tx;'
RX_EXTR = 'SELECT rx_id, x, y, z, rx_set_id FROM rx;'
TX_PAIRS = 'SELECT * FROM (SELECT channel_utd.channel_id, channel_utd.received_power, ' \
           'channel_utd.mean_time_of_arrival, channel_utd.delay_spread ' \
           'FROM channel_utd WHERE channel_utd.channel_id IN ' \
           '(SELECT channel_id FROM channel WHERE tx_id = {})) utd ' \
           'JOIN ' \
           '(SELECT channel.channel_id ,channel.rx_id FROM channel WHERE tx_id = {}) chan ' \
           'ON utd.channel_id = chan.channel_id;'
CHAN_PTH = 'SELECT path_utd_id, received_power, time_of_arrival, departure_phi, departure_theta, arrival_phi,' \
           ' arrival_theta, freespace_path_loss FROM path_utd WHERE path_id IN (SELECT path_id FROM' \
           ' path WHERE channel_id = {});'
INTERS = 'SELECT * FROM interaction_type;'
INTERS_SPEC_CHAN = 'SELECT x, y, z, interaction_type_id, path_id FROM interaction WHERE path_id IN' \
                   '(SELECT path_id FROM path WHERE channel_id = {});'
INTERS_SPEC = 'SELECT x,y,z,interaction_type_id FROM interaction WHERE path_id = {};'


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
        self.clusters = dict()

    def __repr__(self):
        return 'Radio channel'

    def __str__(self):
        return 'Sumpow = {}, delay = {}, ds = {}, Ncl = {}, {} -> {}'.format(self.pow, self.delay, self.ds,
                                                                            self.clusters.__len__(), self.src.node_id,
                                                                            self.dest.node_id)


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
        self.cluster = None
        self.near_field_failed = False

    def __repr__(self):
        return 'Propagation path'

    def __str__(self):
        return 'Inters = {}, pow = {} mW, len = {} m, {} [Az:{},El:{}]-->[Az:{},El:{}] {}'.format(self.interactions.__len__(), self.pow, self.len,
                                                                    self.chan.src.node_id, self.AoD, self.EoD, self.AoA, self.EoA, self.chan.dest.node_id)


class interaction():
    def __init__(self):
        self.typ = 'TX'
        self.coords = [0, 0, 0]
        self.path = None

    def __repr__(self):
        print('Interaction')

    def __str__(self):
        print('{}@{}'.format(self.typ, self.coords))


class data_stor():
    def __init__(self, conf: str = None):
        self.txs = dict()
        self.rxs = dict()
        self.dbname = None
        self.dbconn = None
        self.dbcurs = None
        self.possible_inters = dict()
        self.threshold_chan = -115
        self.threshold_path = -125
        if conf is not None:
            conff = open(conf)
            self.host = conff.readline().strip('\n')
            self.user = conff.readline().strip('\n')
            self.pasw = conff.readline().strip('\n')
            print('Connecting to {} as {}'.format(self.host, self.user))
            conff.close()
            self.nthreads = cpu_count()
            print('Preparing {} threads...'.format(self.nthreads))
            self.pool = cof.ThreadPoolExecutor(max_workers=self.nthreads)

    def __repr__(self):
        return 'SQL data storage'

    def __str__(self):
        if hasattr(self, 'host'):
            return 'Database {} at {}.'.format(self.dbname, self.host)
        else:
            return 'Database in file {}.'.format(self.dbname)

    def load_rxtx(self, dbname: str = None):
        print('Loading TX/RX nodes...', end='', flush=True)
        if self.dbconn is None:
            if not hasattr(self, 'host'):
                self.dbconn = sqlite3.connect(dbname)
                self.dbcurs = self.dbconn.cursor()
                self.dbname = dbname
            else:
                self.dbconn = msqlc.connect(host=self.host, user=self.user, password=self.pasw,
                                            client_flags=[ClientFlag.SSL], database=dbname)
                self.dbcurs = self.dbconn.cursor()
                self.dbcurs.execute('USE {};'.format(dbname))
                self.dbname = dbname

        # Read TX data
        self.dbcurs.execute(TX_EXTR)
        j = self.dbcurs.fetchall()
        for i in j:
            n = Node('TX')
            n.coords = np.asarray([i[1], i[2], i[3]])
            n.node_id = i[0]
            n.setid = i[4]
            self.txs[i[0]] = n

        # Read RX data
        self.dbcurs.execute(RX_EXTR)
        j = self.dbcurs.fetchall()
        for i in j:
            n = Node('RX')
            n.coords = np.asarray([i[1], i[2], i[3]])
            n.node_id = i[0]
            n.setid = i[4]
            self.rxs[i[0]] = n
        print('Success!', flush=True)
        if hasattr(self, 'host'):
            self.dbconn.close()
            self.dbconn = None

    def load_paths(self, npaths: int = 25):
        print('Loading paths...', end='', flush=True)
        self.npaths = npaths
        if self.dbconn is None:
            if os.path.isfile(self.dbname) and not hasattr(self, 'host'):
                print('Error: connect to DB and load txs/rxs first!')
                exit(1)
            else:
                self.dbconn = msqlc.connect(host=self.host, user=self.user, password=self.pasw,
                                            client_flags=[ClientFlag.SSL], database=self.dbname)
                self.dbcurs = self.dbconn.cursor()

        for i in self.txs.keys():
            self.dbcurs.execute(TX_PAIRS.format(i, i))
            k = self.dbcurs.fetchall()
            for j in k:
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
        print('Success!')

        if hasattr(self, 'host'):
            self.dbconn.close()
            self.dbconn = None

    def load_interactions(self, store: bool = True):
        print('Loading interactions...', end='', flush=True)

        if self.dbconn is None:
            if os.path.isfile(self.dbname) and not hasattr(self, 'host'):
                print('Error: connect to DB and load txs/rxs first!')
                exit(1)
            else:
                self.dbconn = msqlc.connect(host=self.host, user=self.user, password=self.pasw,
                                            client_flags=[ClientFlag.SSL], database=self.dbname)
                self.dbcurs = self.dbconn.cursor()

        for i in self.txs.items():
            for j in i[1].chans_to_pairs.items():
                self.dbcurs.execute(INTERS_SPEC_CHAN.format(j[1].chid))
                inters = self.dbcurs.fetchall()
                for k in j[1].paths.items():
                    for l in inters:
                        if l[4] == k[1].pathid:
                            if store:
                                intr = interaction()
                                intr.path = k[1]
                                intr.coords = np.asarray([l[0], l[1], l[2]])
                                intr.typ = l[3]
                            else:
                                intr = False
                            k[1].interactions.append(intr)
        print('Success!', flush=True)

        if hasattr(self, 'host'):
            self.dbconn.close()
            self.dbconn = None

    #def save_procd_to(self, suffix: str = 'procd'):



if __name__ == '__main__':
    DS = data_stor(conf='dbconf.txt')
    DS.load_rxtx(dbname='Human_crawl_TEST_sqlite')
    DS.load_paths()
    DS.load_interactions()
    from phys_path_procs import *
    check_data_NF(DS)
    exit()